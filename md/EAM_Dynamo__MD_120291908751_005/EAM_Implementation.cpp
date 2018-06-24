//
// CDDL HEADER START
//
// The contents of this file are subject to the terms of the Common Development
// and Distribution License Version 1.0 (the "License").
//
// You can obtain a copy of the license at
// http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
// specific language governing permissions and limitations under the License.
//
// When distributing Covered Code, include this CDDL HEADER in each file and
// include the License file in a prominent location with the name LICENSE.CDDL.
// If applicable, add the following below this CDDL HEADER, with the fields
// enclosed by brackets "[]" replaced with your own identifying information:
//
// Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
//
// CDDL HEADER END
//

//
// Copyright (c) 2013--2018, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Ryan S. Elliott
//    Stephen M. Whalen
//


#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#include "EAM_Implementation.hpp"
#include "KIM_ModelDriverHeaders.hpp"

#define IGNORE_RESULT(fn) if(fn){}


//==============================================================================
//
// Implementation of EAM_Implementation public member functions
//
//==============================================================================

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
EAM_Implementation::EAM_Implementation(
    KIM::ModelDriverCreate * const modelDriverCreate,
    KIM::LengthUnit const requestedLengthUnit,
    KIM::EnergyUnit const requestedEnergyUnit,
    KIM::ChargeUnit const requestedChargeUnit,
    KIM::TemperatureUnit const requestedTemperatureUnit,
    KIM::TimeUnit const requestedTimeUnit,
    int * const ier)
    : embeddingData_(0),
      densityData_(0),
      rPhiData_(0),
      publishDensityData_(0),
      publish_rPhiData_(0),
      embeddingCoeff_(0),
      densityCoeff_(0),
      rPhiCoeff_(0),
      cachedNumberOfParticles_(0),
      densityValue_(0),
      embeddingDerivativeValue_(0),
      embeddingSecondDerivativeValue_(0)
{
  // initialize comments to null strings and set pointers for comment fields
  for (int i = 0; i < MAX_PARAMETER_FILES; ++i)
  {
    comments_[i][0] = 0;
    comments_ptr_[i] = comments_[i];
  }

  // set particleNames to null string
  particleNames_[0] = 0;

  // initialize private parameters
  for (int i = 0; i < MAX_NUMBER_OF_SPECIES; ++i)
  {
    particleNumber_[i] = 0;
    particleMass_[i] = 0.0;
    latticeConstant_[i] = 0.0;
    latticeType_[i][0] = 0;
  }

  FILE* parameterFilePointers[MAX_PARAMETER_FILES];
  int numberParameterFiles;
  modelDriverCreate->GetNumberOfParameterFiles(
      &numberParameterFiles);
  *ier = OpenParameterFiles(modelDriverCreate, numberParameterFiles,
                            parameterFilePointers);
  if (*ier) return;

  eamFileType_ = DetermineParameterFileTypes(modelDriverCreate,
                                             parameterFilePointers,
                                             numberParameterFiles);
  if (eamFileType_ == Error)
  {
    *ier = true;
    return;
  }

  SetOfFuncflData funcflData;
  *ier = ProcessParameterFileHeaders(modelDriverCreate,
                                     eamFileType_, parameterFilePointers,
                                     numberParameterFiles, funcflData);
  if (*ier)
  {
    CloseParameterFiles(numberParameterFiles, parameterFilePointers);
    return;
  }

  AllocateParameterMemory();

  *ier = ProcessParameterFileData(modelDriverCreate,
                                  eamFileType_, parameterFilePointers,
                                  numberParameterFiles, funcflData);
  CloseParameterFiles(numberParameterFiles, parameterFilePointers);
  if (*ier) return;

  *ier = ConvertUnits(modelDriverCreate,
                      requestedLengthUnit,
                      requestedEnergyUnit,
                      requestedChargeUnit,
                      requestedTemperatureUnit,
                      requestedTimeUnit);
  if (*ier) return;

  *ier = SetRefreshMutableValues(modelDriverCreate);
  if (*ier) return;

  *ier = RegisterKIMModelSettings(modelDriverCreate);
  if (*ier) return;

  *ier = RegisterKIMParameters(modelDriverCreate, eamFileType_);
  if (*ier) return;

  *ier = RegisterKIMFunctions(modelDriverCreate);
  if (*ier) return;

  // everything is good
  *ier = false;
  return;
}

//******************************************************************************
EAM_Implementation::~EAM_Implementation()
{ // note: it is ok to delete a null pointer and we have ensured that
  // everything is initialized to null

  // Memory that was allocated in AllocateParameterMemory
  Deallocate2DArray(embeddingData_);
  Deallocate3DArray(densityData_);
  Deallocate3DArray(rPhiData_);
  Deallocate2DArray(embeddingCoeff_);
  Deallocate3DArray(densityCoeff_);
  Deallocate3DArray(rPhiCoeff_);
  Deallocate2DArray(publishDensityData_);
  Deallocate2DArray(publish_rPhiData_);

  // Memory that was allocate in SetComputeMutableValues
  delete [] densityValue_;
  delete [] embeddingDerivativeValue_;
  delete [] embeddingSecondDerivativeValue_;

  // Nullify pointers
  densityValue_ = 0;
  embeddingDerivativeValue_ = 0;
  embeddingSecondDerivativeValue_ = 0;
}

//******************************************************************************
#include "KIM_ModelRefreshLogMacros.hpp"
int EAM_Implementation::Refresh(KIM::ModelRefresh * const modelRefresh)
{
  int ier;

  // Check that new cutoff value in KIM API object is still valid
  if (cutoffParameter_ > (numberRPoints_ + 1)*deltaR_)
  { // NOTE: the above check should actually be
    //       (cutoff > (numberRPoints_-1.0)*deltaR_).  However, it appears to
    //       be a defacto standard to set cutoff to numberRPoints*deltaR_.
    //       So, we allow for this slight amount of extrapolation.

    ier = true;
    LOG_ERROR("Model has cutoff value outside of the pair function"
              " interpolation domain");
    return ier;
  }

  for (int i = 0; i < numberModelSpecies_ ; i++)
  {
    for(int j = i; j < numberModelSpecies_; j++)
    {
      int indxPhi = i*numberModelSpecies_ + j - (i*i + i)/2;
      for(int k = 0; k < numberRPoints_; k++)
      {
        rPhiData_[i][j][k] = rPhiData_[j][i][k] = publish_rPhiData_[indxPhi][k];
      }
    }

    for(int j = 0; j < numberModelSpecies_; j++)
    {
      int indxDensity = (eamFileType_ == FinnisSinclair) ?
          i*numberModelSpecies_ + j : i;
      for(int k = 0; k < numberRPoints_; k++)
      {
        densityData_[i][j][k] = publishDensityData_[indxDensity][k];
      }
    }
  }

  ier = SetRefreshMutableValues(modelRefresh);
  if (ier) return ier;

  // nothing else to do for this case

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int EAM_Implementation::Compute(
    KIM::ModelCompute const * const modelCompute,
    KIM::ModelComputeArguments const * const modelComputeArguments)
{
  int ier;

  // KIM API Model Input compute flags
  bool isComputeProcess_dEdr = false;
  bool isComputeProcess_d2Edr2 = false;
  //
  // KIM API Model Output compute flags
  bool isComputeEnergy = false;
  bool isComputeForces = false;
  bool isComputeParticleEnergy = false;
  bool isComputeVirial = false;
  bool isComputeParticleVirial = false;
  //
  // KIM API Model Input
  int const* particleSpeciesCodes = 0;
  int const* particleContributing = 0;
  VectorOfSizeDIM const* coordinates = 0;
  //
  // KIM API Model Output
  double* energy = 0;
  double* particleEnergy = 0;
  VectorOfSizeDIM* forces = 0;
  VectorOfSizeSix* virial = 0;
  VectorOfSizeSix* particleVirial = 0;
  ier = SetComputeMutableValues(modelComputeArguments,
                                isComputeProcess_dEdr,
                                isComputeProcess_d2Edr2, isComputeEnergy,
                                isComputeForces, isComputeParticleEnergy,
                                isComputeVirial, isComputeParticleVirial,
                                particleSpeciesCodes, particleContributing,
                                coordinates, energy, particleEnergy, forces,
                                virial, particleVirial);
  if (ier) return ier;

  // Skip this check for efficiency
  //
  //ier = CheckParticleSpecies(modelComputeArguments, particleSpeciesCodes);
  //if (ier) return ier;

#include "EAM_ImplementationComputeDispatch.cpp"
  return ier;
}

//******************************************************************************
int EAM_Implementation::ComputeArgumentsCreate(
    KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate) const
{
  int ier;

  ier = RegisterKIMComputeArgumentsSettings(modelComputeArgumentsCreate);
  if (ier) return ier;

  // nothing else to do for this case

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int EAM_Implementation::ComputeArgumentsDestroy(
    KIM::ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy)
    const
{
  int ier;

  // nothing to do for this case

  // everything is good
  ier = false;
  return ier;
}


//==============================================================================
//
// Implementation of EAM_Implementation private member functions
//
//==============================================================================

//******************************************************************************
void EAM_Implementation::AllocateParameterMemory()
{
  // allocate memory for data
  AllocateAndInitialize2DArray(embeddingData_, numberModelSpecies_,
                               numberRhoPoints_);
  AllocateAndInitialize3DArray(densityData_, numberModelSpecies_,
                               numberModelSpecies_, numberRPoints_);
  AllocateAndInitialize3DArray(rPhiData_, numberModelSpecies_,
                               numberModelSpecies_, numberRPoints_);
  // allocate memory for non-repeat data
  AllocateAndInitialize2DArray(
      publishDensityData_,
      numberModelSpecies_
      *((eamFileType_ == FinnisSinclair) ? numberModelSpecies_ : 1),
      numberRPoints_);
  AllocateAndInitialize2DArray(publish_rPhiData_, numberUniqueSpeciesPairs_,
                               numberRPoints_);

  // allocate memory for coefficients
  AllocateAndInitialize2DArray(embeddingCoeff_, numberModelSpecies_,
                               numberRhoPoints_ * NUMBER_SPLINE_COEFF);
  AllocateAndInitialize3DArray(densityCoeff_, numberModelSpecies_,
                               numberModelSpecies_,
                               numberRPoints_ * NUMBER_SPLINE_COEFF);
  AllocateAndInitialize3DArray(rPhiCoeff_, numberModelSpecies_,
                               numberModelSpecies_,
                               numberRPoints_ * NUMBER_SPLINE_COEFF);
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int EAM_Implementation::OpenParameterFiles(
    KIM::ModelDriverCreate * const modelDriverCreate,
    int const numberParameterFiles,
    FILE* parameterFilePointers[MAX_PARAMETER_FILES])
{
  int ier;

  if (numberParameterFiles > MAX_PARAMETER_FILES)
  {
    ier = true;
    LOG_ERROR("EAM Dynamo driver given too many parameter files");
  }

  for (int i = 0; i < numberParameterFiles; ++i)
  {
    std::string const * paramFileName;
    ier = modelDriverCreate->GetParameterFileName(
        i,
        &paramFileName);
    if (ier)
    {
      LOG_ERROR("Unable to get parameter file name");
      return ier;
    }
    parameterFilePointers[i] = fopen(paramFileName->c_str(), "r");
    if (parameterFilePointers[i] == 0)
    {
      char message[MAXLINE];
      sprintf(message,
              "EAM parameter file number %d cannot be opened",
              i);
      ier = true;
      LOG_ERROR(message);
      for (int j = i - 1; i <= 0; --i)
      {
        fclose(parameterFilePointers[j]);
      }
      return ier;
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
EAMFileType EAM_Implementation::DetermineParameterFileTypes(
    KIM::ModelDriverCreate * const modelDriverCreate,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
    int const numberParameterFiles)
{
  if ((numberParameterFiles > 1) &&
      (numberParameterFiles <= MAX_PARAMETER_FILES))
  { // should be a set of Funcfl files
    for (int i = 0; i < numberParameterFiles; ++i)
    {
      if (IsFuncflOrSetfl(parameterFilePointers[i]) != Funcfl)
      {
        char message[MAXLINE];
        sprintf(message, "EAM parameter file number %d is not a"
                " funcfl file", i);
        LOG_ERROR(message);
        return Error;
      }
    }

    // everything is good
    return Funcfl;
  }
  else if (numberParameterFiles == 1)
  {
    EAMFileType eamFileType = IsFuncflOrSetfl(parameterFilePointers[0]);

    if (eamFileType == Error)
    {
      LOG_ERROR("Unable to determine parameter file type in EAM Dynamo");
      return Error;
    }

    // distinguish between setfl and Finnis-Sinclair files
    if (eamFileType == Setfl)
    {
      eamFileType = IsSetflOrFinnisSinclair(modelDriverCreate,
                                            parameterFilePointers[0]);
    }

    return eamFileType;
  }
  else
  {
    char message[MAXLINE];
    sprintf(message, "Invalid number (%d) of parameter files in EAM Dynamo",
            numberParameterFiles);
    LOG_ERROR(message);
    return Error;
  }
}

//******************************************************************************
EAMFileType EAM_Implementation::IsFuncflOrSetfl(FILE* const fptr)
{
  int const numberOfLinesToRead = 8;
  // use 1-based counting for line numbers
  bool isInteger[numberOfLinesToRead + 1];
  bool isDouble[numberOfLinesToRead + 1];
  int intValue[numberOfLinesToRead + 1];

  // discard 1st line.  It is a comment in both funcfl and setfl
  char line[MAXLINE];
  char const* cer = fgets(line, MAXLINE, fptr);
  if (cer == 0) return Error;

  // loop over the lines 2--6
  for (int i = 2; i <= numberOfLinesToRead; ++i)
  { // get next line
    cer = fgets(line, MAXLINE, fptr);
    if (cer == 0) return Error;

    // get first token on line
    char const* const tok = strtok(line, " ,\t\n");
    if (tok == 0)
    { // nothing on the line
      isInteger[i] = false;
      isDouble[i] = false;

      // skip to next line
      continue;
    }

    char* endptr;
    intValue[i] = strtol(tok, &endptr, 10);
    if (*endptr == 0)  // entire string used up, thus, tok was an int
    {
      isInteger[i] = true;
      isDouble[i] = false;
    }
    else
    {
      IGNORE_RESULT(strtod(tok, &endptr));
      if (*endptr == 0)  // entire string used up, thus, tok was a double
      {
        isInteger[i] = false;
        isDouble[i] = true;
      }
      else
      {
        isInteger[i] = false;
        isDouble[i] = false;
      }
    }
  }

  // done with the file rewind it for later use
  rewind(fptr);

  bool const isFuncfl = (  // line 2 starts with "ielem" and <= 118
      (isInteger[2] && (intValue[2] <= 118)) &&
      // line 3 starts with "nrho"
      isInteger[3] &&
      // line 4 starts with embedding data
      isDouble[4] &&
      // line 5 is double data
      isDouble[5] &&
      // line 6 is double data
      isDouble[6] &&
      // line 7 is double data
      isDouble[7] &&
      // line 8 is double data
      isDouble[8]);

  if (isFuncfl) return Funcfl;

  bool const isSetfl = (  // line 4 starts with "ntypes"
      isInteger[4] &&
      // line 5 starts with "nrho"
      isInteger[5] &&
      // line 6 starts with "ielem(1)" and <= 118
      (isInteger[6] && (intValue[6] <= 118)) &&
      // line 7 starts with embedding data
      isDouble[7] &&
      // line 8 is double data
      isDouble[8]);

  if (isSetfl) return Setfl;

  return Error;
}

//******************************************************************************
EAMFileType EAM_Implementation::IsSetflOrFinnisSinclair(
    KIM::ModelDriverCreate * const modelDriverCreate, FILE* const fptr)
{ // We are free to assume the file format (of the header, at least)
  // is conforms to the Setfl format

  char line[MAXLINE];
  // discard first three lines.  They are comments
  for (int i = 0; i < 3; ++i)
  {
    char const* cer = fgets(line, MAXLINE, fptr);
    if (cer == 0) return Error;
  }

  // 4th line; read number of elements
  int Nelements;
  {
    char const* const cer = fgets(line, MAXLINE, fptr);
    if (cer == 0) return Error;
    int ier = sscanf(line, "%d", &Nelements);
    if (ier != 1) return Error;
  }

  // 5th line; read Nrho and Nr
  int Nrho;
  int Nr;
  {
    char const* const cer = fgets(line, MAXLINE, fptr);
    if (cer == 0) return Error;
    double dummy;
    int const ier = sscanf(line, "%d %lg %d", &Nrho, &dummy, &Nr);
    if (ier != 3) return Error;
  }

  // 6th line; discard
  {
    char const* cer = fgets(line, MAXLINE, fptr);
    if (cer == 0) return Error;
  }

  // Read element1 data set
  {
    double* const dummy = new double[(Nrho>Nr)? Nrho : Nr];
    int ier = GrabData(modelDriverCreate, fptr, Nrho, dummy);
    if (ier)
    {
      delete[] dummy;
      return Error;
    }

    ier = GrabData(modelDriverCreate, fptr, Nr, dummy);
    if (ier)
    {
      delete[] dummy;
      return Error;
    }

    delete[] dummy;
  }

  // Read next line (Setfl - element header; FinnisSinclair - rho(r) values)
  bool isSetfl;
  {
    char const* cer = fgets(line, MAXLINE, fptr);
    if (cer == 0) return Error;

    // get first token on line
    char const* const tok = strtok(line, " ,\t\n");
    if (tok == 0)
    { // nothing on the line
      return Error;
    }

    char* endptr;
    IGNORE_RESULT(strtol(tok, &endptr, 10));
    if (*endptr == 0)  // entire string used up, thus, tok was an int
    {
      isSetfl = true;
    }
    else
    {
      IGNORE_RESULT(strtod(tok, &endptr));
      if (*endptr == 0)  // entire string used up, thus, tok was a double
      {
        isSetfl = false;
      }
      else
      {
        return Error;
      }
    }
  }

  // done with the file rewind it for later use
  rewind(fptr);

  if (isSetfl)
    return Setfl;
  else
    return FinnisSinclair;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int EAM_Implementation::ProcessParameterFileHeaders(
    KIM::ModelDriverCreate * const modelDriverCreate,
    EAMFileType const eamFileType,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
    int const numberParameterFiles,
    SetOfFuncflData& funcflData)
{
  int ier;

  switch (eamFileType)
  {
    case FinnisSinclair:
    { // FinnisSinclair has same header information as setfl
      ier = ReadSetflHeader(modelDriverCreate, parameterFilePointers[0]);
      if (ier)
      {
        LOG_ERROR("Could not read FinnisSinclair parameter file header.");
        return ier;
      }
      break;
    }
    case Setfl:
    {
      ier = ReadSetflHeader(modelDriverCreate, parameterFilePointers[0]);
      if (ier)
      {
        LOG_ERROR("Could not read Setfl parameter file header");
        return ier;
      }
      break;
    }
    case Funcfl:
    {
      // set number of species to be the number of parameter files
      numberModelSpecies_ = numberParameterFiles;
      numberUniqueSpeciesPairs_
          = ((numberModelSpecies_+1)*numberModelSpecies_)/2;

      // initialize grid values
      deltaRho_ = 0.0;
      deltaR_ = 0.0;
      cutoffParameter_ = 0.0;
      double rhoMax = 0.0;
      double rMax = 0.0;

      for (int i = 0; i < numberParameterFiles; ++i)
      {
        ier = ReadFuncflHeader(modelDriverCreate,
                               parameterFilePointers[i],
                               i,
                               funcflData.numberRhoPoints[i],
                               funcflData.deltaRho[i],
                               funcflData.numberRPoints[i],
                               funcflData.deltaR[i],
                               funcflData.cutoff[i]);
        if (ier)
        {
          LOG_ERROR("Could not read Funcfl parameter file header");
          return ier;
        }

        // update grid values (use maximums)
        deltaRho_ = std::max(deltaRho_, funcflData.deltaRho[i]);
        deltaR_ = std::max(deltaR_, funcflData.deltaR[i]);
        cutoffParameter_ = std::max(cutoffParameter_, funcflData.cutoff[i]);
        rhoMax = std::max(rhoMax, ((funcflData.numberRhoPoints[i] - 1) *
                                   funcflData.deltaRho[i]));
        rMax = std::max(rMax, ((funcflData.numberRPoints[i] - 1) *
                               funcflData.deltaR[i]));
      }

      // determine number of rho and r points in grid.
      // add 1 to account for point at zero
      numberRhoPoints_ = int(rhoMax / deltaRho_ + 0.5) + 1;
      numberRPoints_ = int(rMax / deltaR_ + 0.5) + 1;

      // set particleNames_
      ier = SetParticleNamesForFuncflModels(modelDriverCreate);
      if (ier)
      {
        LOG_ERROR("Could not set particle names");
        return ier;
      }
      break;
    }
    default:  // should never get here
    {
      ier = true;
      LOG_ERROR("Invalid valid parameter files passed to EAM Dynamo");
      return ier;
      break;
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int EAM_Implementation::ReadSetflHeader(
    KIM::ModelDriverCreate * const modelDriverCreate, FILE* const fptr)
{
  int ier;
  char const* cer;
  char line[MAXLINE];

  // read lines 1, 2, and 3 (comment lines)
  for (int i = 0; i < NUMBER_SETFL_COMMENT_LINES; ++i)
  {
    cer = fgets(&comments_[i][0], MAXLINE, fptr);
    if (cer == 0)
    {
      ier = true;
      LOG_ERROR("Error reading comment lines in Setfl parameter file");
      return ier;
    }
    int const cmntlength = strlen(&comments_[i][0]);
    if (comments_[i][cmntlength-1] == '\n') comments_[i][cmntlength-1] = 0;
  }

  // read 4th line (Nelements Element1 Element2 ... ElementN)
  cer = fgets(particleNames_, MAXLINE, fptr);
  int const nameslength = strlen(particleNames_);
  if (particleNames_[nameslength-1] == '\n') particleNames_[nameslength-1] = 0;
  // parse number of particle species
  int number_of_species;
  ier = sscanf(particleNames_, "%d", &number_of_species);
  if ((cer == 0) || (ier != 1))
  {
    ier = true;
    LOG_ERROR("Error reading fourth line of Setfl parameter file");
    return ier;
  }

  numberModelSpecies_ = number_of_species;
  numberUniqueSpeciesPairs_ = ((numberModelSpecies_+1)*numberModelSpecies_)/2;

  // parse the remainder of 4th line for particle species names
  char* const tmpnames = new char[strlen(particleNames_)+1];
  strcpy(tmpnames, particleNames_);
  char** const elems = new char*[numberModelSpecies_];
  char* tmpstring = strtok(tmpnames, " ,\t");  // ignore first token
  int counter = 0;
  while (tmpstring != 0)
  {
    tmpstring = strtok(0, " ,\t\n");
    elems[counter] = tmpstring;
    ++counter;
    if (counter >= numberModelSpecies_)
    {
      break;
    }
  }

  // register species names in the KIM API object and give them indices
  // that we'll use internally
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    modelDriverCreate->SetSpeciesCode(std::string(elems[i]), i);
  }
  delete [] elems;
  delete [] tmpnames;

  // read 5th line (Nrho, deltaRho, Nr, deltaR, cutoff)
  cer = fgets(line, MAXLINE, fptr);
  ier = sscanf(line, "%d %lg %d %lg %lg", &numberRhoPoints_, &deltaRho_,
               &numberRPoints_, &deltaR_, &cutoffParameter_);
  if ((cer == 0) || (ier != 5))
  {
    ier = true;
    LOG_ERROR("Error reading fifth line of Setfl parameter file");
    return ier;
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int EAM_Implementation::ReadFuncflHeader(
    KIM::ModelDriverCreate * const modelDriverCreate,
    FILE* const fptr,
    int const fileIndex,
    int& numberRhoPoints,
    double& deltaRho,
    int& numberRPoints,
    double& deltaR,
    double& cutoffParameter)
{
  int ier;
  char const* cer;
  char line[MAXLINE];

  // read 1st line (comment line)
  cer = fgets(&comments_[fileIndex][0], MAXLINE, fptr);
  if (cer == 0)
  {
    ier = true;
    LOG_ERROR("Error reading first line (the comment line) of Funcfl "
              "parameter file");
    return ier;
  }
  int const cmntlength = strlen(&comments_[fileIndex][0]);
  if (comments_[fileIndex][cmntlength-1] == '\n')
    comments_[fileIndex][cmntlength-1] = 0;

  // read 2nd line (particle number, mass, lattice constant, and lattice type)
  cer = fgets(line, MAXLINE, fptr);
  ier = sscanf(line, "%d %lg %lg %s", &(particleNumber_[fileIndex]),
               &(particleMass_[fileIndex]), &(latticeConstant_[fileIndex]),
               latticeType_[fileIndex]);
  if ((cer == 0) || (ier != 4))
  {
    ier = true;
    LOG_ERROR("Error reading second line of Funcfl parameter file");
    return ier;
  }

  // read 3rd line (Nrho, deltaRho, Nr, deltaR, cutoff)
  cer = fgets(line, MAXLINE, fptr);
  ier = sscanf(line, "%d %lg %d %lg %lg", &numberRhoPoints, &deltaRho,
               &numberRPoints, &deltaR, &cutoffParameter);
  if ((cer == 0) || (ier != 5))
  {
    ier = true;
    LOG_ERROR("Error reading third line of Funcfl parameter file");
    return ier;
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int EAM_Implementation::SetParticleNamesForFuncflModels(
    KIM::ModelDriverCreate * const modelDriverCreate)
{
  int ier;

  // get and correctly order the particle names
  const char** const particleNames = new const char*[numberModelSpecies_];
  KIM::SpeciesName tmp_species;
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    // Our file indexing corresponds to our species indexing, so we can
    // pull up the atomic number using particleNumber_[i] and then ask the KIM
    // API to give us the chemical symbol using that and store it in
    // particleNames
    ier = KIM::SPECIES_NAME::GetSpeciesName(particleNumber_[i], &tmp_species);
    if (ier)
    {
      LOG_ERROR("Error retrieving species names from atomic numbers read from"
                "parameter files");
      delete [] particleNames;
      return ier;
    }
    particleNames[i] = tmp_species.String().c_str();
  }

  // write particleNames_ string and register species names in the
  // KIM API object and give them indices that we'll use internally
  sprintf(particleNames_, "%d ", numberModelSpecies_);
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    strcat(particleNames_, particleNames[i]);
    strcat(particleNames_, " ");

    modelDriverCreate->SetSpeciesCode(std::string(particleNames[i]), i);
  }
  int const nmlength = strlen(particleNames_);
  particleNames_[nmlength - 1] = 0;
  delete [] particleNames;

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int EAM_Implementation::ProcessParameterFileData(
    KIM::ModelDriverCreate * const modelDriverCreate,
    EAMFileType const eamFileType,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
    int const numberParameterFiles,
    SetOfFuncflData& funcflData)
{
  int ier;

  // read data file based on indicated type
  switch (eamFileType)
  {
    case FinnisSinclair:
    {
      ier = ReadFinnisSinclairData(modelDriverCreate, parameterFilePointers[0]);
      if (ier)
      {
        LOG_ERROR("Error reading tabulated data from Finnis-Sinclair"
                  "parameter file");
        return ier;
      }
      break;
    }
    case Setfl:
    {
      ier = ReadSetflData(modelDriverCreate, parameterFilePointers[0]);
      if (ier)
      {
        LOG_ERROR("Error reading tabulated data from Setfl parameter file");
        return ier;
      }
      break;
    }
    case Funcfl:
    {
      for (int i = 0; i < numberParameterFiles; ++i)
      { // allocate memory for Funcfl data
        funcflData.embeddingData[i] = new double[funcflData.numberRhoPoints[i]];
        funcflData.densityData[i] = new double[funcflData.numberRPoints[i]];
        funcflData.ZData[i] = new double[funcflData.numberRPoints[i]];

        ier = ReadFuncflData(modelDriverCreate, parameterFilePointers[i], i,
                             funcflData);
        if (ier)
        {
          LOG_ERROR("Error reading tabulated data from Funcfl parameter file");
          for (int j = 0; j <= i; ++j)
          {
            delete [] funcflData.embeddingData[i];
            delete [] funcflData.densityData[i];
            delete [] funcflData.ZData[i];
          }
          return ier;
        }
      }

      ReinterpolateAndMix(funcflData);
      for (int i = 0; i < numberParameterFiles; ++i)
      {
        delete [] funcflData.embeddingData[i];
        delete [] funcflData.densityData[i];
        delete [] funcflData.ZData[i];
      }

      break;
    }
    default:  // should never get here
    {
      ier = true;
      LOG_ERROR("Invalid valid parameter files passed to EAM Dynamo");
      return ier;
    }
    break;
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int EAM_Implementation::ReadSetflData(
    KIM::ModelDriverCreate * const modelDriverCreate, FILE* const fptr)
{
  int ier;
  char const* cer;
  char line[MAXLINE];

  // loop over each atom type in the data file
  for (int i = 0; i < numberModelSpecies_; ++i)
  { // read header line (partcle number, mass, lattice constant, lattice type)
    cer = fgets(line, MAXLINE, fptr);
    ier = sscanf(line, "%d %lg %lg %s", &(particleNumber_[i]),
                 &(particleMass_[i]), &(latticeConstant_[i]), latticeType_[i]);
    if ((cer == 0) || (ier != 4))
    {
      ier = true;
      LOG_ERROR("Error reading lines of setfl file");
      return ier;
    }

    // read "embed_dat"
    ier = GrabData(modelDriverCreate,
                   fptr, numberRhoPoints_, embeddingData_[i]);
    if (ier)
    {
      LOG_ERROR("Error reading embeddingData lines of setfl file");
      return ier;
    }

    // read "densityData"
    ier = GrabData(modelDriverCreate, fptr, numberRPoints_, densityData_[i][0]);
    if (ier)
    {
      LOG_ERROR("Error reading densityData lines of setfl file");
      return ier;
    }
    // fill in remaining columns
    for (int j = 1; j < numberModelSpecies_; ++j)
    {
      for (int k = 0; k < numberRPoints_; ++k)
      {
        densityData_[i][j][k] = densityData_[i][0][k];
      }
    }
  }

  // read "rPhiData"
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      ier = GrabData(modelDriverCreate, fptr, numberRPoints_, rPhiData_[i][j]);
      if (ier)
      {
        LOG_ERROR("Error reading rPhiData lines of setfl file");
        return ier;
      }
    }
  }

  // filling in upper-triangular part of rPhiData
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    for (int j = i + 1; j < numberModelSpecies_; ++j)
    {
      for (int k = 0; k < numberRPoints_; ++k)
      {
        rPhiData_[i][j][k] = rPhiData_[j][i][k];
      }
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int EAM_Implementation::ReadFinnisSinclairData(
    KIM::ModelDriverCreate * const modelDriverCreate, FILE* const fptr)
{
  int ier;
  char const* cer;
  char line[MAXLINE];

  // loop over each atom species in the data file
  for (int i = 0; i < numberModelSpecies_; ++i)
  { // read header line (partcle number, mass, lattice constant, lattice type)
    cer = fgets(line, MAXLINE, fptr);
    ier = sscanf(line, "%d %lg %lg %s", &(particleNumber_[i]),
                 &(particleMass_[i]), &(latticeConstant_[i]), latticeType_[i]);
    if ((cer == 0) || (ier != 4))
    {
      ier = true;
      LOG_ERROR("Error reading lines of setfl file");
      return ier;
    }

    // read "embed_dat"
    ier = GrabData(modelDriverCreate,
                   fptr, numberRhoPoints_, embeddingData_[i]);
    if (ier)
    {
      LOG_ERROR("Error reading embeddingData lines of setfl file");
      return ier;
    }

    // read "densityData"
    for (int j = 0; j < numberModelSpecies_; ++j)
    {
      ier = GrabData(modelDriverCreate,
                     fptr, numberRPoints_, densityData_[i][j]);
      if (ier)
      {
        LOG_ERROR("Error reading densityData lines of setfl file");
        return ier;
      }
    }
  }

  // read "rPhiData"
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    for (int j = 0; j <= i; ++j)
    {
      ier = GrabData(modelDriverCreate, fptr, numberRPoints_, rPhiData_[i][j]);
      if (ier)
      {
        LOG_ERROR("Error reading rPhiData lines of setfl file");
        return ier;
      }
    }
  }

  // filling in upper-triangular part of rPhiData
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    for (int j = i + 1; j < numberModelSpecies_; ++j)
    {
      for (int k = 0; k < numberRPoints_; ++k)
      {
        rPhiData_[i][j][k] = rPhiData_[j][i][k];
      }
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int EAM_Implementation::ReadFuncflData(
    KIM::ModelDriverCreate * const modelDriverCreate,
    FILE* const fptr,
    int const fileIndex,
    SetOfFuncflData& funcflData)
{
  int ier;

  // read "embed_dat"
  ier = GrabData(modelDriverCreate, fptr, funcflData.numberRhoPoints[fileIndex],
                 funcflData.embeddingData[fileIndex]);
  if (ier)
  {
    LOG_ERROR("Error reading embeddingData lines of funcfl file");
    return ier;
  }

  // read "Z_dat"
  ier = GrabData(modelDriverCreate, fptr, funcflData.numberRPoints[fileIndex],
                 funcflData.ZData[fileIndex]);
  if (ier)
  {
    LOG_ERROR("Error reading Z_dat lines of funcfl file");
    return ier;
  }

  // read "densityData"
  ier = GrabData(modelDriverCreate, fptr, funcflData.numberRPoints[fileIndex],
                 funcflData.densityData[fileIndex]);
  if (ier)
  {
    LOG_ERROR("Error reading densityData lines of funcfl file");
    return ier;
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int EAM_Implementation::GrabData(
    KIM::ModelDriverCreate * const modelDriverCreate,
    FILE* const fptr,
    int const n,
    double* const list)
{ // This function originally obtained under CDDL from Steve Plimpton
  int ier;
  char const* cer;
  char const* ptr;
  char line[MAXLINE];

  int i = 0;
  while (i < n)
  {
    cer = fgets(line, MAXLINE, fptr);
    if (cer == 0)
    {
      ier = true;
      LOG_ERROR("Error reading data from file");
      return ier;
    }

    ptr = strtok(line, " \t\n\r\f");
    list[i] = atof(ptr);
    ++i;
    while ((ptr = strtok(0, " \t\n\r\f")))
    {
      list[i++] = atof(ptr);
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
void EAM_Implementation::ReinterpolateAndMix(
    SetOfFuncflData const& funcflData)
{
  // conversion constants
  double const Hartree = 27.2;  // ev
  double const Bohr = 0.529;    // Angstroms

  // reinterpolate
  double const oneByDrho = ONE / deltaRho_;
  double const oneByDr = ONE / deltaR_;
  if (numberModelSpecies_ > 1)
  {
    for (int i = 0; i < numberModelSpecies_; ++i)
    {
      double* const embeddingCoeff
          = new double[funcflData.numberRhoPoints[i] * NUMBER_SPLINE_COEFF];
      double* const densityCoeff = new double[funcflData.numberRPoints[i] *
                                              NUMBER_SPLINE_COEFF];
      double* const rPhiCoeff = new double[funcflData.numberRPoints[i] *
                                           NUMBER_SPLINE_COEFF];

      SplineInterpolate(funcflData.embeddingData[i],
                        funcflData.deltaRho[i],
                        funcflData.numberRhoPoints[i], embeddingCoeff);
      SplineInterpolate(funcflData.densityData[i], funcflData.deltaR[i],
                        funcflData.numberRPoints[i], densityCoeff);
      SplineInterpolate(funcflData.ZData[i], funcflData.deltaR[i],
                        funcflData.numberRPoints[i], rPhiCoeff);

      for (int j = 0; j < numberRhoPoints_; ++j)
      {
        double densityOffset;
        int densityIndex;
        // compute densityOffset and densityIndex
        double const densityValue = j * deltaRho_;
        GET_DELTAX_AND_INDEX(densityValue, oneByDrho, numberRhoPoints_,
                             densityOffset, densityIndex);
        // interpolate value of embeddingData_[i][j]
        INTERPOLATE_F(embeddingCoeff, densityOffset, densityIndex,
                      embeddingData_[i][j]);
      }

      for (int j = 0; j < numberRPoints_; ++j)
      {
        double rOffset;
        int rIndex;
        double const r = j * deltaR_;
        // compute rOffset and rIndex
        GET_DELTAX_AND_INDEX(r, oneByDr, numberRPoints_, rOffset, rIndex);
        // interpolate value of densityData_[i][j]
        INTERPOLATE_F(densityCoeff, rOffset, rIndex, densityData_[i][0][j]);
        for (int k = 1; k < numberModelSpecies_; ++k)
        {
          densityData_[i][k][j] = densityData_[i][0][j];
        }
        // interpolate value of rPhiData_[i][i][j]
        INTERPOLATE_F(rPhiCoeff, rOffset, rIndex, rPhiData_[i][i][j]);
      }

      delete [] embeddingCoeff;
      delete [] densityCoeff;
      delete [] rPhiCoeff;
    }

    // convert "Z_dat" to r*phi, mix, and store in rPhiData_
    for (int i = 0; i < numberModelSpecies_; ++i)
    {
      for (int j = numberModelSpecies_ - 1; j > i; --j)
      {
        for (int k = 0; k < numberRPoints_; ++k)
          rPhiData_[j][i][k] = rPhiData_[i][j][k]
              = (rPhiData_[i][i][k] * rPhiData_[j][j][k]) * Hartree * Bohr;
      }
      for (int k = 0; k < numberRPoints_; ++k)
        rPhiData_[i][i][k] = (rPhiData_[i][i][k] * rPhiData_[i][i][k])
            * Hartree * Bohr;
    }
  }
  else
  { // if numberModelSpecies_ == 1, don't reinterpolate
    for (int i = 0; i < numberRhoPoints_; ++i)
      embeddingData_[0][i] = funcflData.embeddingData[0][i];
    for (int i = 0; i < numberRPoints_; ++i)
    {
      densityData_[0][0][i] = funcflData.densityData[0][i];
      rPhiData_[0][0][i] = funcflData.ZData[0][i] * funcflData.ZData[0][i]
          * Hartree * Bohr;
    }
  }
}

//******************************************************************************
void EAM_Implementation::CloseParameterFiles(
    int const numberParameterFiles,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES])
{
  for (int i = 0; i < numberParameterFiles; ++i)
    fclose(parameterFilePointers[i]);
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int EAM_Implementation::ConvertUnits(
    KIM::ModelDriverCreate * const modelDriverCreate,
    KIM::LengthUnit const requestedLengthUnit,
    KIM::EnergyUnit const requestedEnergyUnit,
    KIM::ChargeUnit const requestedChargeUnit,
    KIM::TemperatureUnit const requestedTemperatureUnit,
    KIM::TimeUnit const requestedTimeUnit)
{
  int ier;

  // define default base units
  KIM::LengthUnit fromLength = KIM::LENGTH_UNIT::A;
  KIM::EnergyUnit fromEnergy = KIM::ENERGY_UNIT::eV;
  KIM::ChargeUnit fromCharge = KIM::CHARGE_UNIT::e;
  KIM::TemperatureUnit fromTemperature = KIM::TEMPERATURE_UNIT::K;
  KIM::TimeUnit fromTime = KIM::TIME_UNIT::ps;

  // changing units of particle mass and lattice constant
  double convertMass = 1.0;
  ier = modelDriverCreate->ConvertUnit(
      fromLength, fromEnergy, fromCharge, fromTemperature, fromTime,
      requestedLengthUnit, requestedEnergyUnit, requestedChargeUnit,
      requestedTemperatureUnit, requestedTimeUnit,
      -2.0, 1.0, 0.0, 0.0, 2.0,
      &convertMass);
  if (ier)
  {
    LOG_ERROR("Unable to convert mass units");
    return ier;
  }
  double convertLength = 1.0;
  ier = modelDriverCreate->ConvertUnit(
      fromLength, fromEnergy, fromCharge, fromTemperature, fromTime,
      requestedLengthUnit, requestedEnergyUnit, requestedChargeUnit,
      requestedTemperatureUnit, requestedTimeUnit,
      1.0, 0.0, 0.0, 0.0, 0.0,
      &convertLength);
  if (ier)
  {
    LOG_ERROR("Unable to convert length unit");
    return ier;
  }

  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    particleMass_[i] *= convertMass;       // convert to active units
    latticeConstant_[i] *= convertLength;  // convert to active units
  }

  // changing units of embedding function values
  // don't convert the density units (argument of embedding function)
  double convertEnergy = 1.0;
  ier = modelDriverCreate->ConvertUnit(
      fromLength, fromEnergy, fromCharge, fromTemperature, fromTime,
      requestedLengthUnit, requestedEnergyUnit, requestedChargeUnit,
      requestedTemperatureUnit, requestedTimeUnit,
      0.0, 1.0, 0.0, 0.0, 0.0,
      &convertEnergy);
  if (ier)
  {
    LOG_ERROR("Unable to convert energy unit");
    return ier;
  }

  if (convertEnergy != ONE)
  {
    for (int i = 0; i < numberModelSpecies_; ++i)
    {
      for (int j = 0; j < numberRhoPoints_; ++j)
      {
        embeddingData_[i][j] *= convertEnergy;
      }
    }
  }

  //
  // don't convert units of density (rho) --- they are ambiguous
  //

  // changing units of r*phi (stored in rPhiData) function values
  if (convertLength != ONE && convertEnergy != ONE)
  {
    for (int i = 0; i < numberModelSpecies_; ++i)
    {
      for (int j = 0; j < numberModelSpecies_; ++j)
      {
        for (int k = 0; k < numberRPoints_; ++k)
        {
          rPhiData_[i][j][k] *= convertLength*convertEnergy;
        }
      }
    }
  }

  // changing units of cutoff radius and deltaR
  if (convertLength != ONE)
  {
    cutoffParameter_ *= convertLength;
    deltaR_ *= convertLength;

    //
    // don't convert units of deltaRho --- they are ambiguous
    //
  }

  // register units
  ier = modelDriverCreate->SetUnits(
      requestedLengthUnit,
      requestedEnergyUnit,
      KIM::CHARGE_UNIT::unused,
      KIM::TEMPERATURE_UNIT::unused,
      requestedTimeUnit);
  if (ier)
  {
    LOG_ERROR("Unable to set units to requested values");
    return ier;
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int EAM_Implementation::RegisterKIMModelSettings(
    KIM::ModelDriverCreate * const modelDriverCreate) const
{
  // register numbering
  int error = modelDriverCreate->SetModelNumbering(
      KIM::NUMBERING::zeroBased);

  return error;
}

//******************************************************************************
#include "KIM_ModelComputeArgumentsCreateLogMacros.hpp"
int EAM_Implementation::RegisterKIMComputeArgumentsSettings(
    KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate) const
{
  // register arguments
  LOG_INFORMATION("Register argument supportStatus");
  int error =
      modelComputeArgumentsCreate->SetArgumentSupportStatus(
          KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,
          KIM::SUPPORT_STATUS::optional)
      || modelComputeArgumentsCreate->SetArgumentSupportStatus(
          KIM::COMPUTE_ARGUMENT_NAME::partialForces,
          KIM::SUPPORT_STATUS::optional)
      || modelComputeArgumentsCreate->SetArgumentSupportStatus(
          KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
          KIM::SUPPORT_STATUS::optional)
      || modelComputeArgumentsCreate->SetArgumentSupportStatus(
          KIM::COMPUTE_ARGUMENT_NAME::partialVirial,
          KIM::SUPPORT_STATUS::optional)
      || modelComputeArgumentsCreate->SetArgumentSupportStatus(
          KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial,
          KIM::SUPPORT_STATUS::optional);


  // register callbacks
  LOG_INFORMATION("Register callback supportStatus");
  error = error
      || modelComputeArgumentsCreate->SetCallbackSupportStatus(
          KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm,
          KIM::SUPPORT_STATUS::optional)
      || modelComputeArgumentsCreate->SetCallbackSupportStatus(
          KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term,
          KIM::SUPPORT_STATUS::optional);

  return error;
}

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int EAM_Implementation::RegisterKIMParameters(
    KIM::ModelDriverCreate * const modelDriverCreate,
    EAMFileType const eamFileType)
{
  int ier = false;

  ier = modelDriverCreate->SetParameterPointer(1, &cutoffParameter_, "cutoff");
  if (ier)
  {
    LOG_ERROR("Could not set register parameter 'cutoff'");
    return ier;
  }
  ier = modelDriverCreate->SetParameterPointer(1, &deltaRho_, "deltaRho");
  if (ier)
  {
    LOG_ERROR("Could not set register parameter 'deltaRho'");
    return ier;
  }
  ier = modelDriverCreate->SetParameterPointer(1, &deltaR_, "deltaR");
  if (ier)
  {
    LOG_ERROR("Could not set register parameter 'deltaR'");
    return ier;
  }

  for (int i = 0; i < numberModelSpecies_ ; i++)
  {
    for(int j = i; j < numberModelSpecies_; j++)
    {
      int indxPhi = i*numberModelSpecies_ + j - (i*i + i)/2;
      for(int k = 0; k < numberRPoints_; k++)
      {
        publish_rPhiData_[indxPhi][k] = rPhiData_[i][j][k];
      }
    }
  }

  for (int i = 0; i < numberModelSpecies_ ; i++)
  {
    for(int j = 0; j < numberModelSpecies_; j++)
    {
      int indxDensity = (eamFileType_ == FinnisSinclair) ?
          i*numberModelSpecies_ + j : i;
      for(int k = 0; k < numberRPoints_; k++)
      {
        publishDensityData_[indxDensity][k] = densityData_[i][j][k];
      }
      if (eamFileType_ != FinnisSinclair) break;
    }
  }

  ier = modelDriverCreate->SetParameterPointer(
      numberModelSpecies_ * numberRhoPoints_,
      embeddingData_[0], "embeddingData");
  if (ier)
  {
    LOG_ERROR("Could not set register parameter 'embeddingData'");
    return ier;
  }

  ier = modelDriverCreate->SetParameterPointer(
      numberUniqueSpeciesPairs_ * numberRPoints_,
      publish_rPhiData_[0], "rPhiData");
  if (ier)
  {
    LOG_ERROR("Could not set register parameter 'rPhiData'");
    return ier;
  }

  int const publishDensityLength = (eamFileType_ == FinnisSinclair) ?
      numberModelSpecies_*numberModelSpecies_ : numberModelSpecies_;
  ier = modelDriverCreate->SetParameterPointer(
      publishDensityLength * numberRPoints_,
      publishDensityData_[0], "densityData");
  if (ier)
  {
    LOG_ERROR("Could not set register parameter 'densityData'");
    return ier;
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int EAM_Implementation::RegisterKIMFunctions(
    KIM::ModelDriverCreate * const modelDriverCreate)
    const
{
  int error;

  // register the destroy() and reinit() functions
  error = modelDriverCreate->SetDestroyPointer(
      KIM::LANGUAGE_NAME::cpp, (KIM::func*) &(EAM::Destroy))
      || modelDriverCreate->SetRefreshPointer(
          KIM::LANGUAGE_NAME::cpp, (KIM::func*) &(EAM::Refresh))
      || modelDriverCreate->SetComputePointer(
          KIM::LANGUAGE_NAME::cpp, (KIM::func*) &(EAM::Compute))
      || modelDriverCreate->SetComputeArgumentsCreatePointer(
          KIM::LANGUAGE_NAME::cpp,
          (KIM::func*) &(EAM::ComputeArgumentsCreate))
      || modelDriverCreate->SetComputeArgumentsDestroyPointer(
          KIM::LANGUAGE_NAME::cpp,
          (KIM::func*) &(EAM::ComputeArgumentsDestroy));

  return error;
}

//******************************************************************************
template<class ModelObj>
int EAM_Implementation::SetRefreshMutableValues(
    ModelObj * const modelObj)
{ // use (possibly) new values of parameters to compute other quantities
  // NOTE: This function is templated because it's called with both a
  //       modelDriverCreate object during initialization and with a
  //       modelRefresh object when the Model's parameters have been altered
  int ier;

  // Update
  influenceDistance_ = cutoffParameter_;
  modelObj->SetInfluenceDistancePointer(&influenceDistance_);
  modelObj->SetNeighborListCutoffsPointer(1, &influenceDistance_);

  // update EAM_Implementation values
  cutoffSq_ = cutoffParameter_ * cutoffParameter_;
  oneByDr_ = ONE / deltaR_;
  oneByDrho_ = ONE / deltaRho_;

  // calculate spline interpolation coefficients
  SplineInterpolateAllData();

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
void EAM_Implementation::SplineInterpolateAllData()
{
  for (int i = 0; i < numberModelSpecies_; ++i)
  {
    SplineInterpolate(embeddingData_[i], deltaRho_, numberRhoPoints_,
                      embeddingCoeff_[i]);
    for (int j = 0; j < numberModelSpecies_; ++j)
    {
      SplineInterpolate(densityData_[i][j], deltaR_, numberRPoints_,
                        densityCoeff_[i][j]);
      SplineInterpolate(rPhiData_[i][j], deltaR_, numberRPoints_,
                        rPhiCoeff_[i][j]);
    }
  }
}

//******************************************************************************
#include "KIM_ModelComputeArgumentsLogMacros.hpp"
int EAM_Implementation::SetComputeMutableValues(
    KIM::ModelComputeArguments const * const modelComputeArguments,
    bool& isComputeProcess_dEdr,
    bool& isComputeProcess_d2Edr2,
    bool& isComputeEnergy,
    bool& isComputeForces,
    bool& isComputeParticleEnergy,
    bool& isComputeVirial,
    bool& isComputeParticleVirial,
    int const*& particleSpeciesCodes,
    int const*& particleContributing,
    VectorOfSizeDIM const*& coordinates,
    double*& energy,
    double*& particleEnergy,
    VectorOfSizeDIM*& forces,
    VectorOfSizeSix*& virial,
    VectorOfSizeSix*& particleVirial)
{
  int ier = true;

  // get compute flags
  int compProcess_dEdr;
  int compProcess_d2Edr2;

  modelComputeArguments->IsCallbackPresent(
      KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm,
      &compProcess_dEdr);
  modelComputeArguments->IsCallbackPresent(
      KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term,
      &compProcess_d2Edr2);

  isComputeProcess_dEdr = compProcess_dEdr;
  isComputeProcess_d2Edr2 = compProcess_d2Edr2;

  int const* numberOfParticles;
  ier =
      modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::numberOfParticles,
          &numberOfParticles)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes,
          &particleSpeciesCodes)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::particleContributing,
          &particleContributing)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::coordinates,
          (double const ** const) &coordinates)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,
          &energy)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
          &particleEnergy)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialForces,
          (double const ** const) &forces)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialVirial,
          (double const ** const) &virial)
      || modelComputeArguments->GetArgumentPointer(
          KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial,
          (double const ** const) &particleVirial);
  if (ier)
  {
    LOG_ERROR("GetArgumentPointer");
    return ier;
  }

  isComputeEnergy = (energy != 0);
  isComputeParticleEnergy = (particleEnergy != 0);
  isComputeForces = (forces != 0);
  isComputeVirial = (virial != 0);
  isComputeParticleVirial = (particleVirial != 0);

  // allocate memory if needed
  if (cachedNumberOfParticles_ < *numberOfParticles)
  {
    delete [] densityValue_;  // ok to delete null pointer
    densityValue_ = new double[*numberOfParticles];

    delete [] embeddingDerivativeValue_;  // ok to delete null pointer
    embeddingDerivativeValue_ = new double[*numberOfParticles];

    delete [] embeddingSecondDerivativeValue_;  // ok to delete null pointer
    embeddingSecondDerivativeValue_ = new double[*numberOfParticles];
  }

  // update values
  cachedNumberOfParticles_ = *numberOfParticles;

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
#include "KIM_ModelComputeLogMacros.hpp"
int EAM_Implementation::CheckParticleSpeciesCodes(
    KIM::ModelCompute const * const modelCompute,
    int const* const particleSpeciesCodes)
    const
{
  int ier;
  for (int i = 0; i < cachedNumberOfParticles_; ++i)
  {
    if ((particleSpeciesCodes[i] < 0) ||
        (particleSpeciesCodes[i] >= numberModelSpecies_))
    {
      ier = true;
      LOG_ERROR("unsupported particle species codes detected");
      return ier;
    }
  }

  // everything is good
  ier = false;
  return ier;
}

//******************************************************************************
int EAM_Implementation::GetComputeIndex(
    const bool& isComputeProcess_dEdr,
    const bool& isComputeProcess_d2Edr2,
    const bool& isComputeEnergy,
    const bool& isComputeForces,
    const bool& isComputeParticleEnergy,
    const bool& isComputeVirial,
    const bool& isComputeParticleVirial) const
{
  //const int processdE = 2;
  const int processd2E = 2;
  const int energy = 2;
  const int force = 2;
  const int particleEnergy = 2;
  const int virial = 2;
  const int particleVirial = 2;


  int index = 0;

  // processdE
  index += (int(isComputeProcess_dEdr))
      * processd2E * energy * force * particleEnergy * virial
      * particleVirial;

  // processd2E
  index += (int(isComputeProcess_d2Edr2))
      * energy * force * particleEnergy * virial * particleVirial;

  // energy
  index += (int(isComputeEnergy))
      * force * particleEnergy * virial * particleVirial;

  // force
  index += (int(isComputeForces))
      * particleEnergy * virial * particleVirial;

  // particleEnergy
  index += (int(isComputeParticleEnergy))
      * virial * particleVirial;

  // virial
  index += (int(isComputeVirial))
      * particleVirial;

  // particleVirial
  index += (int(isComputeParticleVirial));

  return index;
}

//******************************************************************************
void EAM_Implementation::ProcessVirialTerm(
    const double& dEidr,
    const double& rij,
    const double* const r_ij,
    const int& i,
    const int& j,
    VectorOfSizeSix virial) const
{
  double const v = dEidr/rij;

  virial[0] += v * r_ij[0] * r_ij[0];
  virial[1] += v * r_ij[1] * r_ij[1];
  virial[2] += v * r_ij[2] * r_ij[2];
  virial[3] += v * r_ij[1] * r_ij[2];
  virial[4] += v * r_ij[0] * r_ij[2];
  virial[5] += v * r_ij[0] * r_ij[1];
}

//******************************************************************************
void EAM_Implementation::ProcessParticleVirialTerm(
    const double& dEidr,
    const double& rij,
    const double* const r_ij,
    const int& i,
    const int& j,
    VectorOfSizeSix* const particleVirial) const
{
  double const v = dEidr/rij;
  VectorOfSizeSix vir;

  vir[0] = 0.5 * v * r_ij[0] * r_ij[0];
  vir[1] = 0.5 * v * r_ij[1] * r_ij[1];
  vir[2] = 0.5 * v * r_ij[2] * r_ij[2];
  vir[3] = 0.5 * v * r_ij[1] * r_ij[2];
  vir[4] = 0.5 * v * r_ij[0] * r_ij[2];
  vir[5] = 0.5 * v * r_ij[0] * r_ij[1];

  for (int k = 0; k < 6; ++k)
  {
    particleVirial[i][k] += vir[k];
    particleVirial[j][k] += vir[k];
  }
}

//==============================================================================
//
// Implementation of helper functions
//
//==============================================================================

//******************************************************************************
void AllocateAndInitialize2DArray(double**& arrayPtr, int const extentZero,
                                  int const extentOne)
{ // allocate memory and set pointers
  arrayPtr = new double*[extentZero];
  arrayPtr[0] = new double[extentZero * extentOne];
  for (int i = 1; i < extentZero; ++i)
  {
    arrayPtr[i] = arrayPtr[i-1] + extentOne;
  }

  // initialize
  for (int i = 0; i < extentZero; ++i)
  {
    for (int j = 0; j < extentOne; ++j)
    {
      arrayPtr[i][j] = 0.0;
    }
  }
}

//******************************************************************************
void Deallocate2DArray(double**& arrayPtr)
{ // deallocate memory
  if (arrayPtr != 0) delete [] arrayPtr[0];
  delete [] arrayPtr;

  // nullify pointer
  arrayPtr = 0;
}

//******************************************************************************
void AllocateAndInitialize3DArray(double***& arrayPtr, int const extentZero,
                                  int const extentOne, int const extentTwo)
{ // allocate memory and set pointers
  arrayPtr = new double**[extentZero];
  arrayPtr[0] = new double*[extentZero * extentOne];
  arrayPtr[0][0] = new double[extentZero * extentOne * extentTwo];
  for (int j = 1; j < extentZero; ++j)
  {
    arrayPtr[j] = arrayPtr[j-1] + extentOne;
    arrayPtr[0][j] = arrayPtr[0][j-1] + extentTwo;
  }
  for (int i = 1; i < extentZero; ++i)
  {
    arrayPtr[i][0] = arrayPtr[i-1][extentOne-1] + extentTwo;
    for (int j = 1; j < extentOne; ++j)
    {
      arrayPtr[i][j] = arrayPtr[i][j-1] + extentTwo;
    }
  }

  // initialize
  for (int i = 0; i < extentZero; ++i)
  {
    for (int j = 0; j < extentOne; ++j)
    {
      for (int k = 0; k < extentTwo; ++k)
      {
        arrayPtr[i][j][k] = 0.0;
      }
    }
  }
}

//******************************************************************************
void Deallocate3DArray(double***& arrayPtr)
{ // deallocate memory
  if (arrayPtr != 0)
  {
    delete [] arrayPtr[0][0];
    delete [] arrayPtr[0];
  }
  delete [] arrayPtr;

  // nullify pointer
  arrayPtr = 0;
}
