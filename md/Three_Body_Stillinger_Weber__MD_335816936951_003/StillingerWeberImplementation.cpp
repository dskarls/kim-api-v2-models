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
// Copyright (c) 2018, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Mingjian Wen
//


#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>

#include "StillingerWeberImplementation.hpp"
#include "KIM_ModelDriverHeaders.hpp"

#define MAXLINE 1024


//==============================================================================
//
// Implementation of StillingerWeberImplementation public member functions
//
//==============================================================================

//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
StillingerWeberImplementation::StillingerWeberImplementation(
  KIM::ModelDriverCreate* const modelDriverCreate,
  KIM::LengthUnit const requestedLengthUnit,
  KIM::EnergyUnit const requestedEnergyUnit,
  KIM::ChargeUnit const requestedChargeUnit,
  KIM::TemperatureUnit const requestedTemperatureUnit,
  KIM::TimeUnit const requestedTimeUnit,
  int* const ier)
  : numberModelSpecies_(0),
  numberUniqueSpeciesPairs_(0),
  cutoff_(NULL),
  A_(NULL),
  B_(NULL),
  p_(NULL),
  q_(NULL),
  sigma_(NULL),
  lambda_(NULL),
  gamma_(NULL),
  costheta0_(NULL),
  influenceDistance_(0.0),
  paddingNeighborHints_(1),
  halfListHints_(0),
  cutoffSq_2D_(NULL),
  A_2D_(NULL),
  B_2D_(NULL),
  p_2D_(NULL),
  q_2D_(NULL),
  sigma_2D_(NULL),
  lambda_2D_(NULL),
  gamma_2D_(NULL),
  costheta0_2D_(NULL),
  cachedNumberOfParticles_(0)
{
  FILE* parameterFilePointers[MAX_PARAMETER_FILES];
  int numberParameterFiles;

  modelDriverCreate->GetNumberOfParameterFiles(&numberParameterFiles);
  *ier = OpenParameterFiles(modelDriverCreate, numberParameterFiles,
      parameterFilePointers);
  if (*ier) {
    return;
  }

  *ier = ProcessParameterFiles(modelDriverCreate, numberParameterFiles,
      parameterFilePointers);
  CloseParameterFiles(numberParameterFiles, parameterFilePointers);
  if (*ier) {
    return;
  }

  *ier = ConvertUnits(modelDriverCreate,
      requestedLengthUnit,
      requestedEnergyUnit,
      requestedChargeUnit,
      requestedTemperatureUnit,
      requestedTimeUnit);
  if (*ier) {
    return;
  }

  *ier = SetRefreshMutableValues(modelDriverCreate);
  if (*ier) {
    return;
  }

  *ier = RegisterKIMModelSettings(modelDriverCreate);
  if (*ier) {
    return;
  }

  *ier = RegisterKIMParameters(modelDriverCreate);
  if (*ier) {
    return;
  }

  *ier = RegisterKIMFunctions(modelDriverCreate);
  if (*ier) {
    return;
  }

  // everything is good
  *ier = false;
  return;
}


//******************************************************************************
StillingerWeberImplementation::~StillingerWeberImplementation()
{ // note: it is ok to delete a null pointer and we have ensured that
  // everything is initialized to null

  Deallocate1DArray<double>(A_);
  Deallocate1DArray<double>(B_);
  Deallocate1DArray<double>(p_);
  Deallocate1DArray<double>(q_);
  Deallocate1DArray<double>(sigma_);
  Deallocate1DArray<double>(lambda_);
  Deallocate1DArray<double>(gamma_);
  Deallocate1DArray<double>(costheta0_);
  Deallocate1DArray<double>(cutoff_);

  Deallocate2DArray<double>(A_2D_);
  Deallocate2DArray<double>(B_2D_);
  Deallocate2DArray<double>(p_2D_);
  Deallocate2DArray<double>(q_2D_);
  Deallocate2DArray<double>(sigma_2D_);
  Deallocate2DArray<double>(lambda_2D_);
  Deallocate2DArray<double>(gamma_2D_);
  Deallocate2DArray<double>(costheta0_2D_);
  Deallocate2DArray<double>(cutoffSq_2D_);
}


//******************************************************************************
#include "KIM_ModelRefreshLogMacros.hpp"
int StillingerWeberImplementation::Refresh(KIM::ModelRefresh* const modelRefresh)
{
  int ier;

  ier = SetRefreshMutableValues(modelRefresh);
  if (ier) {
    return ier;
  }

  // nothing else to do for this case

  // everything is good
  ier = false;
  return ier;
}


//******************************************************************************
int StillingerWeberImplementation::Compute(
    KIM::ModelCompute const* const modelCompute,
    KIM::ModelComputeArguments const* const modelComputeArguments)
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
  int const* particleSpeciesCodes = NULL;
  int const* particleContributing = NULL;
  VectorOfSizeDIM const* coordinates = NULL;
  //
  // KIM API Model Output
  double* energy = NULL;
  double* particleEnergy = NULL;
  VectorOfSizeDIM* forces = NULL;
  VectorOfSizeSix* virial = NULL;
  VectorOfSizeSix* particleVirial = NULL;

  ier = SetComputeMutableValues(modelComputeArguments,
      isComputeProcess_dEdr, isComputeProcess_d2Edr2,
      isComputeEnergy, isComputeForces, isComputeParticleEnergy,
      isComputeVirial, isComputeParticleVirial,
      particleSpeciesCodes, particleContributing, coordinates,
      energy, forces, particleEnergy, virial, particleVirial);
  if (ier) {
    return ier;
  }

  // Skip this check for efficiency
  //
  //ier = CheckParticleSpecies(modelComputeArguments, particleSpeciesCodes);
  // if (ier) return ier;


#include "StillingerWeberImplementationComputeDispatch.cpp"
  return ier;
}


//******************************************************************************
int StillingerWeberImplementation::ComputeArgumentsCreate(
    KIM::ModelComputeArgumentsCreate* const modelComputeArgumentsCreate) const
{
  int ier;

  ier = RegisterKIMComputeArgumentsSettings(modelComputeArgumentsCreate);
  if (ier) {
    return ier;
  }

  // nothing else to do for this case

  // everything is good
  ier = false;
  return ier;
}


//******************************************************************************
int StillingerWeberImplementation::ComputeArgumentsDestroy(
    KIM::ModelComputeArgumentsDestroy* const modelComputeArgumentsDestroy)
const
{
  int ier;

  // nothing else to do for this case

  // everything is good
  ier = false;
  return ier;
}


//==============================================================================
//
// Implementation of StillingerWeberImplementation private member functions
//
//==============================================================================

//******************************************************************************
void StillingerWeberImplementation::AllocatePrivateParameterMemory()
{
  // nothing to do for this case
}


//******************************************************************************
void StillingerWeberImplementation::AllocateParameterMemory()
{ // allocate memory for data
  AllocateAndInitialize1DArray<double> (cutoff_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (A_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (B_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (p_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (q_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (sigma_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (lambda_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (gamma_, numberUniqueSpeciesPairs_);
  AllocateAndInitialize1DArray<double> (costheta0_, numberUniqueSpeciesPairs_);

  AllocateAndInitialize2DArray<double> (cutoffSq_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (A_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (B_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (p_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (q_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (sigma_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (lambda_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (gamma_2D_, numberModelSpecies_, numberModelSpecies_);
  AllocateAndInitialize2DArray<double> (costheta0_2D_, numberModelSpecies_, numberModelSpecies_);
}


//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int StillingerWeberImplementation::OpenParameterFiles(
    KIM::ModelDriverCreate* const modelDriverCreate,
    int const numberParameterFiles,
    FILE* parameterFilePointers[MAX_PARAMETER_FILES])
{
  int ier;

  if (numberParameterFiles > MAX_PARAMETER_FILES) {
    ier = true;
    LOG_ERROR("StillingerWeber given too many parameter files");
    return ier;
  }

  for (int i = 0; i < numberParameterFiles; ++i) {
    std::string const* paramFileName;
    ier = modelDriverCreate->GetParameterFileName(i, &paramFileName);
    if (ier) {
      LOG_ERROR("Unable to get parameter file name");
      return ier;
    }

    parameterFilePointers[i] = fopen(paramFileName->c_str(), "r");
    if (parameterFilePointers[i] == 0) {
      char message[MAXLINE];
      sprintf(message,
          "StillingerWeber parameter file number %d cannot be opened",
          i);
      ier = true;
      LOG_ERROR(message);
      for (int j = i - 1; i <= 0; --i) {
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
int StillingerWeberImplementation::ProcessParameterFiles(
    KIM::ModelDriverCreate* const modelDriverCreate,
    int const numberParameterFiles,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES])
{
  int N, ier;
  int endOfFileFlag = 0;
  char spec1[MAXLINE], spec2[MAXLINE], nextLine[MAXLINE];
  int iIndex, jIndex, indx;
  double next_A, next_B, next_p, next_q, next_sigma, next_lambda, next_gamma;
  double next_costheta0, next_cutoff;

  getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  ier = sscanf(nextLine, "%d", &N);
  if (ier != 1) {
    sprintf(nextLine, "unable to read first line of the parameter file");
    ier = true;
    LOG_ERROR(nextLine);
    fclose(parameterFilePointers[0]);
    return ier;
  }
  numberModelSpecies_ = N;
  numberUniqueSpeciesPairs_ = ((numberModelSpecies_ + 1) * numberModelSpecies_) / 2;
  AllocateParameterMemory();

  // set all values of p_ to -1.1e10 for later check that we have read all params
  for (int i = 0; i < ((N + 1) * N / 2); i++) {
    p_[i] = -1.1e10;
  }

  // keep track of known species
  std::map<KIM::SpeciesName const, int, KIM::SPECIES_NAME::Comparator> modelSpeciesMap;
  int index = 0;   // species code integer code starting from 0

  // Read and process data lines
  getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  while (endOfFileFlag == 0)
  {
    ier = sscanf(nextLine, "%s %s %lg %lg %lg %lg %lg %lg %lg %lg %lg",
        spec1, spec2, &next_A, &next_B, &next_p, &next_q, &next_sigma,
        &next_lambda, &next_gamma, &next_costheta0, &next_cutoff);
    if (ier != 11) {
      sprintf(nextLine, "error reading lines of the parameter file");
      LOG_ERROR(nextLine);
      return true;
    }

    // convert species strings to proper type instances
    KIM::SpeciesName const specName1(spec1);
    KIM::SpeciesName const specName2(spec2);
     if ((specName1.String() == "unknown") ||
         (specName2.String() == "unknown") ) {
      sprintf(nextLine, "error parameter file: get unknown species");
      LOG_ERROR(nextLine);
      return true;
    }


    // check for new species
    std::map<KIM::SpeciesName const, int, KIM::SPECIES_NAME::Comparator>::
    const_iterator iIter = modelSpeciesMap.find(specName1);
    if (iIter == modelSpeciesMap.end()) {
      modelSpeciesMap[specName1] = index;
      modelSpeciesCodeList_.push_back(index);

      ier = modelDriverCreate->SetSpeciesCode(specName1, index);
      if (ier) {
        return ier;
      }
      iIndex = index;
      index++;
    }
    else {
      iIndex = modelSpeciesMap[specName1];
    }

    std::map<KIM::SpeciesName const, int, KIM::SPECIES_NAME::Comparator>::
    const_iterator jIter = modelSpeciesMap.find(specName2);
    if (jIter == modelSpeciesMap.end()) {
      modelSpeciesMap[specName2] = index;
      modelSpeciesCodeList_.push_back(index);

      ier = modelDriverCreate->SetSpeciesCode(specName2, index);
      if (ier) {
        return ier;
      }
      jIndex = index;
      index++;
    }
    else {
      jIndex = modelSpeciesMap[specName2];
    }

    if (iIndex >= jIndex) {
      indx = jIndex * N + iIndex - (jIndex * jIndex + jIndex) / 2;
    }
    else {
      indx = iIndex * N + jIndex - (iIndex * iIndex + iIndex) / 2;
    }
    A_[indx] = next_A;
    B_[indx] = next_B;
    p_[indx] = next_p;
    q_[indx] = next_q;
    sigma_[indx] = next_sigma;
    lambda_[indx] = next_lambda;
    gamma_[indx] = next_gamma;
    costheta0_[indx] = next_costheta0;
    cutoff_[indx] = next_cutoff;

    getNextDataLine(parameterFilePointers[0], nextLine, MAXLINE, &endOfFileFlag);
  }

  // check we have read all parameters
  for (int i = 0; i < ((N + 1) * N / 2); i++) {
    if (p_[i] < -1e10) {
      sprintf(nextLine, "error: not enough parameter data.\n");
      sprintf(nextLine, "%d species requires %d data lines.", N, (N + 1) * N / 2);
      LOG_ERROR(nextLine);
      return true;
    }
  }

  // everything is good
  ier = false;
  return ier;
}


//******************************************************************************
void StillingerWeberImplementation::getNextDataLine(
    FILE* const filePtr, char* nextLinePtr, int const maxSize,
    int* endOfFileFlag)
{
  do
  {
    if (fgets(nextLinePtr, maxSize, filePtr) == NULL) {
      *endOfFileFlag = 1;
      break;
    }

    while ((nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t') ||
           (nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r'))
    {
      nextLinePtr = (nextLinePtr + 1);
    }
  } while ((strncmp("#", nextLinePtr, 1) == 0) || (strlen(nextLinePtr) == 0));
}


//******************************************************************************
void StillingerWeberImplementation::CloseParameterFiles(
    int const numberParameterFiles,
    FILE* const parameterFilePointers[MAX_PARAMETER_FILES])
{
  for (int i = 0; i < numberParameterFiles; ++i) {
    fclose(parameterFilePointers[i]);
  }
}


//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int StillingerWeberImplementation::ConvertUnits(
    KIM::ModelDriverCreate* const modelDriverCreate,
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

  // changing units of sigma, gamma, and cutoff
  double convertLength = 1.0;

  ier = modelDriverCreate->ConvertUnit(
      fromLength, fromEnergy, fromCharge, fromTemperature, fromTime,
      requestedLengthUnit, requestedEnergyUnit, requestedChargeUnit,
      requestedTemperatureUnit, requestedTimeUnit,
      1.0, 0.0, 0.0, 0.0, 0.0,
      &convertLength);
  if (ier) {
    LOG_ERROR("Unable to convert length unit");
    return ier;
  }
  // convert to active units
  if (convertLength != ONE) {
    for (int i = 0; i < numberUniqueSpeciesPairs_; ++i) {
      sigma_[i] *= convertLength;
      gamma_[i] *= convertLength;
      cutoff_[i] *= convertLength;
    }
  }

  // changing units of A and lambda
  double convertEnergy = 1.0;
  ier = modelDriverCreate->ConvertUnit(
      fromLength, fromEnergy, fromCharge, fromTemperature, fromTime,
      requestedLengthUnit, requestedEnergyUnit, requestedChargeUnit,
      requestedTemperatureUnit, requestedTimeUnit,
      0.0, 1.0, 0.0, 0.0, 0.0,
      &convertEnergy);
  if (ier) {
    LOG_ERROR("Unable to convert energy unit");
    return ier;
  }
  // convert to active units
  if (convertLength != ONE) {
    for (int i = 0; i < numberUniqueSpeciesPairs_; ++i) {
      A_[i] *= convertEnergy;
      lambda_[i] *= convertEnergy;
    }
  }

  // register units
  ier = modelDriverCreate->SetUnits(
      requestedLengthUnit,
      requestedEnergyUnit,
      KIM::CHARGE_UNIT::unused,
      KIM::TEMPERATURE_UNIT::unused,
      KIM::TIME_UNIT::unused);
  if (ier) {
    LOG_ERROR("Unable to set units to requested values");
    return ier;
  }

  // everything is good
  ier = false;
  return ier;
}


//******************************************************************************
int StillingerWeberImplementation::RegisterKIMModelSettings(
    KIM::ModelDriverCreate* const modelDriverCreate) const
{
  // register numbering
  int error = modelDriverCreate->SetModelNumbering(KIM::NUMBERING::zeroBased);

  return error;
}


//******************************************************************************
#include "KIM_ModelComputeArgumentsCreateLogMacros.hpp"
int StillingerWeberImplementation::RegisterKIMComputeArgumentsSettings(
    KIM::ModelComputeArgumentsCreate* const modelComputeArgumentsCreate) const
{
  // register arguments
  LOG_INFORMATION("Register argument supportStatus");

  int error =
    modelComputeArgumentsCreate->SetArgumentSupportStatus(
        KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,
        KIM::SUPPORT_STATUS::optional) ||
    modelComputeArgumentsCreate->SetArgumentSupportStatus(
        KIM::COMPUTE_ARGUMENT_NAME::partialForces,
        KIM::SUPPORT_STATUS::optional) ||
    modelComputeArgumentsCreate->SetArgumentSupportStatus(
        KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
        KIM::SUPPORT_STATUS::optional) ||
    modelComputeArgumentsCreate->SetArgumentSupportStatus(
        KIM::COMPUTE_ARGUMENT_NAME::partialVirial,
        KIM::SUPPORT_STATUS::optional) ||
    modelComputeArgumentsCreate->SetArgumentSupportStatus(
        KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial,
        KIM::SUPPORT_STATUS::optional);

  // register callbacks
  LOG_INFORMATION("Register callback supportStatus");
  error =
    error ||
    modelComputeArgumentsCreate->SetCallbackSupportStatus(
        KIM::COMPUTE_CALLBACK_NAME::ProcessDEDrTerm,
        KIM::SUPPORT_STATUS::optional) ||
    modelComputeArgumentsCreate->SetCallbackSupportStatus(
        KIM::COMPUTE_CALLBACK_NAME::ProcessD2EDr2Term,
        KIM::SUPPORT_STATUS::optional);

  return error;
}


//******************************************************************************
#include "KIM_ModelDriverCreateLogMacros.hpp"
int StillingerWeberImplementation::RegisterKIMParameters(
    KIM::ModelDriverCreate* const modelDriverCreate)
{
  int ier = false;

  // publish parameters (order is important)
  ier =
    modelDriverCreate->SetParameterPointer(
        numberUniqueSpeciesPairs_, A_, "A") ||
    modelDriverCreate->SetParameterPointer(
        numberUniqueSpeciesPairs_, B_, "B") ||
    modelDriverCreate->SetParameterPointer(
        numberUniqueSpeciesPairs_, p_, "p") ||
    modelDriverCreate->SetParameterPointer(
        numberUniqueSpeciesPairs_, q_, "q") ||
    modelDriverCreate->SetParameterPointer(
        numberUniqueSpeciesPairs_, sigma_, "sigma") ||
    modelDriverCreate->SetParameterPointer(
        numberUniqueSpeciesPairs_, lambda_, "lambda") ||
    modelDriverCreate->SetParameterPointer(
        numberUniqueSpeciesPairs_, gamma_, "gamma") ||
    modelDriverCreate->SetParameterPointer(
        numberUniqueSpeciesPairs_, costheta0_, "costheta0") ||
    modelDriverCreate->SetParameterPointer(
        numberUniqueSpeciesPairs_, cutoff_, "cutoff");
  if (ier) {
    LOG_ERROR("set_parameters");
    return ier;
  }

  // everything is good
  ier = false;
  return ier;
}


//******************************************************************************
int StillingerWeberImplementation::RegisterKIMFunctions(
    KIM::ModelDriverCreate* const modelDriverCreate)
const
{
  int error;

  // register the Destroy(), Refresh(), and Compute() functions
  error =
    modelDriverCreate->SetDestroyPointer(
        KIM::LANGUAGE_NAME::cpp,
        (KIM::func*)&(StillingerWeber::Destroy)) ||
    modelDriverCreate->SetRefreshPointer(
        KIM::LANGUAGE_NAME::cpp,
        (KIM::func*)&(StillingerWeber::Refresh)) ||
    modelDriverCreate->SetComputePointer(
        KIM::LANGUAGE_NAME::cpp,
        (KIM::func*)&(StillingerWeber::Compute)) ||
    modelDriverCreate->SetComputeArgumentsCreatePointer(
        KIM::LANGUAGE_NAME::cpp,
        (KIM::func*)&(StillingerWeber::ComputeArgumentsCreate)) ||
    modelDriverCreate->SetComputeArgumentsDestroyPointer(
        KIM::LANGUAGE_NAME::cpp,
        (KIM::func*)&(StillingerWeber::ComputeArgumentsDestroy));

  return error;
}


//******************************************************************************
template<class ModelObj>
int StillingerWeberImplementation::SetRefreshMutableValues(
    ModelObj* const modelObj)
{ // use (possibly) new values of parameters to compute other quantities
  // NOTE: This function is templated because it's called with both a
  //       modelDriverCreate object during initialization and with a
  //       modelRefresh object when the Model's parameters have been altered
  int ier;

  // update parameters
  for (int i = 0; i < numberModelSpecies_; ++i) {
    for (int j = 0; j <= i; ++j) {
      int const index = j * numberModelSpecies_ + i - (j * j + j) / 2;
      A_2D_[i][j] = A_2D_[j][i] = A_[index];
      B_2D_[i][j] = B_2D_[j][i] = B_[index];
      p_2D_[i][j] = p_2D_[j][i] = p_[index];
      q_2D_[i][j] = q_2D_[j][i] = q_[index];
      sigma_2D_[i][j] = sigma_2D_[j][i] = sigma_[index];
      lambda_2D_[i][j] = lambda_2D_[j][i] = lambda_[index];
      gamma_2D_[i][j] = gamma_2D_[j][i] = gamma_[index];
      costheta0_2D_[i][j] = costheta0_2D_[j][i] = costheta0_[index];
      cutoffSq_2D_[i][j] = cutoffSq_2D_[j][i] = cutoff_[index] * cutoff_[index];
    }
  }

  // update cutoff value in KIM API object
  influenceDistance_ = 0.0;

  for (int i = 0; i < numberModelSpecies_; i++) {
    int indexI = modelSpeciesCodeList_[i];

    for (int j = 0; j < numberModelSpecies_; j++) {
      int indexJ = modelSpeciesCodeList_[j];

      if (influenceDistance_ < cutoffSq_2D_[indexI][indexJ]) {
        influenceDistance_ = cutoffSq_2D_[indexI][indexJ];
      }
    }
  }

  influenceDistance_ = sqrt(influenceDistance_);
  modelObj->SetInfluenceDistancePointer(&influenceDistance_);
  modelObj->SetNeighborListPointers(1,
      &influenceDistance_, &paddingNeighborHints_, &halfListHints_);

  // everything is good
  ier = false;
  return ier;
}


//******************************************************************************
#include "KIM_ModelComputeArgumentsLogMacros.hpp"
int StillingerWeberImplementation::SetComputeMutableValues(
    KIM::ModelComputeArguments const* const modelComputeArguments,
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
    VectorOfSizeDIM*& forces,
    double*& particleEnergy,
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
        &numberOfParticles) ||
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::particleSpeciesCodes,
        &particleSpeciesCodes) ||
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::particleContributing,
        &particleContributing) ||
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::coordinates,
        (double const** const)&coordinates) ||
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::partialEnergy,
        &energy) ||
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::partialForces,
        (double const** const)&forces) ||
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::partialParticleEnergy,
        &particleEnergy) ||
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::partialVirial,
        (double const** const)&virial) ||
    modelComputeArguments->GetArgumentPointer(
        KIM::COMPUTE_ARGUMENT_NAME::partialParticleVirial,
        (double const** const)&particleVirial);
  if (ier) {
    LOG_ERROR("GetArgumentPointer");
    return ier;
  }

  isComputeEnergy = (energy != NULL);
  isComputeForces = (forces != NULL);
  isComputeParticleEnergy = (particleEnergy != NULL);
  isComputeVirial = (virial != NULL);
  isComputeParticleVirial = (particleVirial != NULL);

  // update values
  cachedNumberOfParticles_ = *numberOfParticles;

  // everything is good
  ier = false;
  return ier;
}


//******************************************************************************
// Assume that the particle species interge code starts from 0
#include "KIM_ModelComputeLogMacros.hpp"
int StillingerWeberImplementation::CheckParticleSpeciesCodes(
    KIM::ModelCompute const* const modelCompute,
    int const* const particleSpeciesCodes) const
{
  int ier;

  for (int i = 0; i < cachedNumberOfParticles_; ++i) {
    if ((particleSpeciesCodes[i] < 0) || (particleSpeciesCodes[i] >= numberModelSpecies_)) {
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
int StillingerWeberImplementation::GetComputeIndex(
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
           * processd2E * energy * force * particleEnergy * virial * particleVirial;

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


//==============================================================================
//
// Stillinger-Weber functions
//
//==============================================================================
void StillingerWeberImplementation::CalcPhiTwo(int const ispec, int const jspec,
    double const r, double& phi) const
{
  // get parameters
  double const A = A_2D_[ispec][jspec];
  double const B = B_2D_[ispec][jspec];
  double const p = p_2D_[ispec][jspec];
  double const q = q_2D_[ispec][jspec];
  double const sigma = sigma_2D_[ispec][jspec];
  double const cutoff = sqrt(cutoffSq_2D_[ispec][jspec]);

  double r_cap = r / sigma;

  if (r >= cutoff) {
    phi = 0.0;
  }
  else {
    phi = A * (B * pow(r_cap, -p) - pow(r_cap, -q)) * exp(sigma / (r - cutoff));
  }
}


void StillingerWeberImplementation::CalcPhiDphiTwo(int const ispec, int const jspec,
    double const r, double& phi, double& dphi) const
{
  // get parameters
  double const A = A_2D_[ispec][jspec];
  double const B = B_2D_[ispec][jspec];
  double const p = p_2D_[ispec][jspec];
  double const q = q_2D_[ispec][jspec];
  double const sigma = sigma_2D_[ispec][jspec];
  double const cutoff = sqrt(cutoffSq_2D_[ispec][jspec]);

  double r_cap = r / sigma;

  if (r >= cutoff) {
    phi = 0.0;
    dphi = 0.0;
  }
  else {
    phi = A * (B * pow(r_cap, -p) - pow(r_cap, -q)) * exp(sigma / (r - cutoff));

    dphi = (q * pow(r_cap, -(q + 1)) - p * B * pow(r_cap, -(p + 1)))
           - (B * pow(r_cap, -p) - pow(r_cap, -q)) * pow((r - cutoff) / sigma, -2);
    dphi *= (1 / sigma) * A * exp(sigma / (r - cutoff));
  }
}


void StillingerWeberImplementation::CalcPhiD2phiTwo(int const ispec, int const jspec,
    double const r, double& phi, double& dphi, double& d2phi) const
{
  // get parameters
  double const A = A_2D_[ispec][jspec];
  double const B = B_2D_[ispec][jspec];
  double const p = p_2D_[ispec][jspec];
  double const q = q_2D_[ispec][jspec];
  double const sigma = sigma_2D_[ispec][jspec];
  double const cutoff = sqrt(cutoffSq_2D_[ispec][jspec]);

  double r_cap = r / sigma;

  if (r >= cutoff) {
    phi = 0.0;
    dphi = 0.0;
    d2phi = 0.0;
  }
  else {
    phi = A * (B * pow(r_cap, -p) - pow(r_cap, -q)) * exp(sigma / (r - cutoff));

    dphi = (q * pow(r_cap, -(q + 1)) - p * B * pow(r_cap, -(p + 1)))
           - (B * pow(r_cap, -p) - pow(r_cap, -q)) * pow((r - cutoff) / sigma, -2);
    dphi *= (1 / sigma) * A * exp(sigma / (r - cutoff));

    d2phi = (B * pow(r_cap, -p) - pow(r_cap, -q))
            * (pow((r - cutoff) / sigma, -4) + 2 * pow((r - cutoff) / sigma, -3))
            + 2 * (p * B * pow(r_cap, -(p + 1)) - q * pow(r_cap, -(q + 1)))
            * pow((r - cutoff) / sigma, -2)
            + (p * (p + 1) * B * pow(r_cap, -(p + 2))
               - q * (q + 1) * pow(r_cap, -(q + 2)));
    d2phi *= (1 / (sigma * sigma)) * A * exp(sigma / (r - cutoff));
  }
}


void StillingerWeberImplementation::CalcPhiThree(int const ispec, int const jspec,
    int const kspec, double const rij, double const rik, double const rjk,
    double& phi) const
{
  // get parameters
  double const lambda_ij = lambda_2D_[ispec][jspec];
  double const lambda_ik = lambda_2D_[ispec][kspec];
  double const gamma_ij = gamma_2D_[ispec][jspec];
  double const gamma_ik = gamma_2D_[ispec][kspec];
  double const costheta0_ij = costheta0_2D_[ispec][jspec];
  double const cutoff_ij = sqrt(cutoffSq_2D_[ispec][jspec]);
  double const cutoff_ik = sqrt(cutoffSq_2D_[ispec][kspec]);
  // mix parameters
  double const lambda = sqrt(fabs(lambda_ij) * fabs(lambda_ik));
  double const costheta0 = costheta0_ij;  // do not mix

  if (rij < cutoff_ij && rik < cutoff_ik) {
    double costhetajik = (pow(rij, 2) + pow(rik, 2) - pow(rjk, 2)) / (2 * rij * rik);
    double diff_costhetajik = costhetajik - costheta0;
    double exp_ij_ik = exp(gamma_ij / (rij - cutoff_ij) + gamma_ik / (rik - cutoff_ik));
    phi = lambda * exp_ij_ik * diff_costhetajik * diff_costhetajik;
  }
  else {
    phi = 0.0;
  }
}


void StillingerWeberImplementation::CalcPhiDphiThree(int const ispec, int const jspec,
    int const kspec, double const rij, double const rik, double const rjk,
    double& phi, double* const dphi) const
{
  // get parameters
  double const lambda_ij = lambda_2D_[ispec][jspec];
  double const lambda_ik = lambda_2D_[ispec][kspec];
  double const gamma_ij = gamma_2D_[ispec][jspec];
  double const gamma_ik = gamma_2D_[ispec][kspec];
  double const costheta0_ij = costheta0_2D_[ispec][jspec];
  double const cutoff_ij = sqrt(cutoffSq_2D_[ispec][jspec]);
  double const cutoff_ik = sqrt(cutoffSq_2D_[ispec][kspec]);
  // mix parameters
  double const lambda = sqrt(fabs(lambda_ij) * fabs(lambda_ik));
  double const costheta0 = costheta0_ij;  // do not mix


  if (rij < cutoff_ij && rik < cutoff_ik) {
    double costhetajik = (pow(rij, 2) + pow(rik, 2) - pow(rjk, 2)) / (2 * rij * rik);
    double diff_costhetajik = costhetajik - costheta0;

    /* Derivatives of cosines w.r.t rij, rik, rjk */
    double costhetajik_ij = (pow(rij, 2) - pow(rik, 2) + pow(rjk, 2))
                            / (2 * rij * rij * rik);
    double costhetajik_ik = (pow(rik, 2) - pow(rij, 2) + pow(rjk, 2))
                            / (2 * rij * rik * rik);
    double costhetajik_jk = -rjk / (rij * rik);

    /* Variables for simplifying terms */
    double exp_ij_ik = exp(gamma_ij / (rij - cutoff_ij) + gamma_ik / (rik - cutoff_ik));
    double d_ij = -gamma_ij* pow(rij - cutoff_ij, -2);
    double d_ik = -gamma_ik* pow(rik - cutoff_ik, -2);

    phi = lambda * exp_ij_ik * diff_costhetajik * diff_costhetajik;

    dphi[0] = lambda * diff_costhetajik * exp_ij_ik
              * (d_ij * diff_costhetajik + 2 * costhetajik_ij);
    dphi[1] = lambda * diff_costhetajik * exp_ij_ik
              * (d_ik * diff_costhetajik + 2 * costhetajik_ik);
    dphi[2] = lambda * diff_costhetajik * exp_ij_ik * 2 * costhetajik_jk;
  }
  else {
    phi = 0.0;
    dphi[0] = 0.0;
    dphi[1] = 0.0;
    dphi[2] = 0.0;
  }
}


// Calculate phi_three(rij, rik, rjk) and its 1st & 2nd derivatives
// dphi_three(rij, rik, rjk), d2phi_three(rij, rik, rjk)
//
// dphi has three components as derivatives of phi w.r.t. rij, rik, rjk
//
// d2phi as symmetric Hessian matrix of phi has six components:
//    [0]=(ij,ij), [3]=(ij,ik), [4]=(ij,jk)
//                 [1]=(ik,ik), [5]=(ik,jk)
//                              [2]=(jk,jk)

void StillingerWeberImplementation::CalcPhiD2phiThree(int const ispec, int const jspec,
    int const kspec, double const rij, double const rik, double const rjk,
    double& phi, double* const dphi, double* const d2phi) const
{
  // get parameters
  double const lambda_ij = lambda_2D_[ispec][jspec];
  double const lambda_ik = lambda_2D_[ispec][kspec];
  double const gamma_ij = gamma_2D_[ispec][jspec];
  double const gamma_ik = gamma_2D_[ispec][kspec];
  double const costheta0_ij = costheta0_2D_[ispec][jspec];
  double const cutoff_ij = sqrt(cutoffSq_2D_[ispec][jspec]);
  double const cutoff_ik = sqrt(cutoffSq_2D_[ispec][kspec]);
  // mix parameters
  double const lambda = sqrt(fabs(lambda_ij) * fabs(lambda_ik));
  double const costheta0 = costheta0_ij;  // do not mix


  if (rij < cutoff_ij && rik < cutoff_ik) {
    double costhetajik = (pow(rij, 2) + pow(rik, 2) - pow(rjk, 2)) / (2 * rij * rik);
    double diff_costhetajik = costhetajik - costheta0;
    double diff_costhetajik_2 = diff_costhetajik * diff_costhetajik;

    /* Derivatives of cosines w.r.t. r_ij, r_ik, r_jk */
    double costhetajik_ij = (pow(rij, 2) - pow(rik, 2) + pow(rjk, 2))
                            / (2 * rij * rij * rik);
    double costhetajik_ik = (pow(rik, 2) - pow(rij, 2) + pow(rjk, 2))
                            / (2 * rij * rik * rik);
    double costhetajik_jk = -rjk / (rij * rik);

    /* Hessian matrix of cosine */
    double costhetajik_ij_ij = (pow(rik, 2) - pow(rjk, 2)) / (rij * rij * rij * rik);
    double costhetajik_ik_ik = (pow(rij, 2) - pow(rjk, 2)) / (rij * rik * rik * rik);
    double costhetajik_jk_jk = -1 / (rij * rik);
    double costhetajik_ij_ik = -(pow(rij, 2) + pow(rik, 2) + pow(rjk, 2))
                               / (2 * rij * rij * rik * rik);
    double costhetajik_ij_jk = rjk / (rij * rij * rik);
    double costhetajik_ik_jk = rjk / (rik * rik * rij);

    /* Variables for simplifying terms */
    double exp_ij_ik = exp(gamma_ij / (rij - cutoff_ij) + gamma_ik / (rik - cutoff_ik));
    double d_ij = -gamma_ij* pow(rij - cutoff_ij, -2);
    double d_ik = -gamma_ik* pow(rik - cutoff_ik, -2);
    double d_ij_2 = d_ij * d_ij;
    double d_ik_2 = d_ik * d_ik;
    double dd_ij = 2* gamma_ij* pow(rij - cutoff_ij, -3);
    double dd_ik = 2* gamma_ik* pow(rik - cutoff_ik, -3);

    phi = lambda * exp_ij_ik * diff_costhetajik * diff_costhetajik;

    dphi[0] = lambda * diff_costhetajik * exp_ij_ik
              * (d_ij * diff_costhetajik + 2 * costhetajik_ij);
    dphi[1] = lambda * diff_costhetajik * exp_ij_ik
              * (d_ik * diff_costhetajik + 2 * costhetajik_ik);
    dphi[2] = lambda * diff_costhetajik * exp_ij_ik * 2 * costhetajik_jk;

    d2phi[0] = lambda * exp_ij_ik *
               ((d_ij_2 + dd_ij) * diff_costhetajik_2
                + (4 * d_ij * costhetajik_ij + 2 * costhetajik_ij_ij) * diff_costhetajik
                + 2 * costhetajik_ij * costhetajik_ij);
    d2phi[1] = lambda * exp_ij_ik *
               ((d_ik_2 + dd_ik) * diff_costhetajik_2
                + (4 * d_ik * costhetajik_ik + 2 * costhetajik_ik_ik) * diff_costhetajik
                + 2 * costhetajik_ik * costhetajik_ik);
    d2phi[2] = lambda * 2 * exp_ij_ik *
               (costhetajik_jk_jk * diff_costhetajik
                + costhetajik_jk * costhetajik_jk);
    d2phi[3] = lambda * exp_ij_ik *
               (d_ij * d_ik * diff_costhetajik_2
                + (d_ij * costhetajik_ik + d_ik * costhetajik_ij + costhetajik_ij_ik)
                * 2 * diff_costhetajik + 2 * costhetajik_ij * costhetajik_ik);
    d2phi[4] = lambda * exp_ij_ik *
               ((d_ij * costhetajik_jk + costhetajik_ij_jk)
                * 2 * diff_costhetajik + 2 * costhetajik_ij * costhetajik_jk);
    d2phi[5] = lambda * exp_ij_ik *
               ((d_ik * costhetajik_jk + costhetajik_ik_jk)
                * 2 * diff_costhetajik + 2 * costhetajik_ik * costhetajik_jk);
  }
  else {
    phi = 0.0;
    dphi[0] = dphi[1] = dphi[2] = 0.0;
    d2phi[0] = d2phi[1] = d2phi[2] = d2phi[3] = d2phi[4] = d2phi[5] = 0.0;
  }
}
