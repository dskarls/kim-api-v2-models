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


#ifndef EAM_IMPLEMENTATION_HPP_
#define EAM_IMPLEMENTATION_HPP_

#include <map>
#include <cmath>

#include "KIM_LogVerbosity.hpp"
#include "EAM.hpp"

#define MAXLINE 1024
#define DIMENSION 3
#define ONE 1.0
#define TWO 2.0
#define HALF 0.5

#define MAX_PARAMETER_FILES 20
#define NUMBER_SETFL_COMMENT_LINES 3

#include "EAM_Spline.hpp"

//==============================================================================
//
// Type definitions, enumerations, and helper function prototypes
//
//==============================================================================

// type declaration for get neighbor functions
typedef int (GetNeighborFunction)(void const * const, int const,
                                  int * const, int const ** const);                                  
// type declaration for vector of constant dimension
typedef double VectorOfSizeDIM[DIMENSION];
// type declaration for funcfl data
struct SetOfFuncflData
{
  int numberRhoPoints[MAX_PARAMETER_FILES];
  double deltaRho[MAX_PARAMETER_FILES];
  int numberRPoints[MAX_PARAMETER_FILES];
  double deltaR[MAX_PARAMETER_FILES];
  double cutoff[MAX_PARAMETER_FILES];
  double* embeddingData[MAX_PARAMETER_FILES];
  double* densityData[MAX_PARAMETER_FILES];
  double* ZData[MAX_PARAMETER_FILES];
};

// enumeration for EAMFileType
enum EAMFileType {Setfl, Funcfl, FinnisSinclair, Error};

// helper routine declarations
void AllocateAndInitialize2DArray(double**& arrayPtr, int const extentZero,
                                  int const extentOne);
void Deallocate2DArray(double**& arrayPtr);
void AllocateAndInitialize3DArray(double***& arrayPtr, int const extentZero,
                                  int const extentOne, int const extentTwo);
void Deallocate3DArray(double***& arrayPtr);

//==============================================================================
//
// Declaration of EAM_Implementation class
//
//==============================================================================

//******************************************************************************
class EAM_Implementation
{
 public:
  EAM_Implementation(
      KIM::ModelDriverCreate * const modelDriverCreate,
      KIM::LengthUnit const requestedLengthUnit,
      KIM::EnergyUnit const requestedEnergyUnit,
      KIM::ChargeUnit const requestedChargeUnit,
      KIM::TemperatureUnit const requestedTemperatureUnit,
      KIM::TimeUnit const requestedTimeUnit,
      int * const ier);
  ~EAM_Implementation();  // no explicit Destroy() needed here

  int Refresh(KIM::ModelRefresh * const modelRefresh);
  int Compute(
      KIM::ModelCompute const * const modelCompute,
      KIM::ModelComputeArguments const * const modelComputeArguments);
  int ComputeArgumentsCreate(
      KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate)
      const;
  int ComputeArgumentsDestroy(
      KIM::ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy)
      const;

 private:
  // Constant values that never change
  //   Set in constructor (via SetConstantValues)
  //
  //
  // EAM_Implementation: constants
  int numberModelSpecies_;
  int numberUniqueSpeciesPairs_;
  EAMFileType eamFileType_;


  // Constant values that are read from the input files and never change
  //   Set in constructor (via functions listed below)
  //
  //
  // KIM API: Model Fixed Parameters
  //   Memory allocated in   AllocateFixedParameterMemory()
  //   Memory deallocated in ~EAM_Implementation()
  //   Data set in ReadParameterFile routines
  char* comments_ptr_[MAX_PARAMETER_FILES];
  char comments_[MAX_PARAMETER_FILES][MAXLINE];
  char particleNames_[MAXLINE];
  int* particleNumber_;
  double* particleMass_;
  double* latticeConstant_;
  char** latticeType_;
  int numberRhoPoints_;
  int numberRPoints_;
  //
  // KIM API: Model Free Parameters whose (pointer) values never change
  //   Memory allocated in   AllocateFreeParameterMemory() (from constructor)
  //   Memory deallocated in ~EAM_Implementation()
  //   Data set in ReadParameterFile routines OR by KIM Simulator
  double** embeddingData_;
  double*** densityData_;
  double*** rPhiData_;

  // Free Parameter pointers to be published without repeat data
  double** publishDensityData_;
  double** publish_rPhiData_;

  // Mutable values that only change when reinit() executes
  //   Set in Reinit (via SetReinitMutableValues)
  //
  //
  // KIM API: Model Free Parameters (can be changed directly by KIM Simulator)
  double influenceDistance_;
  double cutoffParameter_;
  double deltaR_;
  double deltaRho_;
  //
  // EAM_Implementation: values (changed only by Refresh())
  double cutoffSq_;
  double oneByDr_;
  double oneByDrho_;
  //   Memory allocated once by AllocateFreeParameterMemory()
  double** embeddingCoeff_;
  double*** densityCoeff_;
  double*** rPhiCoeff_;

  // Mutable values that can change with each call to Refresh() and Compute()
  //   Memory may be reallocated on each call
  //
  //
  // EAM_Implementation: values that change
  int cachedNumberOfParticles_;
  double* densityValue_;
  double* embeddingDerivativeValue_;
  double* embeddingSecondDerivativeValue_;


  // Helper methods
  //
  //
  // Related to constructor
  void AllocateFixedParameterMemory();
  static int OpenParameterFiles(
      KIM::ModelDriverCreate * const modelDriverCreate,
      int const numberParameterFiles,
      FILE* parameterFilePointers[MAX_PARAMETER_FILES]);
  static EAMFileType DetermineParameterFileTypes(
      KIM::ModelDriverCreate * const modelDriverCreate,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
      int const numberParameterFiles);
  static EAMFileType IsFuncflOrSetfl(FILE* const fptr);
  static EAMFileType IsSetflOrFinnisSinclair(
      KIM::ModelDriverCreate * const modelDriverCreate,
      FILE* const fptr);
  int ProcessParameterFileHeaders(
      KIM::ModelDriverCreate * const modelDriverCreate,
      EAMFileType const eamFileType,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
      int const numberParameterFiles, SetOfFuncflData& funcflData);
  int ReadSetflHeader(
      KIM::ModelDriverCreate * const modelDriverCreate,
      FILE* const fptr);
  int ReadFuncflHeader(
      KIM::ModelDriverCreate * const modelDriverCreate,
      FILE* const fptr,
      int const fileIndex, int& numberRhoPoints,
      double& deltaRho, int& numberRPoints,
      double& deltaR, double& cutoffParameter);
  int SetParticleNamesForFuncflModels(
      KIM::ModelDriverCreate * const modelDriverCreate);
  void AllocateFreeParameterMemory();
  int ProcessParameterFileData(
      KIM::ModelDriverCreate * const modelDriverCreate,
      EAMFileType const eamFileType,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
      int const numberParameterFiles, SetOfFuncflData& funcflData);
  int ReadSetflData(
      KIM::ModelDriverCreate * const modelDriverCreate,
      FILE* const fptr);
  int ReadFinnisSinclairData(
      KIM::ModelDriverCreate * const modelDriverCreate,
      FILE* const fptr);
  static int ReadFuncflData(
      KIM::ModelDriverCreate * const modelDriverCreate,
      FILE* const fptr,
      int const fileIndex,
      SetOfFuncflData& funcflData);
  static int GrabData(
      KIM::ModelDriverCreate * const modelDriverCreate,
      FILE* const fptr, int const n, double* const list);
  void ReinterpolateAndMix(SetOfFuncflData const& funcflData);
  static void CloseParameterFiles(
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES],
      int const numberParameterFiles);
  int ConvertUnits(
      KIM::ModelDriverCreate * const modelDriverCreate,
      KIM::LengthUnit const requestedLengthUnit,
      KIM::EnergyUnit const requestedEnergyUnit,
      KIM::ChargeUnit const requestedChargeUnit,
      KIM::TemperatureUnit const requestedTemperatureUnit,
      KIM::TimeUnit const requestedTimeUnit);
  int RegisterKIMModelSettings(
      KIM::ModelDriverCreate * const modelDriverCreate) const;
  int RegisterKIMComputeArgumentsSettings(
      KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate) const;
  int RegisterKIMParameters(
      KIM::ModelDriverCreate * const modelDriverCreate,
      EAMFileType const eamFileType);
  int RegisterKIMFunctions(
      KIM::ModelDriverCreate * const modelDriverCreate) const;
  //
  // Related to Refresh()
  template<class ModelObj>
  int SetReinitMutableValues(ModelObj * const modelObj);
  void SplineInterpolateAllData();
  static void SplineInterpolate(double const* const dat,
                                double const delta, int const n,
                                double* const coe);

  //
  // Related to Compute()
  int SetComputeMutableValues(
      KIM::ModelComputeArguments const * const modelComputeArguments,
      bool& isComputeProcess_dEdr,
      bool& isComputeProcess_d2Edr2,
      bool& isComputeEnergy,
      bool& isComputeForces,
      bool& isComputeParticleEnergy,
      int const*& particleSpeciesCodes,
      int const*& particleContributing,
      VectorOfSizeDIM const*& coordinates,
      double*& energy,
      double*& particleEnergy,
      VectorOfSizeDIM*& forces);
  int CheckParticleSpeciesCodes(KIM::ModelCompute const * const modelCompute,
                                int const* const particleSpeciesCodes) const;
  int GetComputeIndex(const bool& isComputeProcess_dEdr,
                      const bool& isComputeProcess_d2Edr2,
                      const bool& isComputeEnergy,
                      const bool& isComputeForces,
                      const bool& isComputeParticleEnergy) const;

  // compute functions
  template< bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
            bool isComputeEnergy, bool isComputeForces,
            bool isComputeParticleEnergy >
  int Compute(KIM::ModelCompute const * const modelCompute,
              KIM::ModelComputeArguments const * const modelComputeArguments,
              const int* const particleSpeciesCodes,
              const int* const particleContributing,
              const VectorOfSizeDIM* const coordinates,
              double* const energy,
              VectorOfSizeDIM* const forces,
              double* const particleEnergy) const;
};

//==============================================================================
//
// Definition of EAM_Implementation::Compute functions
//
// NOTE: Here we rely on the compiler optimizations to prune dead code
//       after the template expansions.  This provides high efficiency
//       and easy maintenance.
//
//==============================================================================

//******************************************************************************
#include "KIM_ModelComputeLogMacros.hpp"
template< bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
          bool isComputeEnergy, bool isComputeForces,
          bool isComputeParticleEnergy >
int EAM_Implementation::Compute(
    KIM::ModelCompute const * const modelCompute,
    KIM::ModelComputeArguments const * const modelComputeArguments,
    const int* const particleSpeciesCodes,
    const int* const particleContributing,
    const VectorOfSizeDIM* const coordinates,
    double* const energy,
    VectorOfSizeDIM* const forces,
    double* const particleEnergy) const
{
  int ier = false;

  if ((isComputeEnergy == false) &&
      (isComputeParticleEnergy == false) &&
      (isComputeForces == false) &&
      (isComputeProcess_dEdr == false) &&
      (isComputeProcess_d2Edr2 == false))
    return ier;

  // initialize electron density for each contributing particle
  for (int i = 0; i < cachedNumberOfParticles_; ++i)
  {
    if (particleContributing[i])
    {
      densityValue_[i] = 0.0;
      // no need to initialize embeddingDerivativeValue_
    }
  }
  if (isComputeEnergy == true)
  {
    *energy = 0.0;
  }
  if (isComputeForces == true)
  {
    for (int i = 0; i < cachedNumberOfParticles_; ++i)
    {
      for (int j = 0; j < DIMENSION; ++j)
        forces[i][j] = 0.0;
    }
  }

  // compute electron density
  // Setup loop over contributing particles
  int i = 0;
  int numnei = 0;
  int const * n1atom = 0;
  for (i = 0; i < cachedNumberOfParticles_; ++i)
  {
    if (particleContributing[i])
    {
      modelComputeArguments->GetNeighborList(0, i, &numnei, &n1atom);

      // Setup loop over neighbors of current particle
      for (int jj = 0; jj < numnei; ++jj)
      {
        int const j = n1atom[jj];
        double* r_ij;
        double r_ijValue[DIMENSION];
        // Compute r_ij
        r_ij = r_ijValue;
        for (int k = 0; k < DIMENSION; ++k)
          r_ij[k] = coordinates[j][k] - coordinates[i][k];

        // compute distance squared
        double rij2 = 0.0;
        for (int k = 0; k < DIMENSION; ++k)
          rij2 += r_ij[k] * r_ij[k];

        if (rij2 <= cutoffSq_)
        { // compute contribution to electron density
          double rijOffset;
          int rijIndex;
          double const rij = sqrt(rij2);

          // compute rijOffset and rijIndex
          GET_DELTAX_AND_INDEX(rij, oneByDr_, numberRPoints_, rijOffset,
                               rijIndex);

          // interpolate value of rho_beta(r_ij)
          double densityBetaValue;
          double const* const densityBetaCoeff
              = densityCoeff_[particleSpeciesCodes[j]][particleSpeciesCodes[i]];
          INTERPOLATE_F(densityBetaCoeff, rijOffset, rijIndex, densityBetaValue);
          densityValue_[i] += densityBetaValue;
        }
      }  // end of loop over neighbors

      densityValue_[i] = std::max(densityValue_[i], 0.0);  // ensure non-negative
      // Check for density too large
      double const rhoDomainLimit = (numberRhoPoints_ - 1.0)*deltaRho_;
      if (densityValue_[i] > rhoDomainLimit)
      {
        ier = true;
        LOG_ERROR("Particle has density value outside of embedding function"
            " interpolation domain");
        return ier;
      }
    } // end if statement as to whether particle is contributing
  }  // end of loop over contributing particles

  // calculate embedding function and its derivative
  for (i = 0; i < cachedNumberOfParticles_; ++i)
  {
    if (particleContributing[i])
    {
      double densityOffset;
      int densityIndex;
      // compute densityOffset and densityIndex
      GET_DELTAX_AND_INDEX(densityValue_[i], oneByDrho_, numberRhoPoints_,
                           densityOffset, densityIndex);
      double const* const embeddingAlphaCoeff
          = embeddingCoeff_[particleSpeciesCodes[i]];
      if (0 < numnei)
      {
        // interpolate F_i(rho_i)
        double embeddingValue;
        INTERPOLATE_F(embeddingAlphaCoeff, densityOffset, densityIndex,
                      embeddingValue);
        // Contribute embedding term to Energy
        if (isComputeEnergy == true)
        {
          *energy += embeddingValue;
        }
        // Contribute embedding term to ParticleEnergy
        if (isComputeParticleEnergy == true)
        {
          particleEnergy[i] = embeddingValue;
        }
      }
      // Compute embedding derivative
      if ((isComputeForces == true) || (isComputeProcess_dEdr == true) ||
          (isComputeProcess_d2Edr2 == true))
      {
        // interpolate dF_i(rho_i)/d(rho_i)
        INTERPOLATE_DF(embeddingAlphaCoeff, densityOffset, densityIndex,
                       embeddingDerivativeValue_[i]);
      }
      // Compute embedding second derivative
      if (isComputeProcess_d2Edr2 == true)
      {
        // interpolate d^2F_i(rho_i)/d(rho_i)^2
        INTERPOLATE_D2F(embeddingAlphaCoeff, densityOffset, densityIndex,
                        embeddingSecondDerivativeValue_[i]);
      }
    } // end if statement as to whether particle is contributing
  }

  // calculate contribution from electron density to the force part
  // and from pair function
  //
  // Setup loop over contributing particles
  for (i = 0; i < cachedNumberOfParticles_; ++i)
  {
    if (particleContributing[i])
    {
      modelComputeArguments->GetNeighborList(0, i, &numnei, &n1atom);

      // Setup loop over neighbors of current particle
      for (int jj = 0; jj < numnei; ++jj)
      {
        int const j = n1atom[jj];
        double* r_ij;
        double r_ijValue[DIMENSION];
        // Compute r_ij appropriately
        r_ij = r_ijValue;
        for (int k = 0; k < DIMENSION; ++k)
          r_ij[k] = coordinates[j][k] - coordinates[i][k];

        // compute distance squared
        double rij2 = 0.0;
        for (int k = 0; k < DIMENSION; ++k)
          rij2 += r_ij[k] * r_ij[k];

        if (rij2 <= cutoffSq_)
        { // compute contribution to energy and force
          double rijOffset;
          int rijIndex;
          double const rij = sqrt(rij2);

          // compute rijOffset and rijIndex
          GET_DELTAX_AND_INDEX(rij, oneByDr_, numberRPoints_, rijOffset,
                               rijIndex);

          // interpolate r_ij*phi(r_ij)
          double rijPhiValue;
          double const* const rijPhiAlphaBetaCoeff
              = rPhiCoeff_[particleSpeciesCodes[i]][particleSpeciesCodes[j]];
          INTERPOLATE_F(rijPhiAlphaBetaCoeff, rijOffset, rijIndex, rijPhiValue);

          // find phi(r_ij)
          double const oneByRij = ONE / rij;
          double const pairPotentialValue = rijPhiValue * oneByRij;

          // Contribute pair term to Energy as half or full
          if (isComputeEnergy == true)
          {
            *energy += HALF * pairPotentialValue;
          }

          // Contribute pair term to Particle Energy
          if (isComputeParticleEnergy == true)
          {
            particleEnergy[i] += HALF * pairPotentialValue;
          }

          // Compute dEdrByR terms as half or full
          double dEdrByRij = 0.0;
          if ((isComputeForces == true) || (isComputeProcess_dEdr == true))
          {
            // interpolate derivative of r_ij*phi(r_ij) function
            double rijPhiDerivativeValue;
            INTERPOLATE_DF(rijPhiAlphaBetaCoeff, rijOffset, rijIndex,
                           rijPhiDerivativeValue);

            // interpolate derivative of rho_beta(r_ij)
            double densityBetaDerivativeValue;
            double const* const densityBetaCoeff
                = densityCoeff_[particleSpeciesCodes[i]][particleSpeciesCodes[j]];
            INTERPOLATE_DF(densityBetaCoeff, rijOffset, rijIndex,
                           densityBetaDerivativeValue);

            // compute dEdr contribution
            // embedding contribution to dEdr
            double const embeddingContribution
                = (embeddingDerivativeValue_[i] * densityBetaDerivativeValue);

            // pair potential contribution
            double const pairPotentialContribution
                = HALF * (rijPhiDerivativeValue - pairPotentialValue)
                * oneByRij;

            // divide by r so we can multiply by r_ij below
            dEdrByRij = (embeddingContribution + pairPotentialContribution)
                * oneByRij;
          }

          // Contribute dEdrByR to forces
          if (isComputeForces == true)
          {
            for (int k = 0; k < DIMENSION; ++k)
            {
              forces[i][k] += dEdrByRij * r_ij[k];
              forces[j][k] -= dEdrByRij * r_ij[k];
            }
          }

          // Call process_dEdr
          if (isComputeProcess_dEdr == true)
          {
            double const rij = sqrt(rij2);
            double const dEidr = dEdrByRij*rij;
            ier = modelComputeArguments
                ->ProcessDEDrTerm(dEidr, rij, r_ij, i, j);
            if (ier)
            {
              LOG_ERROR("process_dEdr");
              return ier;
            }
          }
        }  // if particles i and j interact
      }  // end of first neighbor loop
    } // end of if statement for whether particle is contributing
  }  // end of loop over contributing particles

  // Separate loop nest for process_d2Edr2
  if (isComputeProcess_d2Edr2 == true)
  {
    // For this potential, the second derivative
    // d^2E/dr_{ij}dr_{kl} is nonzero only if i,j,k,l are
    // not distinct. As such, we need only address all
    // derivatives d^2E/dr_{ij}dr_{ik}.
    //
    // We simply iterate over particle i's neighbor list in a doubly-nested
    // triangular loop.

    // Setup loop over contributing particles
    for (i = 0; i < cachedNumberOfParticles_; ++i)
    {
      if (particleContributing[i])
      {
        modelComputeArguments->GetNeighborList(0, i, &numnei, &n1atom);

        // Setup loop over neighbors of current particle
        for (int jj = 0; jj < numnei; ++jj)
        {
          // adjust index of particle neighbor
          int const j = n1atom[jj];

          // Declare enough space to hold two displacements,
          // one for each neighbor. Ensuring that the two
          // displacements are contiguous in memory allows
          // the use of r_ijValue in the call to process_d2Edr2.
          double r_ijValue[2*DIMENSION];

          // Pointer to the first displacement.
          double* r_ij;

          // Compute r_ij appropriately
          r_ij = r_ijValue;
          for (int k = 0; k < DIMENSION; ++k)
            r_ij[k] = coordinates[j][k] - coordinates[i][k];

          // compute distance squared
          double rij2 = 0.0;
          for (int k = 0; k < DIMENSION; ++k)
            rij2 += r_ij[k] * r_ij[k];

          if (rij2 <= cutoffSq_)
          {
            double rijOffset;
            int rijIndex;
            double rij[2] = { sqrt(rij2), 0.0 };  // rij[1] used for ik pair

            // compute rijOffset and rijIndex
            GET_DELTAX_AND_INDEX(rij[0], oneByDr_, numberRPoints_, rijOffset,
                                 rijIndex);

            // interpolate r_ij*phi(r_ij)
            double rijPhiValue;
            double const* const rijPhiAlphaBetaCoeff
                = rPhiCoeff_[particleSpeciesCodes[i]][particleSpeciesCodes[j]];
            INTERPOLATE_F(rijPhiAlphaBetaCoeff, rijOffset, rijIndex, rijPhiValue);

            // find phi(r_ij)
            double const oneByRij = ONE / rij[0];
            double const pairPotentialValue = rijPhiValue * oneByRij;

            // interpolate derivative of r_ij*phi(r_ij) function
            double rijPhiDerivativeValue;
            INTERPOLATE_DF(rijPhiAlphaBetaCoeff, rijOffset, rijIndex,
                           rijPhiDerivativeValue);

            // find derivative of phi(r_ij)
            double const pairPotentialDerivativeValue
                = (rijPhiDerivativeValue - pairPotentialValue) * oneByRij;

            // interpolate derivative of rho_beta(r_ij)
            double densityBetaDerivativeValue;
            double const* const densityBetaCoeff
                = densityCoeff_[particleSpeciesCodes[i]][particleSpeciesCodes[j]];
            INTERPOLATE_DF(densityBetaCoeff, rijOffset, rijIndex,
                           densityBetaDerivativeValue);

            // Setup second loop over neighbors of current particle
            for (int kk = jj; kk < numnei; ++kk)
            {
              int const k = n1atom[kk];

              // adjust index of particle neighbor
              double* r_ik;

              // Pointer to the second displacement vector
              r_ik = &(r_ijValue[DIMENSION]);

              for (int d = 0; d < DIMENSION; ++d)
                r_ik[d] = coordinates[k][d] - coordinates[i][d];

              // compute distance squared
              double rik2 = 0.0;
              for (int d = 0; d < DIMENSION; ++d)
                rik2 += r_ik[d] * r_ik[d];

              if (rik2 <= cutoffSq_)
              {
                double rikOffset;
                int rikIndex;
                rij[1] = sqrt(rik2);

                // compute rikOffset and rikIndex
                GET_DELTAX_AND_INDEX(rij[1], oneByDr_, numberRPoints_, rikOffset,
                                     rikIndex);

                // interpolate derivative of rho_gamma(r_ik)
                double densityGammaDerivativeValue;
                double const* const densityGammaCoeff
                    = densityCoeff_[particleSpeciesCodes[i]][particleSpeciesCodes[k]];
                INTERPOLATE_DF(densityGammaCoeff, rikOffset, rikIndex,
                               densityGammaDerivativeValue);

                // mixed-index embedding contribution to d2Edr2
                double const mixedEmbeddingContribution
                    = embeddingSecondDerivativeValue_[i]
                    * densityBetaDerivativeValue
                    * densityGammaDerivativeValue;

                double d2Edr2 = mixedEmbeddingContribution;

                if (kk == jj)
                { // interpolate second derivative of r_ij*phi(r_ij) function
                  double rijPhiSecondDerivativeValue;
                  INTERPOLATE_D2F(rijPhiAlphaBetaCoeff, rijOffset, rijIndex,
                                  rijPhiSecondDerivativeValue);

                  // interpolate second derivative of rho_beta(r_ij)
                  double densityBetaSecondDerivativeValue;
                  double const* const densityBetaCoeff
                      = densityCoeff_[particleSpeciesCodes[i]][particleSpeciesCodes[j]];
                  INTERPOLATE_D2F(densityBetaCoeff, rijOffset, rijIndex,
                                  densityBetaSecondDerivativeValue);

                  // pair potential contribution
                  double const pairPotentialSecondDerivativeValue
                      = (rijPhiSecondDerivativeValue
                         - TWO * pairPotentialDerivativeValue) * oneByRij;
                  double const pairPotentialContribution
                      = HALF * pairPotentialSecondDerivativeValue;

                  // second derivative of embedding contribution
                  double const embeddingContribution
                      = embeddingDerivativeValue_[i]
                      * densityBetaSecondDerivativeValue;

                  d2Edr2 += pairPotentialContribution + embeddingContribution;
                }
                else
                { // Process transpose pair
                  const int iis[2] = { i, i };
                  const int jjs[2] = { k, j };
                  int const* const piis = &iis[0];
                  int const* const pjjs = &jjs[0];
                  const double rikrij[2] = { rij[1], rij[0] };
                  double const* const prikrij = &rikrij[0];
                  const double r_ikr_ij[6] = { r_ij[3], r_ij[4], r_ij[5],
                                               r_ij[0], r_ij[1], r_ij[2] };
                  double const* const pr_ikr_ij = &r_ikr_ij[0];

                  ier = modelComputeArguments
                        ->ProcessD2EDr2Term(d2Edr2, prikrij, pr_ikr_ij, piis, pjjs);
                  if (ier)
                  {
                    LOG_ERROR("process_d2Edr2");
                    return ier;
                  }
                }

                int const iis[2] = { i, i };
                int const jjs[2] = { j, k };
                int const* const piis = &iis[0];
                int const* const pjjs = &jjs[0];
                double const* const prij = &rij[0];
                ier = modelComputeArguments
                        ->ProcessD2EDr2Term(d2Edr2, prij, r_ij, piis, pjjs);
                if (ier)
                {
                  LOG_ERROR("process_d2Edr2");
                  return ier;
                }
              }  // if particles i and k interact
            }  // end of second neighbor loop
          }  // if particles i and j interact
        }  // end of first neighbor loop
      } // end of if statement for whether particle is contributing
    }  // end of loop over contributing particles
  }  // if (isComputeProcess_d2Edr2 == true)

  // everything is good
  ier = false;
  return ier;
}

#endif  // EAM_IMPLEMENTATION_HPP_
