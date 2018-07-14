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

#define MAX_NUMBER_OF_SPECIES 150
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
typedef double VectorOfSizeSix[6];
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
//type declaration for deferred neighbors
typedef struct
{
  int index;
} neighbor;
//type declaration for iterating over deferred neighbor lists
typedef std::multimap<int, neighbor>::iterator deferredNeighborIterator;

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
  // Private Model Parameters
  //   Data set in ReadParameterFile routines
  char* comments_ptr_[MAX_PARAMETER_FILES];
  char comments_[MAX_PARAMETER_FILES][MAXLINE];
  char particleNames_[MAXLINE];
  int particleNumber_[MAX_NUMBER_OF_SPECIES];
  double particleMass_[MAX_NUMBER_OF_SPECIES];
  double latticeConstant_[MAX_NUMBER_OF_SPECIES];
  char latticeType_[MAX_NUMBER_OF_SPECIES][MAXLINE];
  int numberRhoPoints_;
  int numberRPoints_;
  //
  // KIM API: Model Parameters whose (pointer) values never change
  //   Memory allocated in AllocateParameterMemory() (from constructor)
  //   Memory deallocated in destructor
  //   Data set in ReadParameterFile routines OR by KIM Simulator
  double** embeddingData_;
  double*** densityData_;
  double*** rPhiData_;

  // Parameter pointers to be published without repeated data
  double** publishDensityData_;
  double** publish_rPhiData_;

  // Mutable values that only change when Refresh() executes
  //   Set in Refresh (via SetRefreshMutableValues)
  //
  //
  // KIM API: Model Parameters (can be changed directly by KIM Simulator)
  double influenceDistance_;
  double cutoffParameter_;
  double deltaR_;
  double deltaRho_;

  // EAM_Implementation: values (changed only by Refresh())
  double cutoffSq_;
  double oneByDr_;
  double oneByDrho_;

  // Memory allocated once by AllocateParameterMemory()
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

  // Hints passed to the simulator by the Model
  int paddingNeighborHints_; // Whether to request that the simulator
                             // omit neighbors of padding atoms
  int halfListHints_; // Whether to request half-lists

  // Helper methods
  //
  //
  // Related to constructor
  void AllocateParameterMemory();
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
      int const numberParameterFiles,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES]);
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
      KIM::ModelComputeArgumentsCreate * const modelComputeArgumentsCreate)
      const;
  int RegisterKIMParameters(
      KIM::ModelDriverCreate * const modelDriverCreate,
      EAMFileType const eamFileType);
  int RegisterKIMFunctions(
      KIM::ModelDriverCreate * const modelDriverCreate) const;
  //
  // Related to Refresh()
  template<class ModelObj>
  int SetRefreshMutableValues(ModelObj * const modelObj);
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
      bool& isComputeVirial,
      bool& isComputeParticleVirial,
      int const*& particleSpeciesCodes,
      int const*& particleContributing,
      VectorOfSizeDIM const*& coordinates,
      double*& energy,
      double*& particleEnergy,
      VectorOfSizeDIM*& forces,
      VectorOfSizeSix*& virial,
      VectorOfSizeSix*& particleViral);
  int CheckParticleSpeciesCodes(KIM::ModelCompute const * const modelCompute,
                                int const* const particleSpeciesCodes) const;
  int GetComputeIndex(const bool& isComputeProcess_dEdr,
                      const bool& isComputeProcess_d2Edr2,
                      const bool& isComputeEnergy,
                      const bool& isComputeForces,
                      const bool& isComputeParticleEnergy,
                      const bool& isComputeVirial,
                      const bool& isComputeParticleVirial) const;
  void ProcessVirialTerm(const double& dEidr,
                         const double& rij,
                         const double* const r_ij,
                         const int& i,
                         const int& j,
                         VectorOfSizeSix virial) const;
  void ProcessParticleVirialTerm(const double& dEidr,
                                 const double& rij,
                                 const double* const r_ij,
                                 const int& i,
                                 const int& j,
                                 VectorOfSizeSix* const particleVirial) const;

  // compute functions
  template< bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
            bool isComputeEnergy, bool isComputeForces,
            bool isComputeParticleEnergy, bool isComputeVirial,
            bool isComputeParticleVirial >
  int Compute(KIM::ModelCompute const * const modelCompute,
              KIM::ModelComputeArguments const * const modelComputeArguments,
              const int* const particleSpeciesCodes,
              const int* const particleContributing,
              const VectorOfSizeDIM* const coordinates,
              double* const energy,
              VectorOfSizeDIM* const forces,
              double* const particleEnergy,
              VectorOfSizeSix virial,
              VectorOfSizeSix* const particleVirial) const;
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
          bool isComputeParticleEnergy, bool isComputeVirial,
          bool isComputeParticleVirial >
int EAM_Implementation::Compute(
    KIM::ModelCompute const * const modelCompute,
    KIM::ModelComputeArguments const * const modelComputeArguments,
    const int* const particleSpeciesCodes,
    const int* const particleContributing,
    const VectorOfSizeDIM* const coordinates,
    double* const energy,
    VectorOfSizeDIM* const forces,
    double* const particleEnergy,
    VectorOfSizeSix virial,
    VectorOfSizeSix* const particleVirial) const
{
  int ier = false;

  if ((isComputeEnergy == false) &&
      (isComputeParticleEnergy == false) &&
      (isComputeForces == false) &&
      (isComputeProcess_dEdr == false) &&
      (isComputeProcess_d2Edr2 == false) &&
      (isComputeVirial == false) &&
      (isComputeParticleVirial == false))
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
  if (isComputeVirial == true)
  {
    for (int i = 0; i < 6; ++i) virial[i] = 0.0;
  }
  if (isComputeForces == true)
  {
    for (int i = 0; i < cachedNumberOfParticles_; ++i)
    {
      for (int j = 0; j < DIMENSION; ++j)
        forces[i][j] = 0.0;
    }
  }
  if (isComputeParticleVirial == true)
  {
    int const cachedNumParticles = cachedNumberOfParticles_;
    for (int i = 0; i < cachedNumParticles; ++i)
    {
      for (int j = 0; j < 6; ++j)
        particleVirial[i][j] = 0.0;
    }
  }

  // compute electron density
  // Setup loop over contributing particles
  int i = 0;
  int numnei = 0;
  int const * n1atom = NULL;
  for (i = 0; i < cachedNumberOfParticles_; ++i)
  {
    if (particleContributing[i])
    {
      modelComputeArguments->GetNeighborList(0, i, &numnei, &n1atom);

      // Setup loop over neighbors of current particle
      for (int jj = 0; jj < numnei; ++jj)
      {
        int const j = n1atom[jj];

        if (i < j) // Effective half list
        {
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
            INTERPOLATE_F(densityBetaCoeff, rijOffset, rijIndex,
                          densityBetaValue);
            densityValue_[i] += densityBetaValue;
            if (particleContributing[j])
            {
              double densityAlphaValue;
              double const* const densityAlphaCoeff
                = densityCoeff_[particleSpeciesCodes[i]][particleSpeciesCodes[j]];
              INTERPOLATE_F(densityAlphaCoeff, rijOffset, rijIndex,
                            densityAlphaValue);
              densityValue_[j] += densityAlphaValue;
            }
          }
        } // end effective half-list check (i < j)
      }  // end of loop over neighbors

      // ensure non-negative
      densityValue_[i] = std::max(densityValue_[i], 0.0);
      // Check for density too large
      double const rhoDomainLimit = (numberRhoPoints_ - 1.0)*deltaRho_;
      if (densityValue_[i] > rhoDomainLimit)
      {
        ier = true;
        LOG_ERROR("Particle has density value outside of embedding function"
                  " interpolation domain");
        return ier;
      }
    }  // end if statement as to whether particle is contributing
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
    }  // end if statement as to whether particle is contributing
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
        if (i < j) // Effective half list
        {
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

            if (particleContributing[j])
            {
              // Contribute pair term to Energy
              if (isComputeEnergy == true)
              {
                *energy += pairPotentialValue;
              }

              // Contribute pair term to Particle Energy
              if (isComputeParticleEnergy == true)
              {
                particleEnergy[i] += HALF * pairPotentialValue;
                particleEnergy[j] += HALF * pairPotentialValue;
              }
            }
            else
            {
              // Contribute pair term to Energy
              if (isComputeEnergy == true)
              {
                *energy += HALF * pairPotentialValue;
              }

              // Contribute pair term to Particle Energy
              if (isComputeParticleEnergy == true)
              {
                particleEnergy[i] += HALF * pairPotentialValue;
              }
            }

            // Compute dEdrByR terms
            double dEdrByRij = 0.0;
            if ((isComputeForces == true) || (isComputeProcess_dEdr == true))
            {
              // interpolate derivative of r_ij*phi(r_ij) function
              double rijPhiDerivativeValue;
              INTERPOLATE_DF(rijPhiAlphaBetaCoeff, rijOffset, rijIndex,
                             rijPhiDerivativeValue);

              // interpolate derivative of rho_beta(r_ij)
              double densityBetaDerivativeValue;
              double const* const densityBetaCoeff =
                  densityCoeff_[particleSpeciesCodes[j]][particleSpeciesCodes[i]];
              INTERPOLATE_DF(densityBetaCoeff, rijOffset, rijIndex,
                             densityBetaDerivativeValue);
              if (particleContributing[j])
              {
                double densityAlphaDerivativeValue;
                double const* const densityAlphaCoeff =
                    densityCoeff_[particleSpeciesCodes[i]][particleSpeciesCodes[j]];
                INTERPOLATE_DF(densityAlphaCoeff, rijOffset, rijIndex,
                             densityAlphaDerivativeValue);

                // compute dEdr contribution
                // embedding contribution to dEdr
                double const embeddingContribution
                    = ((embeddingDerivativeValue_[i] * densityBetaDerivativeValue)
                        +
                       (embeddingDerivativeValue_[j] * densityAlphaDerivativeValue));

                // pair potential contribution
                double const pairPotentialContribution
                    = (rijPhiDerivativeValue - pairPotentialValue) * oneByRij;

                // divide by r so we can multiply by r_ij below
                dEdrByRij = (embeddingContribution + pairPotentialContribution)
                    * oneByRij;
              }
              else
              {
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
            if ((isComputeProcess_dEdr == true) ||
                (isComputeVirial == true) ||
                (isComputeParticleVirial == true))
            {
              double const rij = sqrt(rij2);
              double const dEidr = dEdrByRij*rij;

              if (isComputeProcess_dEdr == true)
              {
                ier = modelComputeArguments
                    ->ProcessDEDrTerm(dEidr, rij, r_ij, i, j);
                if (ier)
                {
                  LOG_ERROR("process_dEdr");
                  return ier;
                }
              }

              if (isComputeVirial == true)
              {
                ProcessVirialTerm(dEidr, rij, r_ij, i, j, virial);
              }

              if (isComputeParticleVirial == true)
              {
                ProcessParticleVirialTerm(dEidr, rij, r_ij, i, j,
                                          particleVirial);
              }
            }
          }  // if particles i and j interact
        } // end effective half-list check (i < j)
      } // end of first neighbor loop
    }  // end of if statement for whether particle is contributing
  }  // end of loop over contributing particles

  // Separate loop nest for process_d2Edr2
  if (isComputeProcess_d2Edr2 == true)
  {
    // For this potential, the second derivative
    // d^2E/dr_{ij}dr_{kl} is nonzero only if i,j,k,l are
    // not distinct. As such, we need only address all
    // derivatives d^2E/dr_{ij}dr_{ik}.
    //
    // If we have a full list, we can simply iterate over
    // particle i's neighbor list in a doubly-nested
    // triangular loop.
    //
    // If we have half lists, then the presence of j in
    // i's neighbor list implies that i is also a neighbor
    // of j, so that we should be processing all derivs
    // d^2E/dr_{ij}dr_{ik} and d^2E/dr_{ji}dr_{jl}. The
    // hurdle to overcome here is the need to access the
    // neighbor list of particle j (and also particle k).
    //
    // We avoid this problem by taking a different
    // approach. For each neighbor (j or k) of i, we save
    // i in a "deferred neighbor list" for each of j and
    // k. For each particle, we process neighbors both in
    // its neighbor list and its deferred neighbor list.
    // (In effect, we are reconstructing full lists for
    // each particle.)

    // Containers for deferred neighbor lists for process_d2Edr2
    std::multimap<int, neighbor> deferredNeighborMap;

    // Iterators over deferred neighbor lists
    std::pair<deferredNeighborIterator, deferredNeighborIterator>
        deferredNeighborRange;
    deferredNeighborIterator deferredNeighborPtrJ;

    // First, do a check to determine whether we've been given a half list or
    // a full list
    bool isHalf = false;
    for (i = 0; i < cachedNumberOfParticles_; ++i)
    {
      modelComputeArguments->GetNeighborList(0, i, &numnei, &n1atom);

      if (numnei > 0)
      {
        // Get the neighbor list of the first neighbor of atom i and check to
        // see if atom i is in its list
        modelComputeArguments->GetNeighborList(0, 0, &numnei, &n1atom);
        if (numnei > 0)
        {
          isHalf = true;
          for (int kk = 0; kk < numnei; ++kk)
          {
            if (n1atom[kk] == i)
            {
              isHalf = false;
              break;
            }
          }
        }
        else
        {
          isHalf = false;
        }
        break; // Stop loop over i
      }
    }

    // Setup loop over contributing particles
    for (i = 0; i < cachedNumberOfParticles_; ++i)
    {
      if (particleContributing[i])
      {
        modelComputeArguments->GetNeighborList(0, i, &numnei, &n1atom);

        // Number of neighbors to visit. May be greater than
        // numnei if previous particles inserted deferred
        // neighbors to this particle, or if this particle
        // is its own neighbor.
        int extendedNumnei;
        extendedNumnei = numnei;

        // Begin by deferring i to each of its neighbors.
        // It is tempting to try to include this in the
        // subsequent loop that also does the calculations,
        // but this cannot be done. The complete deferred
        // list must be available the first time we go
        // through the inner, second-neighbor, loop.
        if (isHalf == true)
        {
          for (int jj = 0; jj < numnei; ++jj)
          {
            neighbor reflectedNeighJ;

            // Index of particle neighbor
            int const j = n1atom[jj];

            if (particleContributing[j])
            {
              reflectedNeighJ.index = i;

              deferredNeighborMap.insert(std::make_pair(j, reflectedNeighJ));
            }
          }
        }

        // No further particles will contribute to i's neighbor list.
        // Safe to finalize the neighbor count.
        extendedNumnei += deferredNeighborMap.count(i);
        deferredNeighborRange = deferredNeighborMap.equal_range(i);
        deferredNeighborPtrJ  = deferredNeighborRange.first;

        // Setup loop over neighbors of current particle
        for (int jj = 0; jj < extendedNumnei; ++jj)
        {
          int j;

          if (isHalf && jj >= numnei)
          { // Look in deferred neighbor list
            j = deferredNeighborPtrJ->second.index;
          }
          else
          { // Not a deferred neighbor
            j = n1atom[jj];
          }

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
            INTERPOLATE_F(rijPhiAlphaBetaCoeff,
                          rijOffset, rijIndex, rijPhiValue);

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
            double const* const densityBetaCoeff =
                densityCoeff_[particleSpeciesCodes[j]][particleSpeciesCodes[i]];
            INTERPOLATE_DF(densityBetaCoeff, rijOffset, rijIndex,
                           densityBetaDerivativeValue);

            // Setup second loop over neighbors of current particle
            deferredNeighborIterator
                  deferredNeighborPtrK = deferredNeighborPtrJ;
            for (int kk = jj; kk < extendedNumnei; ++kk)
            {
              int k;
              if (isHalf && kk >= numnei)
              { // Look in deferred neighbor list
                k = deferredNeighborPtrK->second.index;
              }
              else
              {
                k = n1atom[kk];
              }

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
                GET_DELTAX_AND_INDEX(
                    rij[1], oneByDr_, numberRPoints_, rikOffset, rikIndex);

                // interpolate derivative of rho_gamma(r_ik)
                double densityGammaDerivativeValue;
                double const* const densityGammaCoeff = densityCoeff_[
                    particleSpeciesCodes[k]][particleSpeciesCodes[i]];
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
                  double const* const densityBetaCoeff = densityCoeff_[
                      particleSpeciesCodes[j]][particleSpeciesCodes[i]];
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

                  ier = modelComputeArguments->ProcessD2EDr2Term(
                      d2Edr2, prikrij, pr_ikr_ij, piis, pjjs);
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

              if (isHalf && kk >= numnei)
              {
                ++deferredNeighborPtrK;
              }
            }  // end of second neighbor loop
          }  // if particles i and j interact

          if (isHalf && jj >= numnei)
          {
            ++deferredNeighborPtrJ;
          }
        }  // end of first neighbor loop

        deferredNeighborMap.erase(i);
      }  // end of if statement for whether particle is contributing
    }  // end of loop over contributing particles
  }  // if (isComputeProcess_d2Edr2 == true)

  // everything is good
  ier = false;
  return ier;
}

#endif  // EAM_IMPLEMENTATION_HPP_
