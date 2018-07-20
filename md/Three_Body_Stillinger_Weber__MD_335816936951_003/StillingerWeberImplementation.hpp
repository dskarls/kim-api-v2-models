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


#ifndef STILLINGER_WEBER_IMPLEMENTATION_HPP_
#define STILLINGER_WEBER_IMPLEMENTATION_HPP_

#include <vector>
#include "KIM_LogVerbosity.hpp"
#include "StillingerWeber.hpp"
#include "helper.hpp"

#define DIMENSION 3
#define ONE 1.0
#define HALF 0.5

#define MAX_PARAMETER_FILES 1


//==============================================================================
//
// Declaration of StillingerWeberImplementation class
//
//==============================================================================

//******************************************************************************
class StillingerWeberImplementation
{
public:
  StillingerWeberImplementation(
      KIM::ModelDriverCreate* const modelDriverCreate,
      KIM::LengthUnit const requestedLengthUnit,
      KIM::EnergyUnit const requestedEnergyUnit,
      KIM::ChargeUnit const requestedChargeUnit,
      KIM::TemperatureUnit const requestedTemperatureUnit,
      KIM::TimeUnit const requestedTimeUnit,
      int* const ier);
  ~StillingerWeberImplementation();  // no explicit Destroy() needed here

  int Refresh(KIM::ModelRefresh* const modelRefresh);
  int Compute(KIM::ModelCompute const* const modelCompute,
      KIM::ModelComputeArguments const* const modelComputeArguments);
  int ComputeArgumentsCreate(
      KIM::ModelComputeArgumentsCreate* const modelComputeArgumentsCreate) const;
  int ComputeArgumentsDestroy(
      KIM::ModelComputeArgumentsDestroy* const modelComputeArgumentsDestroy) const;


private:
  // Constant values that never change
  //   Set in constructor (via SetConstantValues)
  //
  //
  // StillingerWeberImplementation: constants
  int numberModelSpecies_;
  std::vector<int> modelSpeciesCodeList_;
  int numberUniqueSpeciesPairs_;


  // Constant values that are read from the input files and never change
  //   Set in constructor (via functions listed below)
  //
  //
  // Private Model Parameters
  //   Memory allocated in AllocatePrivateParameterMemory() (from constructor)
  //   Memory deallocated in destructor
  //   Data set in ReadParameterFile routines
  // none
  //
  // KIM API: Model Parameters whose (pointer) values never change
  //   Memory allocated in AllocateParameterMemory() (from constructor)
  //   Memory deallocated in destructor
  //   Data set in ReadParameterFile routines OR by KIM Simulator
  double* cutoff_;
  double* A_;
  double* B_;
  double* p_;
  double* q_;
  double* sigma_;
  double* lambda_;
  double* gamma_;
  double* costheta0_;


  // Mutable values that only change when Refresh() executes
  //   Set in Refresh (via SetRefreshMutableValues)
  //
  //
  // KIM API: Model Parameters (can be changed directly by KIM Simulator)
  // none
  //
  // StillingerWeberImplementation: values (changed only by Refresh())
  double influenceDistance_;
  int paddingNeighborHints_;
  int halfListHints_;

  double** cutoffSq_2D_;
  double** A_2D_;
  double** B_2D_;
  double** p_2D_;
  double** q_2D_;
  double** sigma_2D_;
  double** lambda_2D_;
  double** gamma_2D_;
  double** costheta0_2D_;


  // Mutable values that can change with each call to Refresh() and Compute()
  //   Memory may be reallocated on each call
  //
  //
  // StillingerWeberImplementation: values that change
  int cachedNumberOfParticles_;


  // Helper methods
  //
  //
  // Related to constructor
  void AllocatePrivateParameterMemory();
  void AllocateParameterMemory();

  static int OpenParameterFiles(
      KIM::ModelDriverCreate * const modelDriverCreate,
      int const numberParameterFiles,
      FILE * parameterFilePointers[MAX_PARAMETER_FILES]);
  int ProcessParameterFiles(
      KIM::ModelDriverCreate* const modelDriverCreate,
      int const numberParameterFiles,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES]);
  void getNextDataLine(
      FILE* const filePtr, char* const nextLine,
      int const maxSize, int* endOfFileFlag);
  static void CloseParameterFiles(
      int const numberParameterFiles,
      FILE* const parameterFilePointers[MAX_PARAMETER_FILES]);
  int ConvertUnits(
      KIM::ModelDriverCreate* const modelDriverCreate,
      KIM::LengthUnit const requestedLengthUnit,
      KIM::EnergyUnit const requestedEnergyUnit,
      KIM::ChargeUnit const requestedChargeUnit,
      KIM::TemperatureUnit const requestedTemperatureUnit,
      KIM::TimeUnit const requestedTimeUnit);
  int RegisterKIMModelSettings(
      KIM::ModelDriverCreate* const modelDriverCreate) const;
  int RegisterKIMComputeArgumentsSettings(
      KIM::ModelComputeArgumentsCreate* const modelComputeArgumentsCreate) const;
  int RegisterKIMParameters(KIM::ModelDriverCreate* const modelDriverCreate);
  int RegisterKIMFunctions(KIM::ModelDriverCreate* const modelDriverCreate) const;

  //
  // Related to Refresh()
  template<class ModelObj>
  int SetRefreshMutableValues(ModelObj* const modelObj);

  //
  // Related to Compute()
  int SetComputeMutableValues(
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
      VectorOfSizeSix*& particleViral);
  int CheckParticleSpeciesCodes(
      KIM::ModelCompute const* const modelCompute,
      int const* const particleSpeciesCodes) const;
  int GetComputeIndex(
      const bool& isComputeProcess_dEdr,
      const bool& isComputeProcess_d2Edr2,
      const bool& isComputeEnergy,
      const bool& isComputeForces,
      const bool& isComputeParticleEnergy,
      const bool& isComputeVirial,
      const bool& isComputeParticleVirial) const;

  // compute functions
  template<bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
      bool isComputeEnergy, bool isComputeForces,
      bool isComputeParticleEnergy, bool isComputeVirial,
      bool isComputeParticleVirial>
  int Compute(
      KIM::ModelCompute const* const modelCompute,
      KIM::ModelComputeArguments const* const modelComputeArguments,
      const int* const particleSpeciesCodes,
      const int* const particleContributing,
      const VectorOfSizeDIM* const coordinates,
      double* const energy,
      VectorOfSizeDIM* const forces,
      double* const particleEnergy,
      VectorOfSizeSix virial,
      VectorOfSizeSix* const particleVirial) const;


  // Stillinger-Weber functions
  void CalcPhiTwo(int const ispec, int const jspec, double const r, double& phi) const;
  void CalcPhiDphiTwo(int const ispec, int const jspec, double const r,
      double& phi, double& dphi) const;
  void CalcPhiD2phiTwo(int const ispec, int const jspec, double const r,
      double& phi, double& dphi, double& d2phi) const;
  void CalcPhiThree(int const ispec, int const jspec, int const kspec,
      double const rij, double const rik, double const rjk,
      double& phi) const;
  void CalcPhiDphiThree(int const ispec, int const jspec, int const kspec,
      double const rij, double const rik, double const rjk,
      double& phi, double* const dphi) const;
  void CalcPhiD2phiThree(int const ispec, int const jspec, int const kspec,
      double const rij, double const rik, double const rjk,
      double& phi, double* const dphi, double* const d2phi) const;
};

//==============================================================================
//
// Definition of StillingerWeberImplementation::Compute functions
//
// NOTE: Here we rely on the compiler optimizations to prune dead code
//       after the template expansions.  This provides high efficiency
//       and easy maintenance.
//
//==============================================================================

#include "KIM_ModelComputeLogMacros.hpp"
template<bool isComputeProcess_dEdr, bool isComputeProcess_d2Edr2,
    bool isComputeEnergy, bool isComputeForces,
    bool isComputeParticleEnergy, bool isComputeVirial,
    bool isComputeParticleVirial>
int StillingerWeberImplementation::Compute(
    KIM::ModelCompute const* const modelCompute,
    KIM::ModelComputeArguments const* const modelComputeArguments,
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
      (isComputeParticleVirial == false)) {
    return ier;
  }

  // initialize energy and forces
  if (isComputeEnergy == true) {
    *energy = 0.0;
  }

  if (isComputeForces == true) {
    for (int i = 0; i < cachedNumberOfParticles_; ++i) {
      for (int j = 0; j < DIMENSION; ++j) {
        forces[i][j] = 0.0;
      }
    }
  }

  if (isComputeParticleEnergy == true) {
    for (int i = 0; i < cachedNumberOfParticles_; ++i) {
      particleEnergy[i] = 0.0;
    }
  }

  if (isComputeVirial == true) {
    for (int i = 0; i < 6; ++i) {
      virial[i] = 0.0;
    }
  }

  if (isComputeParticleVirial == true) {
    for (int i = 0; i < cachedNumberOfParticles_; ++i) {
      for (int j = 0; j < 6; ++j) {
        particleVirial[i][j] = 0.0;
      }
    }
  }


  // calculate contribution from pair function
  //
  // Setup loop over contributing particles
  int i = 0;
  int numnei = 0;
  int const* n1atom = NULL;

  for (i = 0; i < cachedNumberOfParticles_; ++i) {

    if (particleContributing[i]) {
      modelComputeArguments->GetNeighborList(0, i, &numnei, &n1atom);
      int const iSpecies = particleSpeciesCodes[i];

      // Setup loop over neighbors of current particle
      for (int jj = 0; jj < numnei; ++jj) {
        int const j = n1atom[jj];
        int const jSpecies = particleSpeciesCodes[j];

        // Compute rij
        double rij[DIMENSION];
        for (int dim = 0; dim < DIMENSION; ++dim) {
          rij[dim] = coordinates[j][dim] - coordinates[i][dim];
        }

        // compute distance squared
        double const rij_sq = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        if (rij_sq <= cutoffSq_2D_[iSpecies][jSpecies]) {
          double const rij_mag = sqrt(rij_sq);

          // two-body contributions

          if (i < j) {  // effective half list
            double phi_two = 0.0;
            double dphi_two = 0.0;
            double d2phi_two = 0.0;
            double dEidr_two = 0.0;
            double d2Eidr2_two = 0.0;

            // Compute two body potenitals and its derivatives
            if (isComputeProcess_d2Edr2 == true) {
              CalcPhiD2phiTwo(iSpecies, jSpecies, rij_mag, phi_two, dphi_two, d2phi_two);
              if (particleContributing[j] == 1) {
                dEidr_two = dphi_two;
                d2Eidr2_two = d2phi_two;
              }
              else {
                dEidr_two = HALF * dphi_two;
                d2Eidr2_two = HALF * d2phi_two;
              }
            }
            else if ((isComputeProcess_dEdr == true) || (isComputeForces == true) ||
                     (isComputeVirial == true) || (isComputeParticleVirial == true)) {
              CalcPhiDphiTwo(iSpecies, jSpecies, rij_mag, phi_two, dphi_two);
              if (particleContributing[j] == 1) {
                dEidr_two = dphi_two;
              }
              else {
                dEidr_two = HALF * dphi_two;
              }
            }
            else if ((isComputeEnergy == true) || (isComputeParticleEnergy == true)) {
              CalcPhiTwo(iSpecies, jSpecies, rij_mag, phi_two);
            }

            // Contribution to energy
            if (isComputeEnergy == true) {
              if (particleContributing[j] == 1) {
                *energy += phi_two;
              }
              else {
                *energy += HALF * phi_two;
              }
            }

            // Contribution to forces
            if (isComputeForces == true) {
              for (int dim = 0; dim < DIMENSION; ++dim) {
                double const contrib = dEidr_two * rij[dim] / rij_mag;
                forces[i][dim] += contrib;
                forces[j][dim] -= contrib;
              }
            }

            // Contribution to particleEnergy
            if (isComputeParticleEnergy == true) {
              double halfphi = HALF * phi_two;
              particleEnergy[i] += halfphi;
              if (particleContributing[j] == 1) {
                particleEnergy[j] += halfphi;
              }
            }

            // Contribution to virial
            if (isComputeVirial == true) {
              ProcessVirialTerm(dEidr_two, rij_mag, rij, i, j, virial);
            }

            // Contribution to particleVirial
            if (isComputeParticleVirial == true) {
              ProcessParticleVirialTerm(dEidr_two, rij_mag, rij, i, j, particleVirial);
            }

            // Call process_dEdr
            if (isComputeProcess_dEdr == true) {
              ier = modelComputeArguments->ProcessDEDrTerm(dEidr_two, rij_mag, rij, i, j);
              if (ier) {
                LOG_ERROR("ProcessDEdr");
                return ier;
              }
            }

            // Call process_d2Edr2
            if (isComputeProcess_d2Edr2 == true) {
              double const R_pairs[2] = { rij_mag, rij_mag };
              double const* const pRs = &R_pairs[0];
              double const Rij_pairs[6]
                = { rij[0], rij[1], rij[2],
                    rij[0], rij[1], rij[2] };
              double const* const pRijConsts = &Rij_pairs[0];
              int const i_pairs[2] = { i, i };
              int const j_pairs[2] = { j, j };
              int const* const pis = &i_pairs[0];
              int const* const pjs = &j_pairs[0];

              ier = modelComputeArguments->ProcessD2EDr2Term(d2Eidr2_two, pRs,
                  pRijConsts, pis, pjs);
              if (ier) {
                LOG_ERROR("ProcessD2Edr2");
                return ier;
              }
            }
          }  // i < j


          // three-body contribution
          for (int kk = jj + 1; kk < numnei; ++kk) {
            int const k = n1atom[kk];
            int const kSpecies = particleSpeciesCodes[k];

            // Compute rik and rjk vector
            double rik[DIMENSION];
            for (int dim = 0; dim < DIMENSION; ++dim) {
              rik[dim] = coordinates[k][dim] - coordinates[i][dim];
            }

            // compute distance squared and distance
            double const rik_sq = rik[0] * rik[0] + rik[1] * rik[1] + rik[2] * rik[2];


            // compute energy and force
            if (rik_sq <= cutoffSq_2D_[iSpecies][kSpecies]) {
              double const rik_mag = sqrt(rik_sq);

              // Compute rjk
              double rjk[DIMENSION];
              for (int dim = 0; dim < DIMENSION; ++dim) {
                rjk[dim] = coordinates[k][dim] - coordinates[j][dim];
              }
              double const rjk_sq = rjk[0] * rjk[0] + rjk[1] * rjk[1] + rjk[2] * rjk[2];
              double const rjk_mag = sqrt(rjk_sq);



              // three-body contributions
              double phi_three;
              double dphi_three[3];
              double d2phi_three[6];
              double dEidr_three[3];
              double d2Eidr2_three[6];

              // compute three-body potential and its derivatives
              if (isComputeProcess_d2Edr2 == true) {
                CalcPhiD2phiThree(iSpecies, jSpecies, kSpecies,
                    rij_mag, rik_mag, rjk_mag, phi_three, dphi_three, d2phi_three);

                dEidr_three[0] = dphi_three[0];
                dEidr_three[1] = dphi_three[1];
                dEidr_three[2] = dphi_three[2];

                d2Eidr2_three[0] = d2phi_three[0];
                d2Eidr2_three[1] = d2phi_three[1];
                d2Eidr2_three[2] = d2phi_three[2];
                d2Eidr2_three[3] = d2phi_three[3];
                d2Eidr2_three[4] = d2phi_three[4];
                d2Eidr2_three[5] = d2phi_three[5];
              }
              else if ((isComputeProcess_dEdr == true) || (isComputeForces == true) ||
                       (isComputeVirial == true) || (isComputeParticleVirial == true)) {
                CalcPhiDphiThree(iSpecies, jSpecies, kSpecies,
                    rij_mag, rik_mag, rjk_mag, phi_three, dphi_three);

                dEidr_three[0] = dphi_three[0];
                dEidr_three[1] = dphi_three[1];
                dEidr_three[2] = dphi_three[2];
              }
              else if ((isComputeEnergy == true) || (isComputeParticleEnergy == true)) {
                CalcPhiThree(iSpecies, jSpecies, kSpecies,
                    rij_mag, rik_mag, rjk_mag, phi_three);
              }

              // Contribution to energy
              if (isComputeEnergy == true) {
                *energy += phi_three;
              }

              // Contribution to forces
              if (isComputeForces == true) {
                for (int dim = 0; dim < DIMENSION; ++dim) {
                  double const contrib0 = dEidr_three[0] * rij[dim] / rij_mag;
                  double const contrib1 = dEidr_three[1] * rik[dim] / rik_mag;
                  double const contrib2 = dEidr_three[2] * rjk[dim] / rjk_mag;
                  forces[i][dim] += contrib0 + contrib1;
                  forces[j][dim] += -contrib0 + contrib2;
                  forces[k][dim] += -contrib2 - contrib1;
                }
              }

              // Contribution to particleEnergy
              if (isComputeParticleEnergy == true) {
                particleEnergy[i] += phi_three;
              }

              // Contribution to virial
              if (isComputeVirial == true) {
                ProcessVirialTerm(dEidr_three[0], rij_mag, rij, i, j, virial);
                ProcessVirialTerm(dEidr_three[1], rik_mag, rik, i, k, virial);
                ProcessVirialTerm(dEidr_three[2], rjk_mag, rjk, j, k, virial);
              }

              // Contribution to particleVirial
              if (isComputeParticleVirial == true) {
                ProcessParticleVirialTerm(dEidr_three[0], rij_mag, rij, i, j, particleVirial);
                ProcessParticleVirialTerm(dEidr_three[1], rik_mag, rik, i, k, particleVirial);
                ProcessParticleVirialTerm(dEidr_three[2], rjk_mag, rjk, j, k, particleVirial);
              }

              // Call process_dEdr
              if (isComputeProcess_dEdr == true) {
                ier =
                  modelComputeArguments->ProcessDEDrTerm(
                      dEidr_three[0], rij_mag, rij, i, j) ||
                  modelComputeArguments->ProcessDEDrTerm(
                      dEidr_three[1], rik_mag, rik, i, k) ||
                  modelComputeArguments->ProcessDEDrTerm(
                      dEidr_three[2], rjk_mag, rjk, j, k);
                if (ier) {
                  LOG_ERROR("ProcessDEdr");
                  return ier;
                }
              }

              // Call process_d2Edr2
              if (isComputeProcess_d2Edr2 == true) {
                double R_pairs[2];
                double Rij_pairs[6];
                int i_pairs[2];
                int j_pairs[2];
                double* const pRs = &R_pairs[0];
                double* const pRijConsts = &Rij_pairs[0];
                int* const pis = &i_pairs[0];
                int* const pjs = &j_pairs[0];

                R_pairs[0] = R_pairs[1] = rij_mag;
                Rij_pairs[0] = Rij_pairs[3] = rij[0];
                Rij_pairs[1] = Rij_pairs[4] = rij[1];
                Rij_pairs[2] = Rij_pairs[5] = rij[2];
                i_pairs[0] = i_pairs[1] = i;
                j_pairs[0] = j_pairs[1] = j;
                ier = modelComputeArguments
                      ->ProcessD2EDr2Term(d2Eidr2_three[0], pRs, pRijConsts, pis, pjs);
                if (ier) {
                  LOG_ERROR("ProcessD2Edr2");
                  return ier;
                }

                R_pairs[0] = R_pairs[1] = rik_mag;
                Rij_pairs[0] = Rij_pairs[3] = rik[0];
                Rij_pairs[1] = Rij_pairs[4] = rik[1];
                Rij_pairs[2] = Rij_pairs[5] = rik[2];
                i_pairs[0] = i_pairs[1] = i;
                j_pairs[0] = j_pairs[1] = k;
                ier = modelComputeArguments
                      ->ProcessD2EDr2Term(d2Eidr2_three[1], pRs, pRijConsts, pis, pjs);
                if (ier) {
                  LOG_ERROR("ProcessD2Edr2");
                  return ier;
                }

                R_pairs[0] = R_pairs[1] = rjk_mag;
                Rij_pairs[0] = Rij_pairs[3] = rjk[0];
                Rij_pairs[1] = Rij_pairs[4] = rjk[1];
                Rij_pairs[2] = Rij_pairs[5] = rjk[2];
                i_pairs[0] = i_pairs[1] = j;
                j_pairs[0] = j_pairs[1] = k;
                ier = modelComputeArguments
                      ->ProcessD2EDr2Term(d2Eidr2_three[2], pRs, pRijConsts, pis, pjs);
                if (ier) {
                  LOG_ERROR("ProcessD2Edr2");
                  return ier;
                }

                R_pairs[0] = rij_mag;
                R_pairs[1] = rik_mag;
                Rij_pairs[0] = rij[0];
                Rij_pairs[1] = rij[1];
                Rij_pairs[2] = rij[2];
                Rij_pairs[3] = rik[0];
                Rij_pairs[4] = rik[1];
                Rij_pairs[5] = rik[2];
                i_pairs[0] = i;
                j_pairs[0] = j;
                i_pairs[1] = i;
                j_pairs[1] = k;
                ier = modelComputeArguments
                      ->ProcessD2EDr2Term(d2Eidr2_three[3], pRs, pRijConsts, pis, pjs);
                if (ier) {
                  LOG_ERROR("ProcessD2Edr2");
                  return ier;
                }

                R_pairs[0] = rik_mag;
                R_pairs[1] = rij_mag;
                Rij_pairs[0] = rik[0];
                Rij_pairs[1] = rik[1];
                Rij_pairs[2] = rik[2];
                Rij_pairs[3] = rij[0];
                Rij_pairs[4] = rij[1];
                Rij_pairs[5] = rij[2];
                i_pairs[0] = i;
                j_pairs[0] = k;
                i_pairs[1] = i;
                j_pairs[1] = j;
                ier = modelComputeArguments
                      ->ProcessD2EDr2Term(d2Eidr2_three[3], pRs, pRijConsts, pis, pjs);
                if (ier) {
                  LOG_ERROR("ProcessD2Edr2");
                  return ier;
                }

                R_pairs[0] = rij_mag;
                R_pairs[1] = rjk_mag;
                Rij_pairs[0] = rij[0];
                Rij_pairs[1] = rij[1];
                Rij_pairs[2] = rij[2];
                Rij_pairs[3] = rjk[0];
                Rij_pairs[4] = rjk[1];
                Rij_pairs[5] = rjk[2];
                i_pairs[0] = i;
                j_pairs[0] = j;
                i_pairs[1] = j;
                j_pairs[1] = k;
                ier = modelComputeArguments
                      ->ProcessD2EDr2Term(d2Eidr2_three[4], pRs, pRijConsts, pis, pjs);
                if (ier) {
                  LOG_ERROR("ProcessD2Edr2");
                  return ier;
                }


                R_pairs[0] = rjk_mag;
                R_pairs[1] = rij_mag;
                Rij_pairs[0] = rjk[0];
                Rij_pairs[1] = rjk[1];
                Rij_pairs[2] = rjk[2];
                Rij_pairs[3] = rij[0];
                Rij_pairs[4] = rij[1];
                Rij_pairs[5] = rij[2];
                i_pairs[0] = j;
                j_pairs[0] = k;
                i_pairs[1] = i;
                j_pairs[1] = j;
                ier = modelComputeArguments
                      ->ProcessD2EDr2Term(d2Eidr2_three[4], pRs, pRijConsts, pis, pjs);
                if (ier) {
                  LOG_ERROR("ProcessD2Edr2");
                  return ier;
                }

                R_pairs[0] = rik_mag;
                R_pairs[1] = rjk_mag;
                Rij_pairs[0] = rik[0];
                Rij_pairs[1] = rik[1];
                Rij_pairs[2] = rik[2];
                Rij_pairs[3] = rjk[0];
                Rij_pairs[4] = rjk[1];
                Rij_pairs[5] = rjk[2];
                i_pairs[0] = i;
                j_pairs[0] = k;
                i_pairs[1] = j;
                j_pairs[1] = k;
                ier = modelComputeArguments
                      ->ProcessD2EDr2Term(d2Eidr2_three[5], pRs, pRijConsts, pis, pjs);
                if (ier) {
                  LOG_ERROR("ProcessD2Edr2");
                  return ier;
                }

                R_pairs[0] = rjk_mag;
                R_pairs[1] = rik_mag;
                Rij_pairs[0] = rjk[0];
                Rij_pairs[1] = rjk[1];
                Rij_pairs[2] = rjk[2];
                Rij_pairs[3] = rik[0];
                Rij_pairs[4] = rik[1];
                Rij_pairs[5] = rik[2];
                i_pairs[0] = j;
                j_pairs[0] = k;
                i_pairs[1] = i;
                j_pairs[1] = k;
                ier = modelComputeArguments
                      ->ProcessD2EDr2Term(d2Eidr2_three[5], pRs, pRijConsts, pis, pjs);
                if (ier) {
                  LOG_ERROR("ProcessD2Edr2");
                  return ier;
                }
              } // Process_D2Edr2
            }   // if particleContributing
          }     // if particles i and k interact
        }       // if particles i and j interact
      }         // end of first neighbor loop
    }           // if particleContributing
  }             // loop over all particles

  // everything is good
  ier = false;
  return ier;
}


#endif  // STILLINGER_WEBER_IMPLEMENTATION_HPP_
