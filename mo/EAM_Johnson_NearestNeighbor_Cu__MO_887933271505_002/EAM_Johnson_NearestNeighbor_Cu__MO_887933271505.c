/*
*
* CDDL HEADER START
*
* The contents of this file are subject to the terms of the Common Development
* and Distribution License Version 1.0 (the "License").
*
* You can obtain a copy of the license at
* http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
* specific language governing permissions and limitations under the License.
*
* When distributing Covered Code, include this CDDL HEADER in each file and
* include the License file in a prominent location with the name LICENSE.CDDL.
* If applicable, add the following below this CDDL HEADER, with the fields
* enclosed by brackets "[]" replaced with your own identifying information:
*
* Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
*
* CDDL HEADER END
*

*
* Copyright (c) 2013, Regents of the University of Minnesota.
* All rights reserved.
*
* Contributors:
*    Ryan S. Elliott
*    Ellad B. Tadmor
*    Stephen M. Whalen
*
*/

/*******************************************************************************
*
*  EAM_Johnson_NearestNeighbor_Cu
*
*  Johnson pair functional model for Cu
*
*  Reference: R. A. Johnson, "Analytic nearest-neighbor model for fcc metals",
*             Phys. Rev. B, 55(8):4941-4946, 1988.
*
*  Language: C
*
*******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "KIM_ModelHeaders.h"

#define TRUE 1
#define FALSE 0

/******************************************************************************
* Below are the definitions and values of all Model parameters
*******************************************************************************/
#define DIM 3       /* dimensionality of space */
#define SPECCODE 1  /* internal species code */
#define MODEL_CUTOFF 3.5 /* cutoff radius in angstroms */
#define MODEL_CUTSQ  (MODEL_CUTOFF * MODEL_CUTOFF)
#define JEAM_R0  2.556 /* A */
#define JEAM_PHI0 0.59 /* eV */
#define JEAM_GAM  8.00
#define JEAM_G0   0.30 /* eV */
#define JEAM_BET  5.85
#define JEAM_EC   3.54 /* eV/atom */
#define JEAM_ALF  5.09
#define JEAM_RHO0 3.60 /* eV (=12*JEAM_G0) */

/* Model buffer definition */
struct buffer
{
  int paddingNeighborHints;
  int halfListHints;

  double influenceDistance;
  double cutoff;
};
typedef struct buffer buffer;

/* Define prototype for Model create */
int model_create(KIM_ModelCreate * const modelCreate,
                 KIM_LengthUnit const requestedLengthUnit,
                 KIM_EnergyUnit const requestedEnergyUnit,
                 KIM_ChargeUnit const requestedChargeUnit,
                 KIM_TemperatureUnit const requestedTemperatureUnit,
                 KIM_TimeUnit const requestedTimeUnit);

/* Define prototypes for model reinit, compute, and destroy */
/* defined as static to avoid namespace clashes with other Models */
/**/
static int compute(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArguments const * const modelComputeArguments);
static int compute_arguments_create(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsCreate * const modelComputeArgumentsCreate);
static int compute_arguments_destroy(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy);

static int model_refresh(KIM_ModelRefresh * const modelRefresh);
static int model_destroy(KIM_ModelDestroy * const modelDestroy);
/**/
static void calc_phi(double r, double* phi);
static void calc_phi_dphi(double r, double* phi, double* dphi);
static void calc_g(double r, double* g);
static void calc_dg(double r, double* dg);
static void calc_U(double rho, double* U);
static void calc_U_dU(double rho, double* U, double* dU);

/* Calculate pair potential phi(r) */
static void calc_phi(double r, double* phi)
{
   /* local variables */
   double rnorm;

   if (r > MODEL_CUTOFF)
   {
      /* Argument exceeds cutoff radius */
      *phi = 0.0;
   }
   else
   {
      rnorm = r/JEAM_R0;
      *phi = JEAM_PHI0 * exp(-JEAM_GAM*(rnorm - 1.0));
   }

   return;
}

/* Calculate pair potential phi(r) and its derivative dphi(r) */
static void calc_phi_dphi(double r, double* phi, double* dphi)
{
   /* local variables */
   double rnorm;

   if (r > MODEL_CUTOFF)
   {
      /* Argument exceeds cutoff radius */
      *phi  = 0.0;
      *dphi = 0.0;
   }
   else
   {
      rnorm = r/JEAM_R0;
      *phi  = JEAM_PHI0 * exp(-JEAM_GAM*(rnorm - 1.0));
      *dphi = -(JEAM_GAM/JEAM_R0)*(*phi);
   }

   return;
}

/* Calculate electron density g(r) */
static void calc_g(double r, double* g)
{
   /* local variables */
   double rnorm;

   if (r > MODEL_CUTOFF)
   {
      /* Argument exceeds cutoff radius */
      *g = 0.0;
   }
   else
   {
      rnorm = r/JEAM_R0;
      *g = JEAM_G0 * exp(-JEAM_BET*(rnorm - 1.0));
   }

   return;
}

/* Calculate electron density derivative dg(r) */
static void calc_dg(double r, double* dg)
{
   /* local variables */
   double rnorm;
   double g;

   if (r > MODEL_CUTOFF)
   {
      /* Argument exceeds cutoff radius */
      *dg = 0.0;
   }
   else
   {
      rnorm = r/JEAM_R0;
      g = JEAM_G0 * exp(-JEAM_BET*(rnorm - 1.0));
      *dg = -(JEAM_BET/JEAM_R0)*g;
   }

   return;
}

/* Calculate embedding function U(rho) */
static void calc_U(double rho, double* U)
{
   /* local variables */
   double rhonorm;
   double rhonorm_gob;
   double rhonorm_aob;
   double logrhonorm;
   double aob;
   double gob;

   if (rho == 0.0)
   {
      *U = 0.0;

      return;
   }

   aob = JEAM_ALF/JEAM_BET;
   gob = JEAM_GAM/JEAM_BET;
   rhonorm = rho/JEAM_RHO0;
   rhonorm_aob = pow(rhonorm, aob);
   rhonorm_gob = pow(rhonorm, gob);
   logrhonorm = log(rhonorm);

   *U = -JEAM_EC * (1.0 - aob*logrhonorm)*rhonorm_aob - 6.0*JEAM_PHI0*rhonorm_gob;

   return;
}

/* Calculate embedding function U(rho) and first derivative dU(rho) */
static void calc_U_dU(double rho, double* U, double* dU)
{
   /* local variables */
   double rhonorm;
   double rhonorm_gob;
   double rhonorm_aob;
   double logrhonorm;
   double aob;
   double gob;

   if (rho == 0.0)
   {
      *U = 0.0;
      *dU = 0.0;

      return;
   }

   aob = JEAM_ALF/JEAM_BET;
   gob = JEAM_GAM/JEAM_BET;
   rhonorm = rho/JEAM_RHO0;
   rhonorm_aob = pow(rhonorm, aob);
   rhonorm_gob = pow(rhonorm, gob);
   logrhonorm = log(rhonorm);

   *U = -JEAM_EC * (1.0 - aob*logrhonorm)*rhonorm_aob - 6.0*JEAM_PHI0*rhonorm_gob;
   *dU = (JEAM_EC*aob*aob*logrhonorm*rhonorm_aob - 6.0*JEAM_PHI0*gob*rhonorm_gob)/rho;

   return;
}

/* compute function */
#include "KIM_ModelComputeLogMacros.h"
static int compute(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArguments const * const modelComputeArguments)
{
   double r;
   double Rij[DIM];
   double Rsqij;
   double phi;
   double dphi;
   double g;
   double dg;
   double dU;
   double dphieff = 0.0;
   double dphii;
   double dUi;
   double Ei;
   double dphij;
   double dUj;
   double Ej;
   int ier;
   int i;
   int j;
   int jj;
   int k;
   int comp_force;
   int comp_particleEnergy;
   int comp_virial;
   int comp_energy;
   double* rho;
   double U;
   double* derU = 0;

   int const * neighListOfCurrentPart;
   int numOfPartNeigh;
   int* nParts;
   double* energy;
   double* coords;
   double* force;
   double* particleEnergy;

   int* particleContributing;
   int* particleSpeciesCodes;
   double* virial;

   LOG_INFORMATION("Getting data pointers");
   ier =
      KIM_ModelComputeArguments_GetArgumentPointerInteger(
          modelComputeArguments,
          KIM_COMPUTE_ARGUMENT_NAME_numberOfParticles,
          &nParts)
      ||
      KIM_ModelComputeArguments_GetArgumentPointerInteger(
          modelComputeArguments,
          KIM_COMPUTE_ARGUMENT_NAME_particleSpeciesCodes,
          &particleSpeciesCodes)
      ||
      KIM_ModelComputeArguments_GetArgumentPointerInteger(
          modelComputeArguments,
          KIM_COMPUTE_ARGUMENT_NAME_particleContributing,
          &particleContributing)
      ||
      KIM_ModelComputeArguments_GetArgumentPointerDouble(
          modelComputeArguments,
          KIM_COMPUTE_ARGUMENT_NAME_coordinates,
          &coords)
      ||
      KIM_ModelComputeArguments_GetArgumentPointerDouble(
          modelComputeArguments,
          KIM_COMPUTE_ARGUMENT_NAME_partialEnergy,
          &energy)
      ||
      KIM_ModelComputeArguments_GetArgumentPointerDouble(
          modelComputeArguments,
          KIM_COMPUTE_ARGUMENT_NAME_partialForces,
          &force)
      ||
      KIM_ModelComputeArguments_GetArgumentPointerDouble(
          modelComputeArguments,
          KIM_COMPUTE_ARGUMENT_NAME_partialParticleEnergy,
          &particleEnergy)
      ||
      KIM_ModelComputeArguments_GetArgumentPointerDouble(
      modelComputeArguments,
      KIM_COMPUTE_ARGUMENT_NAME_partialVirial,
      &virial);

      if (ier) {
         LOG_ERROR("get data pointers failed");
         return ier; }

   comp_energy = (energy != 0);
   comp_force = (force != 0);
   comp_particleEnergy = (particleEnergy != 0);
   comp_virial = (virial != 0);

   /* Check to be sure that the species are correct */
   /**/
   ier = TRUE; /* assume an error */
   for (i = 0; i < *nParts; ++i) {
      if ( SPECCODE != particleSpeciesCodes[i]) {
         LOG_ERROR("Unexpected species code detected");
         return ier; } }
   ier = FALSE;  /* everything is ok */

   /* initialize potential energies, forces, and virial term */
   /* Note: that the variable `particleEnergy' does not need to be initialized
    * because it's initial value is set during the embedding energy calculation.
    */
   if (comp_energy)
   {
      *energy = 0.0;
   }

   if (comp_force)
   {
      for (i = 0; i < *nParts; ++i)
      {
         for (k = 0; k < DIM; ++k)
         {
            force[i*DIM + k] = 0.0;
         }
      }
   }

   if (comp_virial)
   {
      for (i = 0; i < 6; ++i)
      {
         virial[i] = 0.0;
      }
   }

   /* pair functional electron density */
   rho = (double*) malloc((*nParts)*sizeof(double));
   for (i = 0; i < *nParts; ++i)
   {
      rho[i] = 0.0;
   }

   /* EAM embedded energy deriv */
   if ((comp_force) || (comp_virial))
   {
      derU = (double*) malloc((*nParts)*sizeof(double));
   }

   /* loop over particles in the neighbor list a first time,
    * to compute electron density (=coordination)
    */
   LOG_INFORMATION("Starting first compute loop for elctron density");
   for (i = 0; i< *nParts; ++i)
   {
      if (particleContributing[i])
      {
         ier = KIM_ModelComputeArguments_GetNeighborList(
            modelComputeArguments,
            0, i, &numOfPartNeigh, &neighListOfCurrentPart);
         if (ier) {
            /* some sort of problem, exit */
            LOG_ERROR("GetNeighborList failed");
            ier = TRUE;
            return ier; }

         /* loop over the neighbors of particle i */
         for (jj = 0; jj < numOfPartNeigh; ++ jj)
         {
            j = neighListOfCurrentPart[jj]; /* get neighbor ID */

            if (i < j) /* short-circuit half-list */
            {
              /* compute relative position vector and squared distance */
              Rsqij = 0.0;
              for (k = 0; k < DIM; ++k)
              {
                 Rij[k] = coords[j*DIM + k] - coords[i*DIM + k];

                 /* compute squared distance */
                 Rsqij += Rij[k]*Rij[k];
              }

              /* compute contribution to electron density */
              if (Rsqij < MODEL_CUTSQ) /* particles are interacting ? */
              {
                 r = sqrt(Rsqij); /* compute distance */
                 calc_g(r,&g);     /* compute electron density */
                 rho[i] += g;     /* accumulate electron density */

                 if (particleContributing[j])
                 {
                    rho[j] += g;
                 }

              }
            }
         } /* loop on jj */
      } /* if contributing */
   } /* loop on i */

   /* Now that we know the electron densities, calculate embedding part
    * of energy U and its derivative U' (derU)
    */
   LOG_INFORMATION("Starting second compute loop for embedding part of energy");
   for (i = 0; i < *nParts; ++i)
   {
      if(particleContributing[i])
      {
         if (comp_force || comp_virial)
         {
            calc_U_dU(rho[i], &U, &dU); /* compute embedding energy and its derivative */
            derU[i] = dU;               /* store dU for later use */
         }
         else
         {
            calc_U(rho[i],&U);         /* compute just embedding energy */
         }

         /* accumulate the embedding energy contribution */
         /* Assuming U)rho=0) = 0.0 */
         if (comp_particleEnergy)
         {
            particleEnergy[i] += U;
         }
         else if (comp_energy)
         {
            *energy += U;
         }
      }
   }

   /* Loop over particles in the neighbor list a second time to compute
    * the forces and complete energy calculation
    */

  LOG_INFORMATION("Starting thrid compute loop for forces and energy");
  for (i = 0; i< *nParts; ++i)
  {
      if (particleContributing[i])
      {
         ier = KIM_ModelComputeArguments_GetNeighborList(
             modelComputeArguments,
             0, i, &numOfPartNeigh, &neighListOfCurrentPart);
         if (ier) {
           /* some sort of problem, exit */
           LOG_ERROR("GetNeighborList failed");
           ier = TRUE;
           return ier; }
         /* loop over the neighbors of atom i */
         for (jj = 0; jj < numOfPartNeigh; ++ jj)
         {
            j = neighListOfCurrentPart[jj]; /* get neighbor ID */

            if (i < j) /* short-circuit half-list */
            {
              /* compute relative position vector and squared distance */
              Rsqij = 0.0;
              for (k = 0; k < DIM; ++k)
              {

                 Rij[k] = coords[j*DIM + k] - coords[i*DIM + k];

                 /* compute squared distance */
                 Rsqij += Rij[k]*Rij[k];
              }

              /* compute energy and force */
              if (Rsqij < MODEL_CUTSQ) /* particles are interacting ? */
              {
                 r = sqrt(Rsqij);
                 if (comp_force || comp_virial)
                 {
                    /* compute pair potential and its derivative */
                    calc_phi_dphi(r, &phi, &dphi);

                    /* copmute elect dens first deriv */
                    calc_dg(r, &dg);

                    /* compute dEidr */

                    /* HALF mode -- double contribution */
                    if(particleContributing[j])
                    {
                      dphii = 0.5*dphi;
                      dphij = 0.5*dphi;
                      dUi = derU[i]*dg;
                      dUj = derU[j]*dg;
                    }
                    else
                    {
                      dphii = 0.5*dphi;
                      dphij = 0.0;
                      dUi = derU[i]*dg;
                      dUj = 0.0;
                    }

                    dphieff = dphii + dphij + dUi + dUj;
                 }
                 else
                 {
                    /* compute just pair potential */
                    calc_phi(r, &phi);
                 }


                 /* HALF mode */
                 Ei = 0.5*phi;
                 if(particleContributing[j])
                    Ej = 0.5*phi;
                 else
                    Ej = 0.0;

                 /* contribution to energy */
                 if (comp_particleEnergy)
                 {
                    particleEnergy[i] += Ei; /* accumulate energy Ei */
                    particleEnergy[j] += Ej; /* accumulate energy Ej */
                 }
                 if (comp_energy)
                 {
                    *energy += Ei; /* accumulate energy */
                    *energy += Ej; /* accumulate energy */
                 }

                 /* contribution to virial tensor */
                 if (comp_virial)
                 {
                    /* virial(i,j) = r(i)*r(j)*(dV/dr)/r */
                    virial[0] += Rij[0]*Rij[0]*dphieff/r;
                    virial[1] += Rij[1]*Rij[1]*dphieff/r;
                    virial[2] += Rij[2]*Rij[2]*dphieff/r;
                    virial[3] += Rij[1]*Rij[2]*dphieff/r;
                    virial[4] += Rij[0]*Rij[2]*dphieff/r;
                    virial[5] += Rij[0]*Rij[1]*dphieff/r;
                 }

                 /* contribution to forces */
                 if (comp_force)
                 {  /* Ei contribution */
                    for (k = 0; k < DIM; ++k)
                    {
                       force[i*DIM + k] += dphieff*Rij[k]/r; /* accumulate force on atom i */
                       force[j*DIM + k] -= dphieff*Rij[k]/r; /* accumulate force on atom j */
                    }
                 }
              }
            }
         } /* loop on jj */
      } /* if contributing */
   } /* loop on i */
   LOG_INFORMATION("Finished compute loop");

   /* Free temporary storage */
   free(rho);
   if (comp_force || comp_virial)
   {
      free(derU);
   }

   /* everything is great */
   ier = FALSE;

   return ier;

}


/* Create function */
/* Define prototype for Model create */
#include "KIM_ModelCreateLogMacros.h"
int model_create(KIM_ModelCreate * const modelCreate,
                 KIM_LengthUnit const requestedLengthUnit,
                 KIM_EnergyUnit const requestedEnergyUnit,
                 KIM_ChargeUnit const requestedChargeUnit,
                 KIM_TemperatureUnit const requestedTemperatureUnit,
                 KIM_TimeUnit const requestedTimeUnit)
{
  buffer * bufferPointer;
  int error;

  /* set units */
  LOG_INFORMATION("Set model units");
  error = KIM_ModelCreate_SetUnits(
      modelCreate, /* ignoring requested units */
      KIM_LENGTH_UNIT_A,
      KIM_ENERGY_UNIT_eV,
      KIM_CHARGE_UNIT_unused,
      KIM_TEMPERATURE_UNIT_unused,
      KIM_TIME_UNIT_unused);

  /* register species */
  LOG_INFORMATION("Setting species code");
  error = error ||
      KIM_ModelCreate_SetSpeciesCode(modelCreate,
                                     KIM_SPECIES_NAME_Cu, SPECCODE);

  /* register numbering */
  LOG_INFORMATION("Setting model numbering");
  error = error || KIM_ModelCreate_SetModelNumbering(modelCreate,
                                                     KIM_NUMBERING_zeroBased);

  /* register function pointers */
  LOG_INFORMATION("Register model function pointers");
  error = error ||
      KIM_ModelCreate_SetComputePointer(
          modelCreate,
          KIM_LANGUAGE_NAME_c,
          (func *) &compute);
  error = error ||
      KIM_ModelCreate_SetComputeArgumentsCreatePointer(
          modelCreate,
          KIM_LANGUAGE_NAME_c,
          (func *) &compute_arguments_create);
  error = error ||
      KIM_ModelCreate_SetComputeArgumentsDestroyPointer(
          modelCreate,
          KIM_LANGUAGE_NAME_c,
          (func *) &compute_arguments_destroy);
  error = error ||
      KIM_ModelCreate_SetDestroyPointer(
          modelCreate,
          KIM_LANGUAGE_NAME_c,
          (func *) &model_destroy);
  error = error ||
      KIM_ModelCreate_SetRefreshPointer(
          modelCreate,
          KIM_LANGUAGE_NAME_c,
          (func *) &model_refresh);

  /* allocate buffer */
  bufferPointer = (buffer *) malloc(sizeof(buffer));

  /* store model buffer in KIM object */
  LOG_INFORMATION("Set influence distance and cutoffs");
  KIM_ModelCreate_SetModelBufferPointer(modelCreate,
                                        bufferPointer);

  /* set buffer values */
  bufferPointer->influenceDistance = MODEL_CUTOFF;
  bufferPointer->cutoff = MODEL_CUTOFF;
  bufferPointer->paddingNeighborHints = 1;
  bufferPointer->halfListHints = 1;

  /* register influence distance */
  KIM_ModelCreate_SetInfluenceDistancePointer(
      modelCreate,
      &(bufferPointer->influenceDistance));

  /* register cutoff */
  KIM_ModelCreate_SetNeighborListPointers(modelCreate, 1,
                                                &(bufferPointer->cutoff),
                                                (const int*) &(bufferPointer->paddingNeighborHints),
                                                (const int*) &(bufferPointer->halfListHints));

  if (error)
  {
    free(bufferPointer);
    LOG_ERROR("Unable to successfully initialize model");
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

/* refresh function */
#include "KIM_ModelRefreshLogMacros.h"
static int model_refresh(KIM_ModelRefresh * const modelRefresh)
{
  /* Local variables */
  buffer * bufferPointer;

  /* get model buffer from KIM object */
  LOG_INFORMATION("Getting model buffer");
  KIM_ModelRefresh_GetModelBufferPointer(modelRefresh,
                                         (void **) &bufferPointer);

  LOG_INFORMATION("Resetting influence distance and cutoffs");
  KIM_ModelRefresh_SetInfluenceDistancePointer(
      modelRefresh, &(bufferPointer->influenceDistance));
  KIM_ModelRefresh_SetNeighborListPointers(modelRefresh, 1,
                                                &(bufferPointer->cutoff),
                                                (const int*) &(bufferPointer->paddingNeighborHints),
                                                (const int*) &(bufferPointer->halfListHints));

  return FALSE;
}

/* Initialization function */
#include "KIM_ModelDestroyLogMacros.h"
int model_destroy(KIM_ModelDestroy * const modelDestroy) {
  buffer * bufferPointer;

  LOG_INFORMATION("Getting buffer");
  KIM_ModelDestroy_GetModelBufferPointer(modelDestroy,
                                         (void **) &bufferPointer);
  LOG_INFORMATION("Freeing model memory");
  free(bufferPointer);

  return FALSE; }

/* compute arguments create routine */
#include "KIM_ModelComputeArgumentsCreateLogMacros.h"
static int compute_arguments_create(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsCreate * const modelComputeArgumentsCreate)
{
  int error;
  /* register arguments */
  LOG_INFORMATION("Register argument supportStatus");
  error =
      KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
          modelComputeArgumentsCreate,
          KIM_COMPUTE_ARGUMENT_NAME_partialEnergy, KIM_SUPPORT_STATUS_optional);

  error = error ||
      KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
          modelComputeArgumentsCreate,
          KIM_COMPUTE_ARGUMENT_NAME_partialForces, KIM_SUPPORT_STATUS_optional);

  error = error ||
      KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
          modelComputeArgumentsCreate,
          KIM_COMPUTE_ARGUMENT_NAME_partialParticleEnergy,
          KIM_SUPPORT_STATUS_optional);

  error = error ||
    KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
        modelComputeArgumentsCreate,
        KIM_COMPUTE_ARGUMENT_NAME_partialVirial,
        KIM_SUPPORT_STATUS_optional);

  /* register call backs */
  LOG_INFORMATION("Register call back supportStatus");

  if (error)
  {
    LOG_ERROR("Unable to successfully initialize compute arguments");
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

/* compue arguments destroy routine */
#include "KIM_ModelComputeArgumentsDestroyLogMacros.h"
static int compute_arguments_destroy(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy)
{
  /* nothing to be done */

  return FALSE;

}
