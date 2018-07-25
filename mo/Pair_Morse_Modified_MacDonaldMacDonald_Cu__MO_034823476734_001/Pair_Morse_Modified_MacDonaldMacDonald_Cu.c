/*                                                                            */
/* CDDL HEADER START                                                          */
/*                                                                            */
/* The contents of this file are subject to the terms of the Common           */
/* Development and Distribution License Version 1.0 (the "License").          */
/*                                                                            */
/* You can obtain a copy of the license at                                    */
/* http://www.opensource.org/licenses/CDDL-1.0.  See the License for the      */
/* specific language governing permissions and limitations under the License. */
/*                                                                            */
/* When distributing Covered Code, include this CDDL HEADER in each file and  */
/* include the License file in a prominent location with the name             */
/* LICENSE.CDDL. If applicable, add the following below this CDDL HEADER,     */
/* with the fields enclosed by brackets "[]" replaced with your own           */
/* identifying information:                                                   */
/*                                                                            */
/* Portions Copyright (c) [yyyy] [name of copyright owner].                   */
/* All rights reserved.                                                       */
/*                                                                            */
/* CDDL HEADER END                                                            */
/*                                                                            */

/*                                                                            */
/* Copyright (c) 2013--2014, Regents of the University of Minnesota.          */
/* All rights reserved.                                                       */
/*                                                                            */
/* Contributors:                                                              */
/*    Ryan S. Elliott                                                         */
/*    Ellad B. Tadmor                                                         */
/*    Valeriu Smirichinski                                                    */
/*    Stephen M. Whalen                                                       */
/*                                                                            */
/*                                                                            */
/* Portions Copyright (c) 2014, Regents of the University of Minnesota.       */
/* All rights reserved.                                                       */
/*                                                                            */
/* Contributors:                                                              */
/*    Mingjian Wen                                                            */
/*                                                                            */

/******************************************************************************/
/*                                                                            */
/*  Pair_Morse_Modified_MacDonaldMacDonald_Cu                                 */
/*                                                                            */
/*  Modified Morse pair potential model for Cu                                */
/*                                                                            */
/*  Reference: MacDonald, Rosemary A., and MacDonald, William M.              */
/*             "Thermodynamic properties of fcc metals at high temperatures." */
/*             Physical Review B 24.4 (1981): 1715.                           */
/*                                                                            */
/*  Language: C                                                               */
/*                                                                            */
/*  Release: This file is part of the the paper:                              */
/*           Interpolation effects in tabulated potentials,                   */
/*           submitted to Modelling Simulation Mater. Sci. Eng.               */
/*                                                                            */
/******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "KIM_ModelHeaders.h"

#define TRUE 1
#define FALSE 0

/******************************************************************************/
/* Below are the definitions and values of all Model parameters               */
/******************************************************************************/
#define DIM 3       /* dimensionality of space */
#define SPECCODE 1  /* internal species code */
#define MODEL_CUTOFF 1000  /* cutoff radius in angstroms */
#define MODEL_CUTSQ  (MODEL_CUTOFF * MODEL_CUTOFF)
#define r0 2.5471
#define D0 5.868891232E-1
#define A        1.1857
#define B  2.265

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


/* Calculate pair potential phi(r) */
static void calc_phi(double r, double* phi)
{

  if (r > MODEL_CUTOFF)
  {
    /* Argument exceeds cutoff radius */
    *phi = 0.0 ;
  }
  else
  {
    *phi = (  exp( -2*A*sqrt(B)*(r-r0) )
              - 2.0 *B*exp( -A*(r-r0)/sqrt(B) )  )*D0/(2.0*B-1.0) ;
  }

  return;
}

/* Calculate pair potential phi(r) and its derivative dphi(r) */
static void calc_phi_dphi(double r, double* phi, double* dphi)
{

  if (r > MODEL_CUTOFF)
  {
    /* Argument exceeds cutoff radius */
    *phi  = 0.0;
    *dphi = 0.0;
  }
  else
  {
    *phi = (  exp( -2.0*A*sqrt(B)*(r-r0) )
              - 2.0 *B*exp( -A*(r-r0)/sqrt(B) )  )*D0/(2.0*B-1.0);
    *dphi = (  -2.0*A*sqrt(B)* exp( -2.0*A*sqrt(B)*(r-r0) )
               +2.0*A*sqrt(B)*exp( -A*(r-r0)/sqrt(B) )  )*D0/(2*B-1);
  }

  return;
}

/* compute function */
#include "KIM_ModelComputeLogMacros.h"
static int compute(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArguments const * const modelComputeArguments)
{
  double R;
  double Rsqij;
  double phi;
  double dphi;
  double dEidr = 0.0;
  double Rij[DIM];
  int ier;
  int i;
  int j;
  int jj;
  int k;
  int const * neighListOfCurrentPart;
  int numOfPartNeigh;
  int comp_energy;
  int comp_force;
  int comp_particleEnergy;
  int comp_virial;
  // int comp_particleVirial;

  int* nParts;
  int* particleContributing;
  int* particleSpeciesCodes;
  double* Rij_list;
  double* coords;
  double* energy;
  double* force;
  double* particleEnergy;
  double* virial;
  // double* particleVirial;

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
      // ||
      // KIM_ModelComputeArguments_GetArgumentPointerDouble(
      // modelComputeArguments,
      // KIM_COMPUTE_ARGUMENT_NAME_partialParticleVirial,
      // &particleVirial);

  if (ier) {
    LOG_ERROR("get data pointers failed");
    return ier; }

  comp_energy = (energy != 0);
  comp_force = (force != 0);
  comp_particleEnergy = (particleEnergy != 0);
  comp_virial = (virial != 0);
  // comp_particleVirial = (particleVirial != 0);

  /* Check to be sure that the species are correct */
  /**/
  ier = TRUE; /* assume an error */
  for (i = 0; i < *nParts; ++i) {
    if ( SPECCODE != particleSpeciesCodes[i]) {
      LOG_ERROR("Unexpected species code detected");
      return ier; } }
  ier = FALSE;  /* everything is ok */

  /* initialize potential energies, forces, and virial term */
  if (comp_particleEnergy)
  {
    for (i = 0; i < *nParts; ++i)
    {
      particleEnergy[i] = 0.0;
    }
  }
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

  /* loop over particles and compute enregy and forces */
  LOG_INFORMATION("Starting main compute loop");
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

          /* compute energy and force */
          if (Rsqij < MODEL_CUTSQ) /* particles are interacting ? */
          {
            R = sqrt(Rsqij);
            if (comp_force || comp_virial)
            {
              /* compute pair potential and its derivative */
              calc_phi_dphi(R, &phi, &dphi);

              /* compute dEidr */
              if (particleContributing[j])
              {
                dEidr = dphi;
              }
              else
              {
                dEidr = 0.5*dphi;
              }

            }
            else
            {
              /* compute just pair potential */
              calc_phi(R, &phi);
            }

            /* contribution to energy */
            if (comp_particleEnergy)
            {
              particleEnergy[i] += 0.5*phi;
              if (particleContributing[j])
              {
                particleEnergy[j] += 0.5*phi;
              }
            }
            if (comp_energy)
            {
              if (particleContributing[j])
              {
                *energy += phi;
              }
              else
              {
                *energy += 0.5*phi;
              }
            }
            /* contribution to virial tensor */
            if (comp_virial)
            {
              /* virial(i,j) = r(i)*r(j)*(dV/dr)/r */
              virial[0] += Rij[0]*Rij[0]*dEidr/R;
              virial[1] += Rij[1]*Rij[1]*dEidr/R;
              virial[2] += Rij[2]*Rij[2]*dEidr/R;
              virial[3] += Rij[1]*Rij[2]*dEidr/R;
              virial[4] += Rij[0]*Rij[2]*dEidr/R;
              virial[5] += Rij[0]*Rij[1]*dEidr/R;
            }

            /* contribution to forces */
            if (comp_force)
            {
              for (k = 0; k < DIM; ++k)
              {
                force[i*DIM + k] += dEidr*Rij[k]/R; /* accumulate force on i */
                force[j*DIM + k] -= dEidr*Rij[k]/R; /* accumulate force on j */
              }
            }
          }
        }
      } /* loop on jj */
    } /* if contributing */
  }    /*  loop on i */
  LOG_INFORMATION("Finished compute loop");


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

  // error = error ||
  //     KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
  //         modelComputeArgumentsCreate,
  //         KIM_COMPUTE_ARGUMENT_NAME_partialParticleVirial,
  //         KIM_SUPPORT_STATUS_optional);

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
