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
 * Copyright (c) 2013--2014, Regents of the University of Minnesota.
 * All rights reserved.
 *
 * Contributors:
 *    Ryan S. Elliott
 *    Ellad B. Tadmor
 *    Valeriu Smirichinski
 *    Stephen M. Whalen
 *
 */

/*******************************************************************************
 *
 *  Pair_Morse_Smoothed
 *
 *  Morse pair potential KIM Model Driver
 *  smoothed to have zero energy, force, and stiffness at the cutoff radius
 *
 *  cutoff function is:
 *     g(r) = (r-cutoff)^3 * (g[0] + g[1]*(r-rstar) + 0.5*g[2]*(r-rstar)^2
 *  this gives C2 continuity and zero values (phi, dphi, d2phi) at cutoff.
 *
 *  Language: C
 *
 ******************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "KIM_ModelDriverHeaders.h"

#define TRUE 1
#define FALSE 0
#define ONE 1.0


/******************************************************************************
* Below are the definitions and values of all Model parameters
*******************************************************************************/
#define DIM 3       /* dimensionality of space */
#define SPECCODE 1  /* internal species code */
#define GORDER 3    /* number of coeffs for g(r) */

/* Define prototypes for Model Driver init */
/**/
/* Define prototype for Model Driver initialization */
int model_driver_create(KIM_ModelDriverCreate *const modelDriverCreate,
                        KIM_LengthUnit const requestedLengthUnit,
                        KIM_EnergyUnit const requestedEnergyUnit,
                        KIM_ChargeUnit const requestedChargeUnit,
                        KIM_TemperatureUnit const requestedTemperatureUnit,
                        KIM_TimeUnit const requestedTimeUnit);

/* Define prototype for destroy (defined as static to avoid namespace clashes with other codes) */
static int destroy(KIM_ModelDestroy *const modelDestroy);

/* Define prototype for compute routine */
static int compute(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArguments const * const modelComputeArguments);
static int compute_arguments_create(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsCreate * const modelComputeArgumentsCreate);
static int compute_arguments_destroy(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy);

/* Define prototype for refresh routine */
static int refresh(KIM_ModelRefresh * const modelRefresh);
/**/
static void calc_phi(double* const epsilon,
                     double* const C,
                     double* const Rzero,
                     double const g[GORDER],
                     double* const cutoff,
                     double* const rstar,
                     double const r,
                     double* const phi);
static void calc_phi_dphi(double* const epsilon,
                          double* const C,
                          double* const Rzero,
                          double const g[GORDER],
                          double* const cutoff,
                          double* const rstar,
                          double const r,
                          double* const phi,
                          double* const dphi);
static void calc_phi_d2phi(double* const epsilon,
                           double* const C,
                           double* const Rzero,
                           double const g[GORDER],
                           double* const cutoff,
                           double* const rstar,
                           double const r,
                           double* const phi,
                           double* const dphi,
                           double* const d2phi);
static void compute_gr(double* const epsilon,
                       double* const C,
                       double* const Rzero,
                       double* const cutoff,
                       double* const rstar,
                       double g[GORDER]);

/* Define model_buffer structure */
struct model_buffer {
  int paddingNeighborHints;
  int halfListHints;

  int energy_ind;
  int forces_ind;
  int particleEnergy_ind;
  int process_dEdr_ind;
  int process_d2Edr2_ind;
  int numberOfParticles_ind;
  int particleSpecies_ind;
  int coordinates_ind;
  int numberContributingParticles_ind;
  int get_neigh_ind;
  int cutoff_ind;

  double influenceDistance;
  double transition;
  double Pcutoff;
  double cutsq;
  double rstar;
  double epsilon;
  double C;
  double Rzero;
  double shift;
  double g[GORDER];
};


/* Calculate pair potential phi(r) */
static void calc_phi(double* const epsilon,
                     double* const C,
                     double* const Rzero,
                     double const g[GORDER],
                     double* const cutoff,
                     double* const rstar,
                     double const r,
                     double* const phi)
{
  /* local variables */
  double const dstar = r - *rstar;
  double const dcut = r - *cutoff;
  double const ep = exp(-(*C)*(r-*Rzero));
  double const ep2 = ep*ep;

  if (r < *rstar)
  {
    *phi  = (*epsilon)*( -ep2 + 2.0*ep );
  }
  else if (r < *cutoff)
  {
    *phi  = (dcut*dcut*dcut)*(g[0] + g[1]*dstar + 0.5*g[2]*dstar*dstar);
  }
  else
  {
    /* Argument exceeds cutoff radius */
    *phi  = 0.0;
  }

  return;
}

/* Calculate pair potential phi(r) and its derivative dphi(r) */
static void calc_phi_dphi(double* const epsilon,
                     double* const C,
                     double* const Rzero,
                     double const g[GORDER],
                     double* const cutoff,
                     double* const rstar,
                     double const r,
                     double* const phi,
                     double* const dphi)
{
  /* local variables */
  double const dstar = r - *rstar;
  double const dcut = r - *cutoff;
  double const ep = exp(-(*C)*(r-*Rzero));
  double const ep2 = ep*ep;

  if (r < *rstar)
  {
    *phi  = (*epsilon)*( -ep2 + 2.0*ep );
    *dphi = 2.0*(*epsilon)*(*C)*( -ep + ep2 );
  }
  else if (r < *cutoff)
  {
    *phi  = (dcut*dcut*dcut)*(g[0] + g[1]*dstar + 0.5*g[2]*dstar*dstar);
    *dphi = 3.0*(dcut*dcut)*(g[0] + g[1]*dstar + 0.5*g[2]*dstar*dstar)
        + (dcut*dcut*dcut)*(g[1] + g[2]*dstar);
  }
  else
  {
    /* Argument exceeds cutoff radius */
    *phi  = 0.0;
    *dphi = 0.0;
  }

  return;
}

/*
  Calculate pair potential phi(r) and its 1st & 2nd derivatives dphi(r),
  d2phi(r)
*/
static void calc_phi_d2phi(double* const epsilon,
                           double* const C,
                           double* const Rzero,
                           double const g[GORDER],
                           double* const cutoff,
                           double* const rstar,
                           double const r,
                           double* const phi,
                           double* const dphi,
                           double* const d2phi)
{
  /* local variables */
  double const dstar = r - *rstar;
  double const dcut = r - *cutoff;
  double const ep = exp(-(*C)*(r-*Rzero));
  double const ep2 = ep*ep;

  if (r < *rstar)
  {
    *phi   = (*epsilon)*( -ep2 + 2.0*ep );
    *dphi  = 2.0*(*epsilon)*(*C)*( -ep + ep2 );
    *d2phi = 2.0*(*epsilon)*(*C)*(*C)*(ep - 2.0*ep2);
  }
  else if (r < *cutoff)
  {
    *phi   = (dcut*dcut*dcut)*(g[0] + g[1]*dstar + 0.5*g[2]*dstar*dstar);
    *dphi  = 3.0*(dcut*dcut)*(g[0] + g[1]*dstar + 0.5*g[2]*dstar*dstar)
        + (dcut*dcut*dcut)*(g[1] + g[2]*dstar);
    *d2phi = 6.0*dcut*(g[0] + g[1]*dstar + 0.5*g[2]*dstar*dstar)
        + 6.0*dcut*dcut*(g[1] + g[2]*dstar)
        + dcut*dcut*dcut*g[2];
  }
  else
  {
    /* Argument exceeds cutoff radius */
    *phi   = 0.0;
    *dphi  = 0.0;
    *d2phi = 0.0;
  }

  return;
}

/* ######## Function Definitions ######## */

/* Primary compute function*/
#include "KIM_ModelComputeLogMacros.h"
static int compute(KIM_ModelCompute const * const modelCompute,
                   KIM_ModelComputeArguments const * const modelComputeArguments)
{
  /* local variables */
  double R;
  double R_pairs[2];
  double *pR_pairs = &(R_pairs[0]);
  double Rsqij;
  double phi;
  double dphi;
  double d2phi;
  double dEidr;
  double d2Eidr;
  double Rij[DIM];
  double *pRij = &(Rij[0]);
  double Rij_pairs[2][3];
  double *pRij_pairs = &(Rij_pairs[0][0]);
  int ier;
  int i;
  int i_pairs[2];
  int *pi_pairs = &(i_pairs[0]);
  int j;
  int j_pairs[2];
  int *pj_pairs = &(j_pairs[0]);
  int jj;
  int k;
  int const * neighListOfCurrentAtom;
  struct model_buffer* buffer;
  int comp_energy;
  int comp_force;
  int comp_particleEnergy;
  int comp_process_dEdr;
  int comp_process_d2Edr2;

  int* nAtoms;
  int* particleSpeciesCodes;
  int *particleContributing;
  double* cutoff;
  double* cutsq;
  double* epsilon;
  double* C;
  double* Rzero;
  double* rstar;
  double* coords;
  double* energy;
  double* force;
  double* particleEnergy;
  int numOfAtomNeigh;

  /* Get model buffer */
  KIM_ModelCompute_GetModelBufferPointer(modelCompute, (void **) &buffer);

  /* unpack info from the buffer */
  cutoff = &(buffer->Pcutoff);
  cutsq = &(buffer->cutsq);
  epsilon = &(buffer->epsilon);
  C = &(buffer->C);
  Rzero = &(buffer->Rzero);
  rstar = &(buffer->rstar);

  /* check to see if we have been asked to compute the forces, */
  /* particleEnergy, and d1Edr */
  LOG_INFORMATION("Checking if call backs are present.");
  KIM_ModelComputeArguments_IsCallbackPresent(
    modelComputeArguments,
    KIM_COMPUTE_CALLBACK_NAME_ProcessDEDrTerm,
    &comp_process_dEdr);
  KIM_ModelComputeArguments_IsCallbackPresent(
    modelComputeArguments,
    KIM_COMPUTE_CALLBACK_NAME_ProcessD2EDr2Term,
    &comp_process_d2Edr2);

  LOG_INFORMATION("Getting data pointers");
  ier =
  KIM_ModelComputeArguments_GetArgumentPointerInteger(
     modelComputeArguments,
     KIM_COMPUTE_ARGUMENT_NAME_numberOfParticles,
     &nAtoms)
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
    &particleEnergy);

  if (ier == TRUE){
    LOG_ERROR("Could not retrieve compute arguments");
    return ier;
  }

  /* Check to see what we've been asked to compute */
  comp_energy = (energy != NULL);
  comp_force = (force != NULL);
  comp_particleEnergy = (particleEnergy != NULL);

  /* Check to be sure that the species are correct */
  ier = TRUE;
  for(i = 0; i < *nAtoms; ++i)
  {
    if (SPECCODE != particleSpeciesCodes[i]){
      LOG_ERROR("Unexpected species code detected");
      return ier;
    }
  }
  ier = FALSE;

  if (comp_particleEnergy){
    for (i = 0; i < *nAtoms; ++i)
    {
      particleEnergy[i] = 0.0;
    }
  }
  if (comp_energy){
    *energy = 0.0;
  }

  if (comp_force){
    for (i = 0; i < *nAtoms; ++i)
    {
      for (k = 0; k < DIM; ++k)
      {
        force[i*DIM + k] = 0.0;
      }
    }
  }


  /* Compute energy and forces */

  /* loop over particles and compute enregy and forces */
  LOG_INFORMATION("Starting main compute loop");
  for (i = 0; i< *nAtoms; ++i)
  {
    if (particleContributing[i])
    {
      ier = KIM_ModelComputeArguments_GetNeighborList(
          modelComputeArguments,
          0, i, &numOfAtomNeigh, &neighListOfCurrentAtom);
      if (ier){
        /* some sort of problem, exit */
        LOG_ERROR("GetNeighborList failed");
        ier = TRUE;
        return ier;
      }



      /* loop over the neighbors of atom i */
      for (jj = 0; jj < numOfAtomNeigh; ++ jj)
      {
        /* get neighbor ID */
        j = neighListOfCurrentAtom[jj];

        if (i < j) /* short-circuit half-list */
        {
          /* compute relative position vector and squared distance */
          Rsqij = 0.0;
          for (k = 0; k < DIM; ++k){
            Rij[k] = coords[j*DIM + k] - coords[i*DIM + k];

            /* compute squared distance */
            Rsqij += Rij[k]*Rij[k];
          }

          /* compute energy and force */
          if (Rsqij < *cutsq) /* particles are interacting ? */
          {
            R = sqrt(Rsqij);
            if (comp_process_d2Edr2){
              /* compute pair potential and its derivatives */
              calc_phi_d2phi(epsilon,
                              C,
                              Rzero,
                              buffer->g,
                              cutoff,
                              rstar,
                              R, &phi, &dphi, &d2phi);

              /* compute dEidr */
              /* Half mode -- double contribution */
              if (particleContributing[j])
              {
                dEidr = dphi;
                d2Eidr = d2phi;
              }
              else
              {
                dEidr = 0.5*dphi;
                d2Eidr = 0.5*d2phi;
              }
            }
            else if (comp_force || comp_process_dEdr){
               /* compute pair potential and its derivative */
              calc_phi_dphi(epsilon,
                            C,
                            Rzero,
                            buffer->g,
                            cutoff,
                            rstar,
                            R, &phi, &dphi);

              /* compute dEidr */
              /* Half mode -- double contribution */
              if (particleContributing[j])
              {
                dEidr = dphi;
              }
              else
              {
                dEidr = 0.5*dphi;
              }

            }
            else{
              /* compute just pair potential */
              calc_phi(epsilon,
                        C,
                        Rzero,
                        buffer->g,
                        cutoff,
                        rstar,
                        R, &phi);
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

            /* contribution to process_dEdr */
            if (comp_process_dEdr){
              ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                        modelComputeArguments, dEidr, R, pRij, i, j);
            }

            /* contribution to process_d2Edr2 */
            if (comp_process_d2Edr2){
              R_pairs[0] = R_pairs[1] = R;
              Rij_pairs[0][0] = Rij_pairs[1][0] = Rij[0];
              Rij_pairs[0][1] = Rij_pairs[1][1] = Rij[1];
              Rij_pairs[0][2] = Rij_pairs[1][2] = Rij[2];
              i_pairs[0] = i_pairs[1] = i;
              j_pairs[0] = j_pairs[1] = j;

              ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                    modelComputeArguments, d2Eidr, pR_pairs, pRij_pairs,
                                                     pi_pairs, pj_pairs);
            }

            /* contribution to forces */
            if (comp_force){
              for (k = 0; k < DIM; ++k)
              {  /* accumulate force on atom i */
                force[i*DIM + k] += dEidr*Rij[k]/R;
                /* accumulate force on atom j */
                force[j*DIM + k] -= dEidr*Rij[k]/R;
              }
            }
          }
        }
      } /* loop on jj */
    } /* if contributing */
  }    /* loop on i */
  LOG_INFORMATION("Finished compute loop");

  /* everything is great */
  ier = FALSE;

  return ier;

}

/* Define prototype for Model Driver initialization */
#include "KIM_ModelDriverCreateLogMacros.h"
int model_driver_create(KIM_ModelDriverCreate *const modelDriverCreate,
                        KIM_LengthUnit const requestedLengthUnit,
                        KIM_EnergyUnit const requestedEnergyUnit,
                        KIM_ChargeUnit const requestedChargeUnit,
                        KIM_TemperatureUnit const requestedTemperatureUnit,
                        KIM_TimeUnit const requestedTimeUnit)
{
  // Create Model buffer
  struct model_buffer *buffer;

  // KIM variables
  int numberOfParameterFiles;
  char const *paramfile1name; /* This driver only actually supports a single parameter file */
  char speciesNameString[100];
  KIM_SpeciesName speciesName;

  // Model Parameters
  double a;
  double transition;
  double cutoff;
  double epsilon;
  double C;
  double Rzero;
  int ier;

  KIM_LengthUnit fromLength;
  KIM_EnergyUnit fromEnergy;
  KIM_ChargeUnit fromCharge;
  KIM_TemperatureUnit fromTemperature;
  KIM_TimeUnit fromTime;

  // Misc
  FILE *fid;
  double convertLength = 1.0;
  double convertEnergy = 1.0;

  ier = FALSE;

  /* Get number of parameter files and error out if more than one is listed */
  KIM_ModelDriverCreate_GetNumberOfParameterFiles(modelDriverCreate, &numberOfParameterFiles);
  if (numberOfParameterFiles != 1){
    ier = TRUE;
    LOG_ERROR("Incorrect number of parameter files. Only one parameter file is currently \
             supported by this driver.");
    return ier;
  }

  /* Set paramfile1name */
  ier = KIM_ModelDriverCreate_GetParameterFileName(modelDriverCreate,
                                                 0,
                                                 &paramfile1name);
  if (ier == TRUE){
    LOG_ERROR("Unable to get parameter file name");
    return ier;
  }

  /* Read in Model parameters from parameter file */
  fid = fopen(paramfile1name, "r");
  if(fid == NULL){
    ier = TRUE;
    LOG_ERROR("Unable to open parameter file for Morse parameters");
    return ier;
  }

  ier = fscanf(fid, "%s \n%lf \n%lf \n%lf \n%lf \n%lf",
              speciesNameString,
              &a,  /* cutoff distance in angstroms */
              &transition, /* starting point of tail in fraction of cutoff */
              &epsilon, /* Morse epsilon in eV */
              &C,       /* Morse C in 1/Angstroms */
              &Rzero);  /* Morse Rzero in Angstroms */
  fclose(fid);

  /* check that we read the right number of parameters */
  if (6 != ier){
    ier = TRUE;
    LOG_ERROR("Unable to read all Morse parameters");
    return ier;
  }
  /* Define default base units */
  fromLength = KIM_LENGTH_UNIT_A;
  fromEnergy = KIM_ENERGY_UNIT_eV;
  fromCharge = KIM_CHARGE_UNIT_e;
  fromTemperature = KIM_TEMPERATURE_UNIT_K;
  fromTime = KIM_TIME_UNIT_ps;

  /* Convert units of Model parameters */
  ier = KIM_ModelDriverCreate_ConvertUnit(
        modelDriverCreate,
        fromLength, fromEnergy, fromCharge, fromTemperature, fromTime,
        requestedLengthUnit, requestedEnergyUnit, requestedChargeUnit,
        requestedTemperatureUnit, requestedTimeUnit,
        1.0, 0.0, 0.0, 0.0, 0.0,
        &convertLength);
  if(ier){
    LOG_ERROR("Unable to convert length unit");
    return ier;
  }

  ier = KIM_ModelDriverCreate_ConvertUnit(
        modelDriverCreate,
        fromLength, fromEnergy, fromCharge, fromTemperature, fromTime,
        requestedLengthUnit, requestedEnergyUnit, requestedChargeUnit,
        requestedTemperatureUnit, requestedTimeUnit,
        0.0, 1.0, 0.0, 0.0, 0.0,
        &convertEnergy);
  if(ier){
    LOG_ERROR("Unable to convert energy unit");
    return ier;
  }

  if(convertLength != ONE){
    Rzero *= convertLength;
    a *= convertLength;
    C *= (1.0/convertLength);
  }

  if(convertEnergy != ONE){
    epsilon *= convertEnergy;
  }

   /* Set cutoff equal to 'a' and work with it from now on */
  cutoff = a;

  /* Set units */
  LOG_INFORMATION("Setting units");
  ier = KIM_ModelDriverCreate_SetUnits(
          modelDriverCreate,
          requestedLengthUnit,
          requestedEnergyUnit,
          KIM_CHARGE_UNIT_unused,
          KIM_TEMPERATURE_UNIT_unused,
          KIM_TIME_UNIT_unused);
  if (ier == TRUE){
    LOG_ERROR("Unable to set units to requested set");
  }

  /* Indicate that we use zero-based indexing for particles in this driver */
  LOG_INFORMATION("Setting particle indexing to zero-based");
  ier = KIM_ModelDriverCreate_SetModelNumbering(modelDriverCreate, KIM_NUMBERING_zeroBased);
  if (ier == TRUE){
    LOG_ERROR("Unable to set particle indexing");
    return ier;
  }

  /* Store pointers to functions in the KIM object */
  LOG_INFORMATION("Registering Model function pointers");
  ier = ier || KIM_ModelDriverCreate_SetDestroyPointer(modelDriverCreate,
                                        KIM_LANGUAGE_NAME_c, (func *) &model_driver_create);
  ier = ier || KIM_ModelDriverCreate_SetComputeArgumentsCreatePointer(modelDriverCreate,
                                        KIM_LANGUAGE_NAME_c, (func *) &compute_arguments_create);
  ier = ier || KIM_ModelDriverCreate_SetComputeArgumentsDestroyPointer(modelDriverCreate,
                                        KIM_LANGUAGE_NAME_c, (func *) &compute_arguments_destroy);
  ier = ier || KIM_ModelDriverCreate_SetDestroyPointer(modelDriverCreate,
                                        KIM_LANGUAGE_NAME_c, (func *) &destroy);
  ier = ier || KIM_ModelDriverCreate_SetComputePointer(modelDriverCreate,
                                        KIM_LANGUAGE_NAME_c, (func *) &compute);
  ier = ier || KIM_ModelDriverCreate_SetRefreshPointer(modelDriverCreate,
                                        KIM_LANGUAGE_NAME_c, (func *) &refresh);
  if (ier == TRUE){
    LOG_ERROR("Unable to register Model function pointers");
    return ier;
  }

  /* Register species */
  LOG_INFORMATION("Setting species code");
  speciesName = KIM_SpeciesName_FromString(speciesNameString);
  ier = ier || KIM_ModelDriverCreate_SetSpeciesCode(modelDriverCreate,
                                                  speciesName, SPECCODE);
  if(ier == TRUE){
     LOG_ERROR("Unable to set species code");
     return ier;
  }

  /* Allocate buffer */
  LOG_INFORMATION("Allocating memory for Model buffer");
  buffer = (struct model_buffer*) malloc(sizeof(struct model_buffer));
  if (buffer == NULL){
    ier = TRUE;
    LOG_ERROR("Could not allocate Model buffer");
    return ier;
  }

  /* Store model buffer in KIM object */
  LOG_INFORMATION("Registering Model buffer");
  KIM_ModelDriverCreate_SetModelBufferPointer(modelDriverCreate, (void*) buffer);

  /* Store parameter values in buffer */
  LOG_INFORMATION("Loading Model parameters into buffer");
  buffer->rstar = transition*cutoff;
  buffer->transition = transition;
  buffer->influenceDistance = cutoff;
  buffer->cutsq = cutoff*cutoff;
  buffer->epsilon = epsilon;
  buffer->Rzero = Rzero;
  buffer->C = C;
  buffer->Pcutoff = cutoff;
  buffer->paddingNeighborHints = 1;
  buffer->halfListHints = 1;

  /* calc buf->g[] */
  compute_gr(&(buffer->epsilon),
             &(buffer->C),
             &(buffer->Rzero),
             &(cutoff),
             &(buffer->rstar),
             buffer->g);

  /* Register influence distance pointer */
  LOG_INFORMATION("Registering influence distance pointer");
  KIM_ModelDriverCreate_SetInfluenceDistancePointer(modelDriverCreate, &(buffer->influenceDistance));

  /* Register cutoff pointer */
  LOG_INFORMATION("Registering cutoff pointer");
  KIM_ModelDriverCreate_SetNeighborListPointers(modelDriverCreate, 1,
                                                &(buffer->Pcutoff),
                                                (const int*) &(buffer->paddingNeighborHints),
                                                (const int*) &(buffer->halfListHints));
  /* Register Model parameters */
  LOG_INFORMATION("Registering Model parameters");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->Pcutoff), "cutoff");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->Rzero), "Rzero");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->C), "C");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->epsilon), "epsilon");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->transition), "transition");


  if (ier == TRUE){
    free(buffer);
    LOG_ERROR("Unable to successfully initialize Model");
  }

  return ier;
}

/* Refresh function */
#include "KIM_ModelRefreshLogMacros.h"
static int refresh(KIM_ModelRefresh * const modelRefresh){
  /* Local variables */
  double* cutoff;
  struct model_buffer* buffer;

  /* Get model buffer from KIM object */
  LOG_INFORMATION("Retrieving Model buffer");
  KIM_ModelRefresh_GetModelBufferPointer(modelRefresh, (void **) &buffer);

  LOG_INFORMATION("Resetting influence distance and cutoff, and shift");
  *cutoff = buffer->Pcutoff;

  /* set value of parameter cutsq */
  buffer->cutsq = (*cutoff)*(*cutoff);
  buffer->rstar = (buffer->transition)*(*cutoff);

  /* calc buf->g[] */
  compute_gr(&(buffer->epsilon),
             &(buffer->C),
             &(buffer->Rzero),
             cutoff,
             &(buffer->rstar),
             buffer->g);

  KIM_ModelRefresh_SetInfluenceDistancePointer(
      modelRefresh, &(buffer->influenceDistance));
  KIM_ModelRefresh_SetNeighborListPointers(modelRefresh, 1,
                                              &(buffer->Pcutoff),
                                              (const int*) &(buffer->paddingNeighborHints),
                                              (const int*) &(buffer->halfListHints));
  return FALSE;
}

/* Destroy function */
#include "KIM_ModelDestroyLogMacros.h"
static int destroy(KIM_ModelDestroy * const modelDestroy){
  /* Local variables */
  struct model_buffer* buffer;

  /* Get model buffer from KIM object */
  KIM_ModelDestroy_GetModelBufferPointer(modelDestroy, (void **) &buffer);

  /* Free the buffer */
  LOG_INFORMATION("Deallocating Model buffer");
  free(buffer);

  return FALSE;
}

/* compute_arguments create routine */
#include "KIM_ModelComputeArgumentsCreateLogMacros.h"
static int compute_arguments_create(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsCreate * const modelComputeArgumentsCreate){
  int ier;

  ier = FALSE;

  /* Register arguments */
  ier = ier || KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
                                       modelComputeArgumentsCreate,
                                       KIM_COMPUTE_ARGUMENT_NAME_partialEnergy,
                                       KIM_SUPPORT_STATUS_optional);
  ier = ier || KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
                                       modelComputeArgumentsCreate,
                                       KIM_COMPUTE_ARGUMENT_NAME_partialParticleEnergy,
                                       KIM_SUPPORT_STATUS_optional);
  ier = ier || KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
                                       modelComputeArgumentsCreate,
                                       KIM_COMPUTE_ARGUMENT_NAME_partialForces,
                                       KIM_SUPPORT_STATUS_optional);

    /* register call backs */
  LOG_INFORMATION("Register call back supportStatus");
  ier = ier ||
      KIM_ModelComputeArgumentsCreate_SetCallbackSupportStatus(
          modelComputeArgumentsCreate,
          KIM_COMPUTE_CALLBACK_NAME_ProcessDEDrTerm,
          KIM_SUPPORT_STATUS_optional);
  ier = ier ||
      KIM_ModelComputeArgumentsCreate_SetCallbackSupportStatus(
          modelComputeArgumentsCreate,
          KIM_COMPUTE_CALLBACK_NAME_ProcessD2EDr2Term,
          KIM_SUPPORT_STATUS_optional);

  if (ier == TRUE){
    LOG_ERROR("Unable to set argument supportStatus.");
  }
  return ier;
}

/* compute_arguments destroy routine */
#include "KIM_ModelComputeArgumentsDestroyLogMacros.h"
static int compute_arguments_destroy(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy){
  /* Nothing to be done */
  return FALSE;
}

/* cutoff setup function */
static void compute_gr(double* const epsilon,
                       double* const C,
                       double* const Rzero,
                       double* const cutoff,
                       double* const rstar,
                       double g[GORDER])
{
  double const s = *rstar - (*cutoff);
  double const onebys = 1.0/s;
  double const onebyscube = onebys*onebys*onebys;

  /* compute values at r=rstar */
  double phi, dphi, d2phi;
  calc_phi_d2phi(epsilon,
                 C,
                 Rzero,
                 g,
                 cutoff,
                 cutoff,
                 *rstar,
                 &phi, &dphi, &d2phi);

  /* compute coefficients so that phi is C2 at rstar */
  g[0] = onebyscube*phi;

  g[1] = onebyscube*(dphi - 3.0*s*s*g[0]);

  g[2] = onebyscube*(d2phi - 6.0*s*(g[0] + s*g[1]));
}