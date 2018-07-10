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
 * Copyright (c) 2013--2016, Regents of the University of Minnesota.
 * All rights reserved.
 *
 * Contributors:
 *    Ryan S. Elliott
 *    Andrew Akerson
 *    Ellad B. Tadmor
 *    Valeriu Smirichinski
 *    Stephen M. Whalen
 *
 */


/*******************************************************************************
 *
 *  MorseEIP_GuthikondaElliott_2009
 *
 *  Temperature-dependent Morse pair potential KIM Model Driver
 *  Shifted to have zero energy at the cutoff radius
 *
 *  Language: C
 *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include "KIM_ModelDriverHeaders.h"

#define DIM 3  /* dimensionality of space */
#define MAXLINE 1024  /* max characters in line */
#define TRUE 1
#define FALSE 0
#define ONE 1.0

/* Define prototypes for Model Driver create and unit conversion function */
/**/
int model_driver_create(KIM_ModelDriverCreate * const modelDriverCreate,
  KIM_LengthUnit const requestedLengthUnit,
  KIM_EnergyUnit const requestedEnergyUnit,
  KIM_ChargeUnit const requestedChargeUnit,
  KIM_TemperatureUnit const requestedTemperatureUnit,
  KIM_TimeUnit const requestedTimeUnit);

int ConvertUnits(
  KIM_ModelDriverCreate * const modelDriverCreate,
  KIM_LengthUnit const requestedLengthUnit,
  KIM_EnergyUnit const requestedEnergyUnit,
  KIM_ChargeUnit const requestedChargeUnit,
  KIM_TemperatureUnit const requestedTemperatureUnit,
  KIM_TimeUnit const requestedTimeUnit,
  int numberUniqueSpeciesPairs,
  double * const cutoffs,
  double * const A1s,
  double * const A2s,
  double * const r1s,
  double * const r2s);

/* Define prototypes for Model (Driver) refresh and destroy */
/* defined as static to avoid namespace clashes with other Models    */
/**/
static int refresh(KIM_ModelRefresh * const modelRefresh);
static int destroy(KIM_ModelDestroy * const modelDestroy);

/* Define prototypes for compute, compute_arguments_create, and */
/* compute_arguments_destroy */
/**/
static int compute(KIM_ModelCompute const * const modelCompute,
  KIM_ModelComputeArguments const * const modelComputeArguments);
static int compute_arguments_create(
    KIM_ModelCompute const * const modelCompute,
        KIM_ModelComputeArgumentsCreate * const modelComputeArgumentsCreate);
static int compute_arguments_destroy(
    KIM_ModelCompute const * const modelCompute,
        KIM_ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy);

/* Define prototypes for internal functions used in this driver */
static void calc_phi(double const* const A,
                     double const* const B,
                     double const* const rHat,
                     double const* const shift,
                     double const* const cutoff,
                     double const r,
                     double* const phi);

static void calc_phi_dphi(double const* const A,
                          double const* const B,
                          double const* const rHat,
                          double const* const shift,
                          double const* const cutoff,
                          double const r,
                          double* const phi,
                          double* const dphi);

static void calc_phi_d2phi(double const* const A,
                           double const* const B,
                           double const* const rHat,
                           double const* const shift,
                           double const* const cutoff,
                           double const r,
                           double* const phi,
                           double* const dphi,
                           double* const d2phi);

double calc_A(double const* const A1,
              double const* const A2,
              double const* const A3,
              double const theta);

double calc_B(double const* const B1,
              double const* const B2,
              double const* const B3,
              double const theta);

double calc_rHat(double const* const r1,
                 double const* const r2,
                 double const* const r3,
                 double const theta);

double** AllocateAndInitialize2DArray(int const extentZero,
                                      int const extentOne);

double* AllocateAndInitialize1DArray(int const numberModelSpecies);
void Deallocate2DArrays(int const n, ...);
void Deallocate1DArrays(int const n, ...);
void getNextDataLine(FILE* const filePtr, char* nextLinePtr,
                     int const maxSize, int* const endOfFileFlag);

/* Define model_buffer structure */
struct model_buffer {
  int paddingNeighborHints;
  int halfListHints;

  double influenceDistance;
  double* cutoffs;

  int numberModelSpecies;

  double temperature;

  double* A1s;
  double* A2s;
  double* A3s;
  double* B1s;
  double* B2s;
  double* B3s;
  double* r1s;
  double* r2s;
  double* r3s;

  double** cutsq2D;
  double** As2D;
  double** Bs2D;
  double** rHats2D;
  double** shifts2D;
};

/* Calculate Guthikonda Elliott Paramaters */
double calc_A(double const* const A1,
              double const* const A2,
              double const* const A3,
              double const theta)
{
  double A;

  A = *A1 + (*A2)*(pow(theta, *A3) - 1.0);
  return A;
}

double calc_B(double const* const B1,
              double const* const B2,
              double const* const B3,
              double const theta)
{
  double B;

  B = *B1 + (*B2)*(pow(theta, *B3) - 1.0);
  return B;
}

double calc_rHat(double const* const r1,
                 double const* const r2,
                 double const* const r3,
                 double const theta)
{
  double rHat;

  rHat = *r1 + (*r2)*(exp((*r3)*(theta - 1.0)) - 1.0);
  return rHat;
}

/* Calculate pair potential phi(r) */
static void calc_phi(double const* const A,
                     double const* const B,
                     double const* const rHat,
                     double const* const shift,
                     double const* const cutoff,
                     double const r,
                     double* const phi)
{
  /* local variables */
  double ep;
  double ep2;
  double epsilon;
  double C;

  epsilon = -(*A);
  C = (*B)/(*rHat);

  ep  = exp(-(C)*(r-*rHat));
  ep2 = ep*ep;

  if (r > *cutoff)
  {
    /* Argument exceeds cutoff radius */
    *phi = 0.0;
  }
  else
  {
    *phi   = (epsilon)*( -ep2 + 2.0*ep ) + *shift;
  }

  return;
}

/* Calculate pair potential phi(r) and its derivative dphi(r) */
static void calc_phi_dphi(double const* const A,
                          double const* const B,
                          double const* const rHat,
                          double const* const shift,
                          double const* const cutoff,
                          double const r,
                          double* const phi,
                          double* const dphi)
{
  /* local variables */
  double ep;
  double ep2;
  double epsilon;
  double C;

  epsilon = -(*A);
  C = (*B)/(*rHat);

  ep  = exp(-(C)*(r-*rHat));
  ep2 = ep*ep;

  if (r > *cutoff)
  {
    /* Argument exceeds cutoff radius */
    *phi  = 0.0;
    *dphi = 0.0;
  }
  else
  {
    *phi  = (epsilon)*( -ep2 + 2.0*ep ) + *shift;
    *dphi = 2.0*(epsilon)*(C)*( -ep + ep2 );
  }

  return;
}

/*
  Calculate pair potential phi(r) and its 1st & 2nd derivatives dphi(r),
  d2phi(r)
*/
static void calc_phi_d2phi(double const* const A,
                           double const* const B,
                           double const* const rHat,
                           double const* const shift,
                           double const* const cutoff,
                           double const r,
                           double* const phi,
                           double* const dphi,
                           double* const d2phi)
{
  /* local variables */
  double ep;
  double ep2;
  double epsilon;
  double C;

  epsilon = -(*A);
  C = (*B)/(*rHat);
  ep  = exp(-(C)*(r-*rHat));
  ep2 = ep*ep;

  if (r > *cutoff)
  {
    /* Argument exceeds cutoff radius */
    *phi   = 0.0;
    *dphi  = 0.0;
    *d2phi = 0.0;
  }
  else
  {
    *phi   = (epsilon)*( -ep2 + 2.0*ep ) + *shift;
    *dphi  = 2.0*(epsilon)*(C)*( -ep + ep2 );
    *d2phi = 2.0*(epsilon)*(C)*(C)*(ep - 2.0*ep2);
  }

  return;
}

/*****************************************************************************/
/* Compute function                                                          */
/*****************************************************************************/
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
  int iSpecies, jSpecies;
  int const *neighListOfCurrentAtom;
  struct model_buffer* buffer;
  int comp_energy;
  int comp_force;
  int comp_particleEnergy;
  int comp_process_dEdr;
  int comp_process_d2Edr2;

  int* nAtoms;
  int* particleSpeciesCodes;
  int* particleContributing;
  double ijcutoff;
  double** cutsq2D;
  double** As2D;
  double** Bs2D;
  double** rHats2D;
  double** shifts2D;
  double* coords;
  double* energy;
  double* force;
  double* particleEnergy;
  int numOfAtomNeigh;

  /* get buffer from KIM object */
  KIM_ModelCompute_GetModelBufferPointer(modelCompute, (void **) &buffer);

  /* unpack info from the buffer */
  cutsq2D = (buffer->cutsq2D);
  As2D = (buffer->As2D);
  Bs2D = (buffer->Bs2D);
  rHats2D = (buffer->rHats2D);
  shifts2D = (buffer->shifts2D);

  /*
    check to see if we have been asked to compute the forces, particleEnergy,
    process_dEdr, and process_d2Edr2
  */
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
  if (ier == TRUE)
  {
    LOG_ERROR("GetArgumentPointer");
    return ier;
  }

  comp_energy = (energy != NULL);
  comp_force = (force != NULL);
  comp_particleEnergy = (particleEnergy != NULL);

  ier =
    KIM_ModelComputeArguments_IsCallbackPresent(
      modelComputeArguments,
      KIM_COMPUTE_CALLBACK_NAME_ProcessDEDrTerm,
      &comp_process_dEdr)
    ||
    KIM_ModelComputeArguments_IsCallbackPresent(
      modelComputeArguments,
      KIM_COMPUTE_CALLBACK_NAME_ProcessD2EDr2Term,
      &comp_process_d2Edr2);
  if (ier == TRUE)
  {
    LOG_ERROR("IsCallbackPresent");
    return ier;
  }

  /* initialize potential energies and forces */
  if (comp_particleEnergy)
  {
    for (i = 0; i < *nAtoms; ++i)
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
    for (i = 0; i < *nAtoms; ++i)
    {
      for (k = 0; k < DIM; ++k)
      {
        force[i*DIM + k] = 0.0;
      }
    }
  }

  /* Compute energy and forces */

  /* loop over particles and compute energy and forces */
  for (i=0; i<*nAtoms; i++)
  {
    if (particleContributing[i])
    {
      ier = KIM_ModelComputeArguments_GetNeighborList(
          modelComputeArguments,
          0, i, &numOfAtomNeigh, &neighListOfCurrentAtom);
      if (ier)
      {
        /* some sort of problem, exit */
        LOG_ERROR("KIM_get_neigh");
        ier = TRUE;
        return ier;
      }

      iSpecies = particleSpeciesCodes[i];

      /* loop over the neighbors of atom i */
      for (jj = 0; jj < numOfAtomNeigh; ++ jj)
      {
        /* get neighbor ID */
        j = neighListOfCurrentAtom[jj];

        if (i < j)
        {
          jSpecies = particleSpeciesCodes[j];
          /* compute relative position vector and squared distance */
          Rsqij = 0.0;
          for (k = 0; k < DIM; ++k)
          {
            Rij[k] = coords[j*DIM + k] - coords[i*DIM + k];

            /* compute squared distance */
            Rsqij += Rij[k]*Rij[k];
          }

          /* compute energy and force */
          if (Rsqij < cutsq2D[iSpecies][jSpecies])  /* particles are interacting? */
          {
            ijcutoff = sqrt(cutsq2D[iSpecies][jSpecies]);
            R = sqrt(Rsqij);
            if (comp_process_d2Edr2)
            {
              /* compute pair potential and its derivatives */
              calc_phi_d2phi(&As2D[iSpecies][jSpecies],
                             &Bs2D[iSpecies][jSpecies],
                             &rHats2D[iSpecies][jSpecies],
                             &shifts2D[iSpecies][jSpecies],
                             &ijcutoff, R, &phi, &dphi, &d2phi);

              /* compute dEidr and d2Eidr */
              dEidr = 0.5*dphi;
              d2Eidr = 0.5*d2phi;
            }
            else if (comp_force || comp_process_dEdr)
            {
              /* compute pair potential and its derivative */
              calc_phi_dphi(&(As2D[iSpecies][jSpecies]),
                            &(Bs2D[iSpecies][jSpecies]),
                            &(rHats2D[iSpecies][jSpecies]),
                            &(shifts2D[iSpecies][jSpecies]),
                            &ijcutoff, R, &phi, &dphi);

              /* compute dEidr */
              dEidr = 0.5*dphi;
            }
            else
            {
              /* compute just pair potential */
              calc_phi(&(As2D[iSpecies][jSpecies]),
                       &(Bs2D[iSpecies][jSpecies]),
                       &(rHats2D[iSpecies][jSpecies]),
                       &(shifts2D[iSpecies][jSpecies]),
                       &ijcutoff, R, &phi);
            }

            /* contribution to energy */
            if (comp_particleEnergy)
            {
              particleEnergy[i] += 0.5*phi;
            }
            if (comp_energy)
            {
              /* Full mode -- add half v to total energy */
              *energy += 0.5*phi;
            }

            /* contribution to process_dEdr */
            if (comp_process_dEdr)
            {
              ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                      modelComputeArguments,
                      dEidr, R, pRij, i, j);
            }

            /* contribution to process_d2Edr2 */
            if (comp_process_d2Edr2)
            {
              R_pairs[0] = R_pairs[1] = R;
              Rij_pairs[0][0] = Rij_pairs[1][0] = Rij[0];
              Rij_pairs[0][1] = Rij_pairs[1][1] = Rij[1];
              Rij_pairs[0][2] = Rij_pairs[1][2] = Rij[2];
              i_pairs[0] = i_pairs[1] = i;
              j_pairs[0] = j_pairs[1] = j;

              ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                      modelComputeArguments,
                      d2Eidr, pR_pairs, pRij_pairs, pi_pairs, pj_pairs);
            }

            /* contribution to forces */
            if (comp_force)
            {
              for (k = 0; k < DIM; ++k)
              { /* accumulate force on atom i */
                force[i*DIM + k] += dEidr*Rij[k]/R;
                /* accumulate force on atom j */
                force[j*DIM + k] -= dEidr*Rij[k]/R;
              }
            }
          }
        } /* End effective half-list check (i < j) */
      }  /* loop on jj */
    }  /* Check on whether particle i is contributing */
  } /* Outer loop over all atoms */

  /* everything is great */
  ier = FALSE;
  return ier;
}


/*****************************************************************************/
/* Miscellaneous helper functions                                            */
/*****************************************************************************/
void getNextDataLine( FILE* const filePtr, char* nextLinePtr,
                      int const maxSize, int* const endOfFileFlag)
{
  *endOfFileFlag = 0;
  do
  {
    if(fgets(nextLinePtr, maxSize, filePtr) == NULL)
    {
      *endOfFileFlag = 1;
      break;
    }
    while ((nextLinePtr[0] == ' ' || nextLinePtr[0] == '\t') ||
           (nextLinePtr[0] == '\n' || nextLinePtr[0] == '\r' ))
    {
      nextLinePtr = (nextLinePtr + 1);
    }
  }
  while ((strncmp("#", nextLinePtr, 1) == 0) || (strlen(nextLinePtr) == 0));
}

/*****************************************************************************/
double**  AllocateAndInitialize2DArray(int const extentZero,
                                       int const extentOne)
{
  double** arrayPtr;
  int i, j;
  /* allocate memory and set pointers */
  arrayPtr = malloc(extentZero*sizeof(double*));
  arrayPtr[0] = malloc((extentZero * extentOne)*sizeof(double));
  for (i = 1; i < extentZero; ++i)
  {
    arrayPtr[i] = arrayPtr[i-1] + extentOne;
  }

  /* initialize */
  for (i = 0; i < extentZero; ++i)
  {
    for (j = 0; j < extentOne; ++j)
    {
      arrayPtr[i][j] = 0.0;
    }
  }
  return arrayPtr;
}

/*****************************************************************************/
double* AllocateAndInitialize1DArray(int const numberModelSpecies)
{ double* arrayPtr;
  int numberUniqueSpeciesPairs;

  /* allocate memory and set pointers */
  numberUniqueSpeciesPairs = ((numberModelSpecies+1)*numberModelSpecies)/2;
  arrayPtr = calloc(numberUniqueSpeciesPairs, sizeof(double));
  return arrayPtr;
}

/*****************************************************************************/
void Deallocate1DArrays(int const n, ...)
{
  int i;
  va_list pointerArgs;

  va_start(pointerArgs, n);
  for (i = 0; i < n; i++)
  {
    free(va_arg(pointerArgs, double*));
  }
  va_end(pointerArgs);
}

/*****************************************************************************/
void Deallocate2DArrays(int const n, ...)
{
  double** nextArray;
  int i;
  va_list doublePointerArgs;

  va_start(doublePointerArgs, n);
  for (i = 0; i < n; i++)
  {
    nextArray = va_arg(doublePointerArgs, double**);
    free(nextArray[0]);
    free(nextArray);
  }
  va_end(doublePointerArgs);
}

/*****************************************************************************/
/* Unit conversion function                                                  */
/*****************************************************************************/
#include "KIM_ModelDriverCreateLogMacros.h"
int ConvertUnits(
  KIM_ModelDriverCreate * const modelDriverCreate,
  KIM_LengthUnit const requestedLengthUnit,
  KIM_EnergyUnit const requestedEnergyUnit,
  KIM_ChargeUnit const requestedChargeUnit,
  KIM_TemperatureUnit const requestedTemperatureUnit,
  KIM_TimeUnit const requestedTimeUnit,
  int numberUniqueSpeciesPairs,
  double * const cutoffs,
  double * const A1s,
  double * const A2s,
  double * const r1s,
  double * const r2s)
{
  int ier;
  int i;

  double convertLength = 1.0;
  double convertEnergy = 1.0;

  /* define default base units */
  KIM_LengthUnit const fromLength = KIM_LENGTH_UNIT_A;
  KIM_EnergyUnit const fromEnergy = KIM_ENERGY_UNIT_eV;
  KIM_ChargeUnit const fromCharge = KIM_CHARGE_UNIT_e;
  KIM_TemperatureUnit const fromTemperature = KIM_TEMPERATURE_UNIT_K;
  KIM_TimeUnit const fromTime = KIM_TIME_UNIT_ps;

  ier = KIM_ModelDriverCreate_ConvertUnit(modelDriverCreate,
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

  if (convertLength != ONE)
  {
    for (i=0; i < numberUniqueSpeciesPairs; ++i)
    {
      cutoffs[i] *= convertLength;
      r1s[i] *= convertLength;
      r2s[i] *= convertLength;
    }
  }

  /* Energy */
  ier = KIM_ModelDriverCreate_ConvertUnit(modelDriverCreate,
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
    for (i=0; i < numberUniqueSpeciesPairs; ++i)
    {
      A1s[i] *= convertEnergy;
      A2s[i] *= convertEnergy;
    }
  }

  /* register units */
  ier = KIM_ModelDriverCreate_SetUnits(modelDriverCreate,
      requestedLengthUnit,
      requestedEnergyUnit,
      KIM_CHARGE_UNIT_unused,
      KIM_TEMPERATURE_UNIT_unused,
      requestedTimeUnit);
  if (ier)
  {
    LOG_ERROR("Unable to set units to requested values");
    return ier;
  }

  ier = FALSE;
  return ier;
}


/*****************************************************************************/
/* Refresh function                                                          */
/*****************************************************************************/
static int refresh(KIM_ModelRefresh * const modelRefresh)
{
  /* Local variables */
  int ier;
  double influenceDistance;
  double numberModelSpecies;
  double *cutoffs;
  double *A1s, *A2s, *A3s;
  double *B1s, *B2s, *B3s;
  double *r1s, *r2s, *r3s;
  double temperature;
  double theta;

  double nextShift;
  double dummy;
  int indx;
  int i, j;

  struct model_buffer* buffer;

  /* get buffer from KIM object */
  KIM_ModelRefresh_GetModelBufferPointer(modelRefresh, (void **) &buffer);

  /* set value for 2D parameters */
  numberModelSpecies = buffer->numberModelSpecies;
  cutoffs = buffer->cutoffs;
  temperature = buffer->temperature;
  A1s = buffer->A1s;
  A2s = buffer->A2s;
  A3s = buffer->A3s;
  B1s = buffer->B1s;
  B2s = buffer->B2s;
  B3s = buffer->B3s;
  r1s = buffer->r1s;
  r2s = buffer->r2s;
  r3s = buffer->r3s;

  theta = temperature/333.15;
  influenceDistance = 0.0;
  for (i = 0; i < numberModelSpecies; ++i)
  {
    for (j = i; j < numberModelSpecies ; ++j)
    {
      indx = i*numberModelSpecies + j - (i*i + i)/2;

      if (influenceDistance < cutoffs[indx])
      {
        influenceDistance = cutoffs[indx];
      }

      buffer->cutsq2D[i][j] = buffer->cutsq2D[j][i]
          = (cutoffs[indx]*cutoffs[indx]);

      buffer->As2D[i][j] = buffer->As2D[j][i]
          = calc_A(&A1s[indx],
                   &A2s[indx],
                   &A3s[indx],
                   theta);

      buffer->Bs2D[i][j] = buffer->Bs2D[j][i]
          = calc_B(&B1s[indx],
                   &B2s[indx],
                   &B3s[indx],
                   theta);

      buffer->rHats2D[i][j] = buffer->rHats2D[j][i]
          = calc_rHat(&r1s[indx],
                      &r2s[indx],
                      &r3s[indx],
                      theta);
    }
  }
  buffer->influenceDistance = influenceDistance;

  /* Set influence distance pointer and cutoffs pointer */
  KIM_ModelRefresh_SetInfluenceDistancePointer(modelRefresh,
    &(buffer->influenceDistance));
  KIM_ModelRefresh_SetNeighborListPointers(modelRefresh,
    1, /* Use only a single neighbor list */
    &(buffer->influenceDistance),
    (const int*) &(buffer->paddingNeighborHints),
    (const int*) &(buffer->halfListHints));

  /* Set Values for Shifts */
  dummy = 0.0;
  influenceDistance += 1.0;  /* add a bit to avoid rounding problem */
  for (i = 0 ; i < numberModelSpecies; i++)
  {
    for (j = i; j < numberModelSpecies; j++)
    {
      indx = i*numberModelSpecies + j - (i*i + i)/2;
      /* call calc_phi with r=cutoff and shift=0.0 */
      calc_phi(&(buffer->As2D[i][j]),
               &(buffer->Bs2D[i][j]),
               &(buffer->rHats2D[i][j]),
               &dummy,
               &influenceDistance, cutoffs[indx], &nextShift);
      /* set shift to -shift */
      buffer->shifts2D[i][j] = buffer->shifts2D[j][i] = -nextShift;
    }
  }

  ier = FALSE;
  return ier;
}

/*****************************************************************************/
/* Destroy function                                                          */
/*****************************************************************************/
static int destroy(KIM_ModelDestroy * const modelDestroy)
{
  /* Local variables */
  int ier;
  struct model_buffer *buffer;

  KIM_ModelDestroy_GetModelBufferPointer(modelDestroy, (void **) &buffer);

  /* destroy the buffer */
  Deallocate2DArrays(5,
                     buffer->cutsq2D,
                     buffer->As2D,
                     buffer->Bs2D,
                     buffer->rHats2D,
                     buffer->shifts2D);
  Deallocate1DArrays(10,
                     buffer->cutoffs,
                     buffer->A1s,
                     buffer->A2s,
                     buffer->A3s,
                     buffer->B1s,
                     buffer->B2s,
                     buffer->B3s,
                     buffer->r1s,
                     buffer->r2s,
                     buffer->r3s);
  free(buffer);

  ier = FALSE;
  return ier;
}

/*****************************************************************************/
/* Compute arguments create function                                         */
/*****************************************************************************/
#include "KIM_ModelComputeArgumentsCreateLogMacros.h"
static int compute_arguments_create(
  KIM_ModelCompute const * const modelCompute,
  KIM_ModelComputeArgumentsCreate * const modelComputeArgumentsCreate)
{
  int ier;
  /* register arguments */
  ier = KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
      modelComputeArgumentsCreate,
      KIM_COMPUTE_ARGUMENT_NAME_partialEnergy,
      KIM_SUPPORT_STATUS_optional)
      ||
      KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
          modelComputeArgumentsCreate,
          KIM_COMPUTE_ARGUMENT_NAME_partialParticleEnergy,
          KIM_SUPPORT_STATUS_optional)
      ||
      KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
          modelComputeArgumentsCreate,
          KIM_COMPUTE_ARGUMENT_NAME_partialForces,
          KIM_SUPPORT_STATUS_optional);

  /* register callbacks */
  ier = ier ||
    KIM_ModelComputeArgumentsCreate_SetCallbackSupportStatus(
      modelComputeArgumentsCreate,
      KIM_COMPUTE_CALLBACK_NAME_ProcessDEDrTerm,
      KIM_SUPPORT_STATUS_optional) ||
   KIM_ModelComputeArgumentsCreate_SetCallbackSupportStatus(
      modelComputeArgumentsCreate,
      KIM_COMPUTE_CALLBACK_NAME_ProcessD2EDr2Term,
      KIM_SUPPORT_STATUS_optional);

  if (ier == TRUE)
  {
    LOG_ERROR("Unable to set argument supportStatus.");
    return TRUE;
  }
  else
  {
    return FALSE;
  }
}

/*****************************************************************************/
/* Compute arguments destroy function                                        */
/*****************************************************************************/
/* compue arguments destroy routine */
#include "KIM_ModelComputeArgumentsDestroyLogMacros.h"
static int compute_arguments_destroy(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy)
{
  /* nothing to be done */

  return FALSE;
}

/*****************************************************************************/
/* Create  function                                                          */
/*****************************************************************************/
#include "KIM_ModelDriverCreateLogMacros.h"
int model_driver_create(KIM_ModelDriverCreate * const modelDriverCreate,
  KIM_LengthUnit const requestedLengthUnit,
  KIM_EnergyUnit const requestedEnergyUnit,
  KIM_ChargeUnit const requestedChargeUnit,
  KIM_TemperatureUnit const requestedTemperatureUnit,
  KIM_TimeUnit const requestedTimeUnit)
{
  /* KIM variables */
  int numberOfParameterFiles;
  char const * paramfile1name;

  /* Local variables */
  FILE* fid;
  KIM_SpeciesName speciesName;

  double influenceDistance;
  double *cutoffs;
  double *A1s, *A2s, *A3s;
  double *B1s, *B2s, *B3s;
  double *r1s, *r2s, *r3s;

  double theta;
  double dbldummy;
  double nextShift;

  int ier;
  int i, j;
  struct model_buffer* buffer;

  int numberModelSpecies;
  int numberUniqueSpeciesPairs;

  int endOfFileFlag;
  int indx;
  char spec1[MAXLINE], spec2[MAXLINE], nextLine[MAXLINE];
  char dummy[12];
  char *nextLinePtr;
  double initialTemp;
  double nextCutoff;
  double nextA1, nextA2, nextA3;
  double nextB1, nextB2, nextB3;
  double nextr1, nextr2, nextr3;

  nextLinePtr = nextLine;

  /* using fixed units */
  ier = KIM_ModelDriverCreate_SetUnits(modelDriverCreate,
                                       KIM_LENGTH_UNIT_A,
                                       KIM_ENERGY_UNIT_eV,
                                       KIM_CHARGE_UNIT_unused,
                                       KIM_TEMPERATURE_UNIT_unused,
                                       KIM_TIME_UNIT_unused);
  if (ier == TRUE)
  {
    LOG_ERROR("Problem setting units");
    return ier;
  }

  ier = KIM_ModelDriverCreate_SetModelNumbering(modelDriverCreate,
                                                KIM_NUMBERING_zeroBased);
  if (ier == TRUE)
  {
    LOG_ERROR("Unable to set numbering");
    return ier;
  }

  /* store pointer to functions in KIM object */
  KIM_ModelDriverCreate_SetDestroyPointer(modelDriverCreate,
                                          KIM_LANGUAGE_NAME_c,
                                          (func *) destroy);
  KIM_ModelDriverCreate_SetComputeArgumentsCreatePointer(modelDriverCreate,
                                          KIM_LANGUAGE_NAME_c,
                                          (func *) compute_arguments_create);
  KIM_ModelDriverCreate_SetComputeArgumentsDestroyPointer(modelDriverCreate,
                                          KIM_LANGUAGE_NAME_c,
                                          (func *) compute_arguments_destroy);
  KIM_ModelDriverCreate_SetComputePointer(modelDriverCreate,
                                          KIM_LANGUAGE_NAME_c,
                                          (func *) compute);
  KIM_ModelDriverCreate_SetRefreshPointer(
      modelDriverCreate, KIM_LANGUAGE_NAME_c, (func *) refresh);

  /* get number of parameter files */
  KIM_ModelDriverCreate_GetNumberOfParameterFiles(
      modelDriverCreate, &numberOfParameterFiles);
  if (numberOfParameterFiles != 1)
  {
    ier = TRUE;
    LOG_ERROR("Incorrect number of parameter files. This driver expects only "
              "one Model parameter file.");
    return ier;
  }

  /* get parameter file name */
  ier = KIM_ModelDriverCreate_GetParameterFileName(
      modelDriverCreate,
      0,
      &paramfile1name);
  if (ier == TRUE)
  {
    LOG_ERROR("Unable to get Model parameter file name.");
    return ier;
  }

  fid = fopen(paramfile1name, "r");
  if (fid == NULL)
  {
    ier = TRUE;
    LOG_ERROR("Unable to open parameter file for Morse parameters");
    return ier;
  }

  /* Read line 0 of parameter file */
  getNextDataLine(fid, nextLinePtr, MAXLINE, &endOfFileFlag);
  ier = sscanf(nextLine, "%d %lg", &numberModelSpecies, &initialTemp);
  if (ier != 2)
  {
    sprintf(nextLine, "unable to read first line of the parameter file");
    ier = TRUE;
    LOG_ERROR(nextLine);
    fclose(fid);
    return ier;
  }

  /* Allocate memory used to store parameters read in */
  cutoffs = AllocateAndInitialize1DArray(numberModelSpecies);
  A1s = AllocateAndInitialize1DArray(numberModelSpecies);
  A2s = AllocateAndInitialize1DArray(numberModelSpecies);
  A3s = AllocateAndInitialize1DArray(numberModelSpecies);
  B1s = AllocateAndInitialize1DArray(numberModelSpecies);
  B2s = AllocateAndInitialize1DArray(numberModelSpecies);
  B3s = AllocateAndInitialize1DArray(numberModelSpecies);
  r1s = AllocateAndInitialize1DArray(numberModelSpecies);
  r2s = AllocateAndInitialize1DArray(numberModelSpecies);
  r3s = AllocateAndInitialize1DArray(numberModelSpecies);

  numberUniqueSpeciesPairs = ((numberModelSpecies+1)*numberModelSpecies)/2;
  /* set all values of cutoffs1D to -1 for check later */
  for (i = 0; i<numberUniqueSpeciesPairs; i++)
  {
    cutoffs[i] = -1.0;
  }

  /* Read all pure species lines. The species codes used corresponding to */
  /* the order in which they are read (starting with 0) */
  for (i = 0; i<numberModelSpecies; i++)
  {
    getNextDataLine(fid, nextLinePtr, MAXLINE, &endOfFileFlag);
    ier = sscanf(nextLine, "%s  %s  %lf  %lf  %lf  %lf",
                   spec1, spec2, &nextCutoff, &nextA1, &nextA2, &nextA3);
    if (ier != 6)
    {
      ier = TRUE;
      sprintf(nextLine, "error reading lines of the parameter file");
      LOG_ERROR(nextLine);
      fclose(fid);
      Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                         r1s, r2s, r3s);
      return ier;
    }

    if (strcmp(spec1, spec2) != 0)
    {
      ier = TRUE;
      sprintf(nextLine, "not all pure species interactions were given at "
        "the top of the parameter file");
      LOG_ERROR(nextLine);
      fclose(fid);
      Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                         r1s, r2s, r3s);
      return ier;
    }
    else
    {
      /* Register species code in API */
      speciesName = KIM_SpeciesName_FromString(spec1);
      ier = KIM_ModelDriverCreate_SetSpeciesCode(modelDriverCreate,
                                                        speciesName,
                                                        i);
      if (ier == TRUE)
      {
        LOG_ERROR("Unable to set species code");
        fclose(fid);
        Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                           r1s, r2s, r3s);
        return ier;
      }

      /* Calculate the appropriate array indices to store these params at */
      indx = i*numberModelSpecies + i - (i*i + i)/2;

      /* Store this line of parameters in the appropriate places */
      cutoffs[indx] = nextCutoff;
      A1s[indx] = nextA1;
      A2s[indx] = nextA2;
      A3s[indx] = nextA3;

      /* Now read the second line of parameters in this block and store */
      getNextDataLine(fid, nextLinePtr, MAXLINE, &endOfFileFlag);
      ier = sscanf(nextLine, "%lf %lf %lf",
                   &nextB1, &nextB2, &nextB3);
      if (ier != 3)
      {
        ier = TRUE;
        sprintf(nextLine, "error reading lines of the parameter file");
        LOG_ERROR(nextLine);
        fclose(fid);
        Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                           r1s, r2s, r3s);
        return ier;
      }
      B1s[indx] = nextB1;
      B2s[indx] = nextB2;
      B3s[indx] = nextB3;

      /* Now read the third line of parameters in this block and store */
      getNextDataLine(fid, nextLinePtr, MAXLINE, &endOfFileFlag);
      ier = sscanf(nextLine, "%lf %lf %lf",
                   &nextr1, &nextr2, &nextr3);
      if (ier != 3)
      {
        ier = TRUE;
        sprintf(nextLine, "error reading lines of the parameter file");
        LOG_ERROR(nextLine);
        fclose(fid);
        Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                           r1s, r2s, r3s);
        return ier;
      }
      r1s[indx] = nextr1;
      r2s[indx] = nextr2;
      r3s[indx] = nextr3;
    } /* End check on whether species listed are the same */
  } /* End loop over number of species */


  /* Now read all of the interspecies parameters, assuming they are */
  /* ordered correctly:                                             */
  /* (00, 01, ... ,0N, 10, 11, ... 1N, ..., N0, N1, ... NN)         */
  for (i = 0; i<numberModelSpecies; i++)
  {
    for (j = i+1; j < numberModelSpecies; j++)
    {
      getNextDataLine(fid, nextLinePtr, MAXLINE, &endOfFileFlag);
      ier = sscanf(nextLine, "%s  %s  %lf  %lf  %lf  %lf",
                     spec1, spec2, &nextCutoff, &nextA1, &nextA2, &nextA3);
      if (ier != 6)
      {
        ier = TRUE;
        sprintf(nextLine, "error reading lines of the parameter file");
        LOG_ERROR(nextLine);
        fclose(fid);
        Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                           r1s, r2s, r3s);
        return ier;
      }

      if (strcmp(spec1, spec2) == 0)
      {
        ier = TRUE;
        sprintf(nextLine, "more than one listing found for pure species "
          "interactions for species code %s", spec1);
        LOG_ERROR(nextLine);
        fclose(fid);
        Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                           r1s, r2s, r3s);
        return ier;
      }
      else
      {
        /* Calculate the appropriate array indices to store these params at */
        indx = i*numberModelSpecies + j - (i*i + i)/2;

        /* Store this line of parameters in the appropriate places */
        cutoffs[indx] = nextCutoff;
        A1s[indx] = nextA1;
        A2s[indx] = nextA2;
        A3s[indx] = nextA3;

        /* Now read the second line of parameters in this block and store */
        getNextDataLine(fid, nextLinePtr, MAXLINE, &endOfFileFlag);
        ier = sscanf(nextLine, "%lf %lf %lf",
                     &nextB1, &nextB2, &nextB3);
        if (ier != 3)
        {
          ier = TRUE;
          sprintf(nextLine, "error reading lines of the parameter file");
          LOG_ERROR(nextLine);
          fclose(fid);
          Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                             r1s, r2s, r3s);
          return ier;
        }
        B1s[indx] = nextB1;
        B2s[indx] = nextB2;
        B3s[indx] = nextB3;

        /* Now read the third line of parameters in this block and store */
        getNextDataLine(fid, nextLinePtr, MAXLINE, &endOfFileFlag);
        ier = sscanf(nextLine, "%lf %lf %lf",
                     &nextr1, &nextr2, &nextr3);
        if (ier != 3)
        {
          ier = TRUE;
          sprintf(nextLine, "error reading lines of the parameter file");
          LOG_ERROR(nextLine);
          fclose(fid);
          Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                             r1s, r2s, r3s);
          return ier;
        }
        r1s[indx] = nextr1;
        r2s[indx] = nextr2;
        r3s[indx] = nextr3;

      } /* End check on whether species listed are the same */
    } /* End inner loop over unique species combinations */
  } /* End outer loop over number of species */

  /* close parameter file */
  fclose(fid);

  /* Check that we got parameters for all pairs */
  ier = FALSE;
  sprintf(nextLine, "There are not values for the following pairs: \n");
  for (i = 0; i<numberModelSpecies; i++)
  {
    for (j = i; j < numberModelSpecies; j++)
    {
      indx = i*numberModelSpecies + j - (i*i + i)/2;
      if (cutoffs[indx] == -1.0)
      {
        sprintf(dummy, "%d and %d\n", i, j);
        strcat(nextLine, dummy);
        ier = TRUE;
      }
    }
  }
  if (ier == TRUE)
  {
    LOG_ERROR(nextLine);
    Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                       r1s, r2s, r3s);
    return ier;
  }

  /* convert parameters to appropriate units (in-place) */
  ier = ConvertUnits(
          modelDriverCreate,
          requestedLengthUnit,
          requestedEnergyUnit,
          requestedChargeUnit,
          requestedTemperatureUnit,
          requestedTimeUnit,
          numberUniqueSpeciesPairs,
          cutoffs,
          A1s,
          A2s,
          r1s,
          r2s);
  if (ier == TRUE)
  {
    sprintf(nextLine, "failed to convert units");
    LOG_ERROR(nextLine);
    Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                       r1s, r2s, r3s);
    return ier;
  }

  /* allocate buffer */
  buffer = (struct model_buffer*) malloc(sizeof(struct model_buffer));
  if (NULL == buffer)
  {
    ier = TRUE;
    LOG_ERROR("Coul not allocate memory for Model buffer");
    Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                       r1s, r2s, r3s);
    return ier;
  }

  /* register model buffer */
  KIM_ModelDriverCreate_SetModelBufferPointer(modelDriverCreate, (void*) buffer);

  /* Set number of species in buffer */
  buffer->numberModelSpecies = numberModelSpecies;

  /* Allocate arrays within buffer */
  buffer->cutsq2D
      = AllocateAndInitialize2DArray(numberModelSpecies, numberModelSpecies);
  buffer->As2D
      = AllocateAndInitialize2DArray(numberModelSpecies, numberModelSpecies);
  buffer->Bs2D
      = AllocateAndInitialize2DArray(numberModelSpecies, numberModelSpecies);
  buffer->rHats2D
      = AllocateAndInitialize2DArray(numberModelSpecies, numberModelSpecies);
  buffer->shifts2D
      = AllocateAndInitialize2DArray(numberModelSpecies, numberModelSpecies);

  /* set buffer parameters to those read in */
  buffer->temperature = initialTemp;
  buffer->cutoffs = cutoffs;
  buffer->A1s = A1s;
  buffer->A2s = A2s;
  buffer->A3s = A3s;
  buffer->B1s = B1s;
  buffer->B2s = B2s;
  buffer->B3s = B3s;
  buffer->r1s = r1s;
  buffer->r2s = r2s;
  buffer->r3s = r3s;

  /* Request that simulator not provide neighbors of padding atoms */
  buffer->paddingNeighborHints = 1;

  /* Request half lists from the simulator */
  buffer->halfListHints = 1;

  /* Register params for mutability by simulator */
  ier =
    KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate,
      1,
      &(buffer->temperature), "temperature")
    ||
    KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate,
      numberUniqueSpeciesPairs,
      buffer->cutoffs, "cutoffs")
    ||
    KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate,
      numberUniqueSpeciesPairs,
      buffer->A1s, "A1s")
    ||
    KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate,
      numberUniqueSpeciesPairs,
      buffer->A2s, "A2s")
    ||
    KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate,
      numberUniqueSpeciesPairs,
      buffer->A3s, "A3s")
    ||
    KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate,
      numberUniqueSpeciesPairs,
      buffer->B1s, "B1s")
    ||
    KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate,
      numberUniqueSpeciesPairs,
      buffer->B2s, "B2s")
    ||
    KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate,
      numberUniqueSpeciesPairs,
      buffer->B3s, "B3s")
    ||
    KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate,
      numberUniqueSpeciesPairs,
      buffer->r1s, "r1s")
    ||
    KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate,
      numberUniqueSpeciesPairs,
      buffer->r2s, "r2s")
    ||
    KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate,
      numberUniqueSpeciesPairs,
      buffer->r3s, "r3s");
  if (ier == TRUE)
  {
    LOG_ERROR("Could not register parameter pointers");
    Deallocate1DArrays(10, cutoffs, A1s, A2s, A3s, B1s, B2s, B3s,
                       r1s, r2s, r3s);
    Deallocate2DArrays(5, buffer->cutsq2D, buffer->As2D, buffer->Bs2D,
                       buffer->rHats2D, buffer->shifts2D);
    free(buffer);
    return ier;
  }

  /* Do some processing to set up 2D arrays (same stuff as is done in */
  /* refresh()                                                        */
  theta = buffer->temperature/333.15;
  influenceDistance = 0.0;
  for (i = 0; i < numberModelSpecies; ++i)
  {
    for (j = i; j < numberModelSpecies ; ++j)
    {
      indx = i*numberModelSpecies + j - (i*i + i)/2;

      if (influenceDistance < cutoffs[indx])
      {
        influenceDistance = cutoffs[indx];
      }

      buffer->cutsq2D[i][j] = buffer->cutsq2D[j][i]
          = (cutoffs[indx]*cutoffs[indx]);

      buffer->As2D[i][j] = buffer->As2D[j][i]
          = calc_A(&A1s[indx],
                   &A2s[indx],
                   &A3s[indx],
                   theta);

      buffer->Bs2D[i][j] = buffer->Bs2D[j][i]
          = calc_B(&B1s[indx],
                   &B2s[indx],
                   &B3s[indx],
                   theta);

      buffer->rHats2D[i][j] = buffer->rHats2D[j][i]
          = calc_rHat(&r1s[indx],
                      &r2s[indx],
                      &r3s[indx],
                      theta);
    }
  }
  /* Set influence distance */
  buffer->influenceDistance = influenceDistance;

  /* Set influence distance pointer and cutoff list pointer*/
  KIM_ModelDriverCreate_SetInfluenceDistancePointer(modelDriverCreate,
    &(buffer->influenceDistance));
  KIM_ModelDriverCreate_SetNeighborListPointers(modelDriverCreate,
    1, /* Use only a single neighbor list */
    &(buffer->influenceDistance),
    (const int*) &(buffer->paddingNeighborHints),
    (const int*) &(buffer->halfListHints));

  /* Set Values for Shifts */
  dbldummy = 0.0;
  influenceDistance += 1.0;  /* add a bit to avoid rounding problem */
  for (i = 0 ; i < numberModelSpecies; i++)
  {
    for (j = i; j < numberModelSpecies; j++)
    {
      indx = i*numberModelSpecies + j - (i*i + i)/2;
      /* call calc_phi with r=cutoff and shift=0.0 */
      calc_phi(&(buffer->As2D[i][j]),
               &(buffer->Bs2D[i][j]),
               &(buffer->rHats2D[i][j]),
               &dbldummy,
               &influenceDistance, cutoffs[indx], &nextShift);
      /* set shift to -shift */
      buffer->shifts2D[i][j] = buffer->shifts2D[j][i] = -nextShift;
    }
  }

  ier = FALSE;
  return ier;
}
