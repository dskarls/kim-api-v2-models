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
*/

/*
*
* Author: Daniel S. Karls (implementation adapted from Martin Z. Bazant's original code)
*
* Copyright (c) 2012-2018, Regents of the University of Minnesota.  All rights reserved.
*
*/

/*******************************************************************************
*
*  EDIP (Environmentally-Dependent Interatomic Potential) for Silicon
*
*  Language: C
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "KIM_ModelDriverHeaders.h"

#define TRUE 1
#define FALSE 0
#define ONE 1.0

#define MAX_NBRS 1024
#define DIM 3 /* dimensionality of the atomic coordinate space */
#define SPECCODE 1 /* internal species code */

/* A note on how this model is implemented here:

I adapted this model from Prof. Martin Z. Bazant's (bazant@math.mit.edu) code
(thanks also to Prof. Joao F. Justo for an additional FORTRAN reference code).
Because of the coordination-dependency of the 2-body and 3-body terms, we're
forced to do multiple loops over all of the atoms.  The first inner loop (the
"Pre-pass loop") goes over all of the other atoms and determines which ones are
inside the cutoff radius of the atom being considered by the outer loop.  For
those within the cutoff radius, it computes part of the 2-body interaction (the
part that doesn't require the coordination to be known) and the radial parts of
the 3-body interaction.  Finally, the coordination of each of the atoms within
the cutoff is added onto the total coordination, Z, for that atom.  Having Z,
we now have to loop over all other atoms once more in order to calculate the
coordination-dependent portion of the 2-body energy and get the final 2-body
energy for the atom being considered in the outermost loop, V2.  Note also that
we are unable to perform a "half-summation" for the two-body terms (with an
outer loop over i and an inner loop over j>i) because of the asymmetry the
coordination introduces into the 2-body energy (V2(R_ij,Z_i) != V2(R_ji,Z_j)).
We also need to use Z to perform another set of nested loops in order to
calculate the coordination-dependent portion of the 3-body energy, h(l,Z).
This will give us the total 3-body energy for the atom being considered in the
outermost loop, V3.  Since we are forced to loop over the same atoms multiple
times (first to get the coordination and then to compute the energies), it
makes sense to reduce the computational expense of the subsequent looping by
only looping over the atoms which have already been determined to fall within
the cutoff radius.  To this end, this code generates its own *internal*
neighbor list even if it is not provided by the Test.

The original EDIP publications can be found at:

     1.  M. Z. Bazant and E. Kaxiras, Phys. Rev. Lett. 77, 4370 (1996).
     2.  M. Z. Bazant, E. Kaxiras, J. F. Justo, Phys. Rev. B 56, 8542 (1997).
     3.  J. F. Justo, M. Z. Bazant, E. Kaxiras, V. V. Bulatov, and S. Yip,
     Phys. Rev. B 58, 2539 (1998).

If you are interested in an explanation of the looping structure, etc., in this
code and why it was chosen, I highly recommend looking at chapter 6 of Prof.
Bazant's 1997 Ph.D. thesis, available at

    http://web.mit.edu/bazant/www/thesis/md.html

Page 156 includes a pseudocode of what has been implemented here.

Written by Daniel S. Karls (University of Minnesota). Contact: karl0100 |AT| umn.edu.
*/

/* Model buffer definition */
struct model_buffer{
  double influenceDistance;
  double cutoff;
  double A;
  double B;
  double rh;
  double sig;
  double lam;
  double gam;
  double b;
  double c;
  double mu;
  double Qo;
  double eta;
  double bet;
  double alp;
  double u1;
  double u2;
  double u3;
  double u4;

  int paddingNeighborHints;
  int halfListHints;
};

/* ######## Function Prototypes ######## */

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




/* ######## Function Definitions ######## */

/* Primary compute function*/
#include "KIM_ModelComputeLogMacros.h"
static int compute(KIM_ModelCompute const * const modelCompute,
                   KIM_ModelComputeArguments const * const modelComputeArguments){
  /* Potential parameters */
  double A;
  double B;
  double rh;
  double a; /* Two-body interaction cutoff */
  double sig;
  double lam;
  double gam;
  double b; /* Three-body interaction cutoff (taken to be equal to parameter 'a' in this code) */
  double c;
  double mu;
  double Qo;
  double eta;
  double bet;
  double alp;
  double u1; /* 'u1', 'u2', 'u3', and 'u4' are parameters for the interpolation function tau(z) (Ismail & Kaxiras, 1993) */
  double u2;
  double u3;
  double u4;

  /* Local vars */
  struct model_buffer *buffer;
  int *numAtoms;
  int *particleSpeciesCodes;
  int *particleContributing;
  double *coords;
  double *energy;
  double *forces;
  double *virial;
  int const *neighListOfCurrentAtom;
  int numNeigh;
  int ier;
  int compute_energy, compute_forces, compute_virial;
  int i, j,k,jj,j2,j3,l,nk,nl;
  int num2[MAX_NBRS], num3[MAX_NBRS], numz[MAX_NBRS];
  int n2, n3, nz, dummycount;
  double V2, V3;
  double dx,dy,dz,r,rsqr,asqr;
  double rinv,rmainv,xinv,xinv3,den,Z,fZ;
  double dV2j,dV2ijx,dV2ijy,dV2ijz,pZ,dp;
  double temp0,temp1;
  double Qort,muhalf,u5;
  double rmbinv,winv,dwinv,tau,dtau,lcos,x,H, dHdx, dhdl;
  double dV3rij,dV3rijx,dV3rijy,dV3rijz;
  double dV3rik,dV3rikx,dV3riky,dV3rikz;
  double dV3l,dV3ljx,dV3ljy,dV3ljz,dV3lkx,dV3lky,dV3lkz;
  double dV2dZ,dxdZ,dV3dZ,dVdZ_sum;
  double dEdrl,dEdrlx,dEdrly,dEdrlz;
  double bmc,cmbinv;
  double fjx,fjy,fjz,fkx,fky,fkz;

  typedef struct{
    double t0, t1, t2, t3;
    double dx, dy, dz;
    double r;
  } store2;

  typedef struct{
    double g,dg;         /* 3-body radial function and its derivative */
    double rinv;         /* 1/r */
    double dx,dy,dz;     /* unit separation vector */
    double r;     /* bond length (only needed for virial) */
  } store3;

  typedef struct{
    double df;
    double dx,dy,dz;
    double r;
  } storez;

  store2 s2[MAX_NBRS];
  store3 s3[MAX_NBRS];
  storez sz[MAX_NBRS];

  /* Get model buffer */
  KIM_ModelCompute_GetModelBufferPointer(modelCompute, (void **) &buffer);

  /* Unpack parameters from buffer */
  a = buffer->cutoff;
  A = buffer->A;
  B = buffer->B;
  rh = buffer->rh;
  sig = buffer->sig;
  lam = buffer->lam;
  gam = buffer->gam;
  b = buffer->b;
  c = buffer->c;
  mu = buffer->mu;
  Qo = buffer->Qo;
  eta = buffer->eta;
  bet = buffer->bet;
  alp = buffer->alp;
  u1 = buffer->u1;
  u2 = buffer->u2;
  u3 = buffer->u3;
  u4 = buffer->u4;

  /* Combine Coefficients */
  asqr = a*a;
  Qort = sqrt(Qo);
  muhalf = mu*0.5;
  u5 = u2*u4;
  bmc = b-c;
  cmbinv = 1.0/(c-b);


  ier =
      KIM_ModelComputeArguments_GetArgumentPointerInteger(
          modelComputeArguments,
          KIM_COMPUTE_ARGUMENT_NAME_numberOfParticles,
          &numAtoms)
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
          &forces)
      ||
        KIM_ModelComputeArguments_GetArgumentPointerDouble(
          modelComputeArguments,
          KIM_COMPUTE_ARGUMENT_NAME_partialVirial,
          &virial);

  if (ier == TRUE){
    LOG_ERROR("Could not retrieve compute arguments");
    return ier;
  }

  /* Check to see what we've been asked to compute */
  compute_energy = (energy != 0);
  compute_forces = (forces != 0);
  compute_virial = (virial != 0);

  /* If the virial requested, we have to compute forces whether they were requested separately or not */
  if(compute_virial){
    compute_forces = 1;
  }

  /* Check to be sure that the species are correct */
  ier = TRUE;
  for(i = 0; i < *numAtoms; ++i){
    if (SPECCODE != particleSpeciesCodes[i]){
      LOG_ERROR("Unexpected species code detected");
      return ier;
    }
  }
  ier = FALSE;

  /* Initialize energy/forces/virial if they were requested */
  if(compute_energy){
    *energy = 0.0;
  }
  if(compute_forces){
    for(i = 0; i < *numAtoms; ++i){
      for(k = 0; k < DIM; ++k){
        forces[i*DIM + k] = 0.0;
      }
    }
  }
  if(compute_virial){
    for(i = 0; i < 6; i++){
      virial[i] = 0.0;
    }
  }

  V2 = 0.0;
  V3 = 0.0;

  /* --- Level 1: Outer loop over atoms --- */
  for(i = 0; i < *numAtoms; i++){

    if (particleContributing[i]){
      /* Get neighbor list for atom i */
      ier = KIM_ModelComputeArguments_GetNeighborList(modelComputeArguments,
          0, i, &numNeigh, &neighListOfCurrentAtom);
      if(ier){
        LOG_ERROR("Could not access neighbor list");
        ier = TRUE;
        return ier;
      }

      /* Reset coordination Z and neighbor numbers */
      Z = 0.0;
      n2 = 0;
      n3 = 0;
      nz = 0;

      for(dummycount = 0; dummycount < MAX_NBRS; dummycount++){
        numz[dummycount]=0;
      }

      /* --- Level 2: Prepass loop over neighbors of atom i --- */
      for (jj = 0; jj < numNeigh; ++jj){
        j = neighListOfCurrentAtom[jj]; /* get neighbor ID */

        if(j != i){
          dx = coords[j*DIM+0]-coords[i*DIM+0];
          dy = coords[j*DIM+1]-coords[i*DIM+1];
          dz = coords[j*DIM+2]-coords[i*DIM+2];

          if(dx < a){
            if(dy < a){
              if(dz < a){
                rsqr = dx*dx + dy*dy + dz*dz;
                if(rsqr < asqr){
                  num2[n2] = j;
                  r = sqrt(rsqr);
                  rinv=1.0/r;
                  dx *= rinv;
                  dy *= rinv;
                  dz *= rinv;
                  rmainv=1.0/(r - a);
                  s2[n2].t0 = A*exp(sig*rmainv);
                  s2[n2].t1 = pow(B*rinv,rh);
                  s2[n2].t2 = rh*rinv;
                  s2[n2].t3 = sig*rmainv*rmainv;
                  s2[n2].dx = dx;
                  s2[n2].dy = dy;
                  s2[n2].dz = dz;
                  s2[n2].r = r;
                  n2++;

                  /* Also compute the radial part of the 3-body energy term */
                  if(r < a){
                    num3[n3] = j;
                    rmbinv = 1.0/(r - a);
                    temp1 = gam*rmbinv;
                    temp0 = exp(temp1);
                    s3[n3].g = temp0;
                    s3[n3].dg = -rmbinv*temp1*temp0;
                    s3[n3].dx = dx;
                    s3[n3].dy = dy;
                    s3[n3].dz = dz;
                    s3[n3].rinv = rinv;
                    s3[n3].r = r;
                    n3++;

                    /* Coordination and neighbor function c < r < a */
                    if(r < b){
                      if(r <= c){
                        Z += 1.0;
                      }else{
                        xinv = bmc/(r - c);
                        xinv3 = xinv*xinv*xinv;
                        den = 1.0/(1 - xinv3);
                        temp1 = alp*den;
                        fZ = exp(temp1);
                        Z += fZ;
                        numz[nz] = j;
                        sz[nz].df = fZ*temp1*den*3.0*xinv3*xinv*cmbinv;   /* df/dr */
                        sz[nz].dx = dx;
                        sz[nz].dy = dy;
                        sz[nz].dz = dz;
                        sz[nz].r = r;
                        nz++;
                      }
                    } /* End if-statement of r < b */
                  } /* End if-statement r < a */
                }/*  End if-statement rsqrt < asqr */
              } /* End if-statement dz < a */
            } /* End if-statement dy < a */
          } /* End if-statement dx < a */
        } /* End if-statement j!=i */
      } /*  End loop over j */

      dVdZ_sum=0.0;
      temp0 = bet*Z;
      pZ = exp(-temp0*Z);
      dp = -2.0*temp0*pZ; /* derivative of bond order */

      /* --- Level 2: Second loop over pairs to get the final 2-body energy --- */
      for(j2=0; j2 < n2; j2++){
        temp0 = s2[j2].t1 - pZ;

        /* Two-body energy */
        V2 += temp0*s2[j2].t0;

        /* Two-body forces */
        if(compute_forces){
          dV2j = - (s2[j2].t0) * ((s2[j2].t1)*(s2[j2].t2) + temp0 * (s2[j2].t3)); /* dV2/dr */
          dV2ijx = dV2j * s2[j2].dx; /* x component of force on atom i */
          dV2ijy = dV2j * s2[j2].dy; /* y component of force on atom i */
          dV2ijz = dV2j * s2[j2].dz; /* z component of force on atom i */

          forces[i*DIM+0] += dV2ijx;
          forces[i*DIM+1] += dV2ijy;
          forces[i*DIM+2] += dV2ijz;

          j = num2[j2];
          forces[j*DIM+0] -= dV2ijx;
          forces[j*DIM+1] -= dV2ijy;
          forces[j*DIM+2] -= dV2ijz;

          /* Accumulation of pair coordination forces */
          dV2dZ = - dp * s2[j2].t0;
          dVdZ_sum += dV2dZ;
        }

        /* Two-body contribution to virial */
        if(compute_virial){
          virial[0] += s2[j2].r * dV2ijx * s2[j2].dx;
          virial[1] += s2[j2].r * dV2ijy * s2[j2].dy;
          virial[2] += s2[j2].r * dV2ijz * s2[j2].dz;
          virial[3] += s2[j2].r * dV2ijy * s2[j2].dz;
          virial[4] += s2[j2].r * dV2ijz * s2[j2].dx;
          virial[5] += s2[j2].r * dV2ijx * s2[j2].dy;
        }
      } /* End loop over j2 */

      /* Coordination-dependence of three-body interactions */
      winv = Qort*exp(-muhalf*Z); /* inverse width of angular function */
      dwinv = -muhalf*winv;       /* its derivative */
      temp0 = exp(-u4*Z);
      tau = u1+u2*temp0*(u3-temp0); /* -cosine of angular minimum */
      dtau = u5*temp0*(2*temp0-u3); /* its derivative */

      /* --- Level 2: First loop for three-body interactions --- */
      for(j3 = 0; j3 < (n3 - 1); j3++){
        j = num3[j3];

        /* --- Level 3: Second loop for three-body interactions --- */
        for(nk = j3 + 1; nk < n3; nk++){
          k = num3[nk];

          /* Angular function h(l,Z) */
          lcos = s3[j3].dx * s3[nk].dx + s3[j3].dy * s3[nk].dy + s3[j3].dz * s3[nk].dz;
          x = (lcos + tau)*winv;
          temp0 = exp(-x*x);

          H = lam*(1 - temp0 + eta*x*x);
          dHdx = 2*lam*x*(temp0 + eta);
          dhdl = dHdx*winv;

          /* Three-body energy */
          temp1 = s3[j3].g * s3[nk].g;
          V3 += temp1*H;

          if(compute_forces){
            /* (-) radial force on atom j */
            dV3rij = s3[j3].dg * s3[nk].g * H;
            dV3rijx = dV3rij * s3[j3].dx;
            dV3rijy = dV3rij * s3[j3].dy;
            dV3rijz = dV3rij * s3[j3].dz;
            fjx = dV3rijx;
            fjy = dV3rijy;
            fjz = dV3rijz;

            /* (-) radial force on atom k */
            dV3rik = s3[j3].g * s3[nk].dg * H;
            dV3rikx = dV3rik * s3[nk].dx;
            dV3riky = dV3rik * s3[nk].dy;
            dV3rikz = dV3rik * s3[nk].dz;
            fkx = dV3rikx;
            fky = dV3riky;
            fkz = dV3rikz;

            /* (-) angular force on j */
            dV3l = temp1*dhdl;
            dV3ljx = dV3l * (s3[nk].dx - lcos * s3[j3].dx) * s3[j3].rinv;
            dV3ljy = dV3l * (s3[nk].dy - lcos * s3[j3].dy) * s3[j3].rinv;
            dV3ljz = dV3l * (s3[nk].dz - lcos * s3[j3].dz) * s3[j3].rinv;
            fjx += dV3ljx;
            fjy += dV3ljy;
            fjz += dV3ljz;

            /* (-) angular force on k */
            dV3lkx = dV3l * (s3[j3].dx - lcos * s3[nk].dx) * s3[nk].rinv;
            dV3lky = dV3l * (s3[j3].dy - lcos * s3[nk].dy) * s3[nk].rinv;
            dV3lkz = dV3l * (s3[j3].dz - lcos * s3[nk].dz) * s3[nk].rinv;
            fkx += dV3lkx;
            fky += dV3lky;
            fkz += dV3lkz;

            /* Apply radial + angular forces to i, j, k */
            forces[j*DIM+0] -= fjx;
            forces[j*DIM+1] -= fjy;
            forces[j*DIM+2] -= fjz;
            forces[k*DIM+0] -= fkx;
            forces[k*DIM+1] -= fky;
            forces[k*DIM+2] -= fkz;
            forces[i*DIM+0] += fjx + fkx;
            forces[i*DIM+1] += fjy + fky;
            forces[i*DIM+2] += fjz + fkz;
          }

          /* dV3/dR contributions to virial */
          if(compute_virial){
            virial[0] += s3[j3].r * (fjx*s3[j3].dx) + s3[nk].r * (fkx*s3[nk].dx);
            virial[1] += s3[j3].r * (fjy*s3[j3].dy) + s3[nk].r * (fky*s3[nk].dy);
            virial[2] += s3[j3].r * (fjz*s3[j3].dz) + s3[nk].r * (fkz*s3[nk].dz);
            virial[3] += s3[j3].r * (fjy*s3[j3].dz) + s3[nk].r * (fky*s3[nk].dz);
            virial[4] += s3[j3].r * (fjz*s3[j3].dx) + s3[nk].r * (fkz*s3[nk].dx);
            virial[5] += s3[j3].r * (fjx*s3[j3].dy) + s3[nk].r * (fkx*s3[nk].dy);
          }

          dxdZ = dwinv*(lcos + tau) + winv*dtau;
          dV3dZ = temp1*dHdx*dxdZ;

          /* Accumulation of 3-body coordination forces */
          dVdZ_sum += dV3dZ;
        } /* End loop over k */
      } /* End loop over j3 */


      /* --- Level 2: Loop to apply coordination forces --- */
      if(compute_forces){
        for(nl = 0; nl < nz; nl++){
          dEdrl = dVdZ_sum * sz[nl].df;
          dEdrlx = (dEdrl * sz[nl].dx);
          dEdrly = (dEdrl * sz[nl].dy);
          dEdrlz = (dEdrl * sz[nl].dz);
          forces[i*DIM+0] += dEdrlx;
          forces[i*DIM+1] += dEdrly;
          forces[i*DIM+2] += dEdrlz;
          l = numz[nl];
          forces[l*DIM+0] -= dEdrlx;
          forces[l*DIM+1] -= dEdrly;
          forces[l*DIM+2] -= dEdrlz;

          /* dE/dZ*dZ/dr contributions to virial */
          if(compute_virial){ /*  Recall that, above, we set compute_forces to 1 if compute_virial is 1 so this nested if-statement is ok */
            virial[0] += sz[nl].r * (dEdrlx * sz[nl].dx);
            virial[1] += sz[nl].r * (dEdrly * sz[nl].dy);
            virial[2] += sz[nl].r * (dEdrlz * sz[nl].dz);
            virial[3] += sz[nl].r * (dEdrly * sz[nl].dz);
            virial[4] += sz[nl].r * (dEdrlz * sz[nl].dx);
            virial[5] += sz[nl].r * (dEdrlx * sz[nl].dy);
          }
        }
      }
    } /* Check on whether particle i is contributing */
  } /* End loop over i */

  if(compute_energy){
    *energy = V2+V3;
  }

  /* Return the status code indicating everything is fine */
  ier = FALSE;
  return ier;
}




/* Create function */
#include "KIM_ModelDriverCreateLogMacros.h"
int model_driver_create(KIM_ModelDriverCreate *const modelDriverCreate,
                        KIM_LengthUnit const requestedLengthUnit,
                        KIM_EnergyUnit const requestedEnergyUnit,
                        KIM_ChargeUnit const requestedChargeUnit,
                        KIM_TemperatureUnit const requestedTemperatureUnit,
                        KIM_TimeUnit const requestedTimeUnit){
  /* Create Model buffer */
  struct model_buffer *buffer;

  /* KIM variables */
  int numberOfParameterFiles;
  char const *paramfile1name; /* This driver only actually supports a single parameter file */
  char speciesNameString[100];
  KIM_SpeciesName speciesName;

  /* Model parameters */
  double cutoff;
  double A;
  double B;
  double rh;
  double a;
  double sig;
  double lam;
  double gam;
  double b;
  double c;
  double mu;
  double Qo;
  double eta;
  double bet;
  double alp;
  double u1;
  double u2;
  double u3;
  double u4;
  int ier;
  KIM_LengthUnit fromLength;
  KIM_EnergyUnit fromEnergy;
  KIM_ChargeUnit fromCharge;
  KIM_TemperatureUnit fromTemperature;
  KIM_TimeUnit fromTime;

  /* Miscellaneous */
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
    LOG_ERROR("Unable to open parameter file for EDIP parameters");
    return ier;
  }
  ier = fscanf(fid, "%s \n%lf \n%lf \n%lf \n%lf %lf \n%lf \n%lf \n%lf \n%lf \n%lf \n%lf %lf \n%lf \
                     \n%lf \n%lf %lf \n%lf \n%lf",
               speciesNameString,
               &a, /* cutoff in Angstroms */
               &A, &B, &rh, &sig, &lam, &gam, &b, &c, &mu, &Qo, &eta, &bet, &alp, &u1, &u2, &u3, &u4);
  fclose(fid);

  /* Check that we read the right number of parameters */
  if (19 != ier){
    ier = TRUE;
    LOG_ERROR("Unable to read all EDIP parameters");
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
    B *= convertLength;
    a *= convertLength;
    sig *= convertLength;
    gam *= convertLength;
    b *= convertLength;
    c *= convertLength;
  }

  if(convertEnergy != ONE){
    A *= convertEnergy;
    lam *= convertEnergy;
  }
  /*  NOTE: rh, mu, Qo, eta, bet, alp, u1, u2, u3, u4 are unitless so skip them*/

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
                                          KIM_LANGUAGE_NAME_c, (func *) &destroy);
  ier = ier || KIM_ModelDriverCreate_SetComputeArgumentsCreatePointer(modelDriverCreate,
                                          KIM_LANGUAGE_NAME_c, (func *) &compute_arguments_create);
  ier = ier || KIM_ModelDriverCreate_SetComputeArgumentsDestroyPointer(modelDriverCreate,
                                          KIM_LANGUAGE_NAME_c, (func *) &compute_arguments_destroy);
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
  buffer->influenceDistance = cutoff;
  buffer->cutoff = cutoff;
  buffer->A = A;
  buffer->B = B;
  buffer->rh = rh;
  buffer->sig = sig;
  buffer->lam = lam;
  buffer->gam = gam;
  buffer->b = b;
  buffer->c = c;
  buffer->mu = mu;
  buffer->Qo = Qo;
  buffer->eta = eta;
  buffer->bet = bet;
  buffer->alp = alp;
  buffer->u1 = u1;
  buffer->u2 = u2;
  buffer->u3 = u3;
  buffer->u4 = u4;

  /* Request that simulator omit neighbors of padding atoms */
  buffer->paddingNeighborHints = 1;

  /* Request a full list rather than a half list (still more efficient due to
   * the particular looping structure adopted from Bazant here) */
  buffer->halfListHints = 0;

  /* Register influence distance pointer */
  LOG_INFORMATION("Registering influence distance pointer");
  KIM_ModelDriverCreate_SetInfluenceDistancePointer(modelDriverCreate, &(buffer->influenceDistance));

  /* Register cutoff pointer */
  LOG_INFORMATION("Registering cutoff pointer");
  KIM_ModelDriverCreate_SetNeighborListPointers(modelDriverCreate, 1, &(buffer->cutoff),
      &(buffer->paddingNeighborHints), &(buffer->halfListHints));

  /* Register Model parameters */
  LOG_INFORMATION("Registering Model parameters");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->cutoff), "a");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->A), "A");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->B), "B");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->rh), "rh");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->sig), "sig");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->lam), "lam");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->gam), "gam");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->b), "b");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->c), "c");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->mu), "mu");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->Qo), "Qo");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->eta), "eta");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->bet), "bet");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->alp), "alp");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->u1), "u1");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->u2), "u2");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->u3), "u3");
  KIM_ModelDriverCreate_SetParameterPointerDouble(modelDriverCreate, 1, &(buffer->u4), "u4");

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
  struct model_buffer* buffer;

  /* Get model buffer from KIM object */
  LOG_INFORMATION("Retrieving Model buffer");
  KIM_ModelRefresh_GetModelBufferPointer(modelRefresh, (void **) &buffer);

  LOG_INFORMATION("Resetting influence distance and cutoff");
  KIM_ModelRefresh_SetInfluenceDistancePointer(
      modelRefresh, &(buffer->influenceDistance));
  KIM_ModelRefresh_SetNeighborListPointers(
      modelRefresh, 1, &(buffer->influenceDistance),
      &(buffer->paddingNeighborHints), &(buffer->halfListHints));

  return FALSE;
}




/* Destroy function */
#include "KIM_ModelDestroyLogMacros.h"
static int destroy(KIM_ModelDestroy * const modelDestroy){
  /* Local variables */
  struct model_buffer* buffer;
  int ier;

  /* Get model buffer from KIM object */
  KIM_ModelDestroy_GetModelBufferPointer(modelDestroy, (void **) &buffer);

  /* Free the buffer */
  LOG_INFORMATION("Deallocating Model buffer");
  free(buffer);

  ier = FALSE;
  return ier;
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
                                       KIM_SUPPORT_STATUS_notSupported);
  ier = ier || KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
                                       modelComputeArgumentsCreate,
                                       KIM_COMPUTE_ARGUMENT_NAME_partialForces,
                                       KIM_SUPPORT_STATUS_optional);
  ier = ier || KIM_ModelComputeArgumentsCreate_SetArgumentSupportStatus(
                                       modelComputeArgumentsCreate,
                                       KIM_COMPUTE_ARGUMENT_NAME_partialVirial,
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
