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
* Copyright (c) 2012, Regents of the University of Minnesota.  All rights reserved.
*
* Contributors:
*    Amit Singh
*/

/******************************************************************************
*
*  Four_Body_Mistriotis_Flytzanis_Farantos
*
*  Mistriotis-Flytzanis-Farantos Four-Body potential KIM Model Driver
*
*  Language: C
*
*  Release:
*
******************************************************************************/
/******************************************************************************

* For four-body potential any potential energy function for N-particle system can be written in the
* following form of two-body, three-body and four-body terms:
     F(1, 2, ...., N)   = sum_{i; 1 <= i \neq j <= N} 0.5*phi_two(r; r = r_ij)
                          + sum_{i; 1 <= i \neq j < k <= N} phi_three(r_ij, r_ik, r_jk);
                          + sum_{i; 1 <= i \neq j < k < l <= N} phi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl);

* For modified Stillinger-Weber proposed by "Mistriotis et al, Physical Review B 39, 1212 (1989)",
* the two-body term phi_two can be written as:

     phi_two(r) = epsilon * f_2(r_cap; r_cap = r/sigma); where f_2 is

     f_2(r_cap) = A * ( B*r_cap^(-p) - r_cap^-q ) * exp(1/(r_cap - a)); when  r_cap <  a, where a = cutoff/sigma,
                = 0   when r_cap >= a

* And three-body term phi_three can be written as:

     phi_three(r_ij, r_ik, r_jk) = epsilon * f_3(r1_cap, r2_cap, r3_cap); where

     r1_cap = r_ij/sigma, r2_cap = r_ik/sigma, r3_cap = r_jk/sigma,      and

     f_3(r1_cap, r2_cap, r3_cap) = lambda * (exp(gamma((1/(r1_cap - a)) + (1/(r2_cap - a)))))
                                 * {1 - exp[-Q((costheta_jik - costheta_0)^2)]};  when r1_cap < a && r2_cap < a
                                 = 0;                                             otherwise

     costheta_jik = (r_ij^2 + r_ik^2 - r_jk^2) / (2*r_ij*r_ik).

* and four-body term phi_four can be written as:

     phi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl) = epsilon * f_4(r1_cap, r2_cap, r3_cap, r4_cap, r5_cap, r6_cap); where

     r1_cap = r_ij/sigma, r2_cap = r_ik/sigma, r3_cap = r_il/sigma, r4_cap = r_jk/sigma, r5_cap = r_jl/sigma, r6_cap = r_kl/sigma    and

     f_4(r1_cap, r2_cap, r3_cap, r4_cap, r5_cap, r6_cap) = g(r1_cap, r2_cap, r3_cap, costheta_jik, costheta_jil, costheta_kil); where

     g(r1_cap, r2_cap, r3_cap, costheta_jik, costheta_jil, costheta_kil) = lambda_2 * (exp(gamma((1/(r1_cap - a)) + (1/(r2_cap - a)) + (1/(r3_cap - a)))))
                                                                         * {1 - exp[-Q((costheta_jik - costheta_0)^2 + (costheta_jil - costheta_0)^2 + (costheta_kil - costheta_0)^2)]}; when r1_cap < a && r2_cap < a && r3_cap < a
                                                                         = 0;                                                                                                       otherwise

     costheta_jik = (r_ij^2 + r_ik^2 - r_jk^2) / (2*r_ij*r_ik);
     costheta_jil = (r_ij^2 + r_il^2 - r_jl^2) / (2*r_ij*r_il);
     costheta_kil = (r_ik^2 + r_il^2 - r_kl^2) / (2*r_ik*r_il).

******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "KIM_ModelDriverHeaders.h"

#define TRUE 1
#define FALSE 0
#define ONE 1.0
#define DIM 3       /* dimensionality of space */
#define SPEC1 1  /* internal species code */
#define SPEC2 2  /* internal species code */

/* Define prototype for Model Driver initialization */
int model_driver_create(KIM_ModelDriverCreate *const modelDriverCreate,
                        KIM_LengthUnit const requestedLengthUnit,
                        KIM_EnergyUnit const requestedEnergyUnit,
                        KIM_ChargeUnit const requestedChargeUnit,
                        KIM_TemperatureUnit const requestedTemperatureUnit,
                        KIM_TimeUnit const requestedTimeUnit);

/* Define prototype for compute routine */
/* (defined as static to avoid namespace clashes with other codes) */
static int compute(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArguments const * const modelComputeArguments);
static int compute_arguments_create(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsCreate * const modelComputeArgumentsCreate);
static int compute_arguments_destroy(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy);

int ConvertUnits(
  KIM_ModelDriverCreate * const modelDriverCreate,
  KIM_LengthUnit const requestedLengthUnit,
  KIM_EnergyUnit const requestedEnergyUnit,
  KIM_ChargeUnit const requestedChargeUnit,
  KIM_TemperatureUnit const requestedTemperatureUnit,
  KIM_TimeUnit const requestedTimeUnit,
  int num_interactions,
  double * const A,
  double * const sigma,
  double * const lambda,
  double * const lambda_2,
  double * const epsilon);

/* Define prototype for refresh routine */
/* (defined as static to avoid namespace clashes with other codes) */
static int refresh(KIM_ModelRefresh * const modelRefresh);

/* Define prototype for destroy */
/* (defined as static to avoid namespace clashes with other codes) */
static int destroy(KIM_ModelDestroy *const modelDestroy);

static void calc_phi_two(double* A, double* B, double* p, double* q, double* a,
                         double* sigma, double* epsilon, double r, double* phi);
static void calc_phi_dphi_two(double* A, double* B, double* p, double* q,
                              double* a, double* sigma, double* epsilon,
                              double r, double* phi, double* dphi);
static void calc_phi_d2phi_two(double* A, double* B, double* p, double* q,
                               double* a, double* sigma, double* epsilon,
                               double r, double* phi, double* dphi,
                               double* d2phi);

static void calc_phi_three(double* a, double* lambda, double* gamma,
                           double* sigma, double* epsilon, double* Q,
                           double* costhetat, double rij, double rik,
                           double rjk, double* phi);
static void calc_phi_dphi_three(double* a, double* lambda, double* gamma,
                                double* sigma, double* epsilon, double* Q,
                                double* costhetat, double rij, double rik,
                                double rjk, double* phi, double* dphi);
static void calc_phi_d2phi_three(double* a, double* lambda, double* gamma,
                                 double* sigma, double* epsilon, double* Q,
                                 double* costhetat, double rij, double rik,
                                 double rjk, double* phi, double* dphi,
                                 double* d2phi);

static void calc_phi_four(double* a, double* lambda_2, double* gamma,
                          double* sigma, double* epsilon, double* Q,
                          double* costhetat, double rij, double rik, double ril,
                          double rjk, double rjl, double rkl, double* phi);
static void calc_phi_dphi_four(double* a, double* lambda_2, double* gamma,
                               double* sigma, double* epsilon, double* Q,
                               double* costhetat, double rij, double rik,
                               double ril, double rjk, double rjl, double rkl,
                               double* phi, double* dphi);
static void calc_phi_d2phi_four(double* a, double* lambda_2, double* gamma,
                                double* sigma, double* epsilon, double* Q,
                                double* costhetat, double rij, double rik,
                                double ril, double rjk, double rjl, double rkl,
                                double* phi, double* dphi, double* d2phi);

/* Define model_buffer structure */
struct model_buffer {
  double influenceDistance;
  double* cutsq;
  double* A;
  double* B;
  double* p;
  double* q;
  double* a;
  double* lambda;
  double* lambda_2;
  double* gamma;
  double* sigma;
  double* epsilon;
  double* Q;
  double* costhetat;

  double num_interactions;

  int paddingNeighborHints;
  int halfListHints;
};

/* Calculate two-body term phi_two(r) */
static void calc_phi_two(double* A, double* B, double* p, double* q, double* a,
                         double* sigma, double* epsilon, double r, double* phi)
{
  /* Local variables */
  double r_cap;

  r_cap = r/(*sigma);

  if (r_cap >=  *a)
  {
    *phi = 0.0;
  }
  else
  {
    *phi = (*epsilon) * (*A) * ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) )
           * exp(1/(r_cap - *a));
  }

  return;
}

/* Calculate two-body term phi_two(r) and its derivative dphi_two(r) */
static void calc_phi_dphi_two(double* A, double* B, double* p, double* q,
                              double* a, double* sigma, double* epsilon,
                              double r, double* phi, double* dphi)
{
  /* Local variables */
  double r_cap;

  r_cap = r/(*sigma);

  if (r_cap >= *a)
  {
    *phi = 0.0;
    *dphi = 0.0;
  }
  else
  {
    *phi = (*epsilon) * (*A) * ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) )
           * exp(1/(r_cap - *a));

    *dphi =
      ( (*q) * pow(r_cap,-((*q)+1)) - (*p * (*B)) * pow(r_cap,-((*p)+1)) )
      - ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) * pow((r_cap - *a),-2);
    *dphi *= (*epsilon / *sigma) * (*A) * exp(1/(r_cap - *a));
   }

   return;
}

/* Calculate two-body term phi_two(r) and its 1st & 2nd derivatives dphi_two(r), d2phi_two(r) */
static void calc_phi_d2phi_two(double* A, double* B, double* p, double* q,
                               double* a, double* sigma, double* epsilon,
                               double r, double* phi, double* dphi,
                               double* d2phi)
{
  /* Local variables */
  double r_cap;

  r_cap = r/(*sigma);

  if (r_cap >= *a)
  {
    *phi = 0.0;
    *dphi = 0.0;
    *d2phi = 0.0;
  }
  else
  {
    *phi = (*epsilon) * (*A) * ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) )
           * exp(1/(r_cap - *a));

    *dphi =  ( (*q) * pow(r_cap,-((*q)+1)) - (*p * *B) * pow(r_cap,-((*p)+1)) )
             - ( (*B) * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) *
             pow((r_cap - *a),-2);
    *dphi *= (*epsilon / *sigma) * (*A) * exp(1/(r_cap - *a));

    *d2phi = ( *B * pow(r_cap,-(*p)) - pow(r_cap,-(*q)) ) *
             ( pow((r_cap - *a),-4) + 2*pow((r_cap - *a),-3) )
             + 2 * ( *p * *B *pow(r_cap,-(*p + 1)) - *q * pow(r_cap,-(*q + 1)) )
             * pow((r_cap - *a),-2) + ( *p * (*p + 1) *  *B * pow(r_cap,-(*p + 2))
             - *q * (*q + 1) * pow(r_cap,-(*q + 2)) );
    *d2phi *= (*epsilon / (*sigma * *sigma)) * (*A) * exp(1/(r_cap - *a));
  }

  return;

}

/* Calculate  three-body term phi_three(r_ij, r_ik, r_jk) */
static void calc_phi_three(double* a, double* lambda, double* gamma,
                           double* sigma, double* epsilon, double* Q,
                           double* costhetat, double rij, double rik,
                           double rjk, double* phi)
{
  /* local variables */
  double c1;

  double rij_cap;
  double rik_cap;

  double costhetajik;

  double diff_costhetajik;

  double exp_ij_ik;

  c1 = *lambda * *epsilon;

  rij_cap = rij/(*sigma);
  rik_cap = rik/(*sigma);

  costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);

  diff_costhetajik = costhetajik - *costhetat;

  exp_ij_ik = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a)));

  if ((rij_cap < *a) && (rik_cap < *a))
  {
    *phi  = c1 * exp_ij_ik * (1 - exp(- *Q *pow(diff_costhetajik,2)));
  }

  return;
}

/* Calculate three-body term phi_three(r_ij, r_ik, r_jk) and its 1st derivative
  dphi_three(r_ij, r_ik, r_jk)

  dphi has three components as derivatives of phi w.r.t. r_ij, r_ik, r_jk
*/
static void calc_phi_dphi_three(double* a, double* lambda, double* gamma,
                                double* sigma, double* epsilon, double* Q,
                                double* costhetat, double rij, double rik,
                                double rjk, double* phi, double* dphi)
{
  /* local variables */
  double c1;
  double c2;

  double rij_cap;
  double rik_cap;

  double costhetajik;

  double diff_costhetajik;

  double costhetajik_ij;
  double costhetajik_ik;
  double costhetajik_jk;

  double exp_ij_ik;

  double d_ij;
  double d_ik;

  double g1, g2, g3, g4;

  c1 = *lambda * *epsilon;
  c2 = c1 / *sigma;

  rij_cap = rij/(*sigma);
  rik_cap = rik/(*sigma);

  costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);

  diff_costhetajik = costhetajik - *costhetat;

  /* Derivatives of cosines w.r.t rij, rik, rjk */
  costhetajik_ij = (*sigma)*(pow(rij,2) - pow(rik,2) +
                   pow(rjk,2))/(2*rij*rij*rik);
  costhetajik_ik = (*sigma)*(pow(rik,2) - pow(rij,2) +
                   pow(rjk,2))/(2*rij*rik*rik);
  costhetajik_jk = -(*sigma)*rjk/(rij*rik);

  /* Variables for simplifying terms */
  exp_ij_ik = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a)));

  d_ij = -(*gamma)*pow((rij_cap - *a),-2);
  d_ik = -(*gamma)*pow((rik_cap - *a),-2);


  if ((rij_cap < *a) && (rik_cap < *a))
  {
    g1  =  pow(diff_costhetajik,2);
    g2  =  exp(- *Q * g1);
    g3  =  1 - g2;
    g4  =  2 * *Q * g2;

    *phi  = c1 * exp_ij_ik * g3;

    /* w.r.t. ij */
    dphi[0] = c2 * exp_ij_ik *
              (g3 * d_ij + g4 * (diff_costhetajik * costhetajik_ij));
    /* w.r.t. ik */
    dphi[1] = c2 * exp_ij_ik *
              (g3 * d_ik + g4 * (diff_costhetajik * costhetajik_ik));
    /* w.r.t. jk */
    dphi[2] = c2 * exp_ij_ik *
              (g4 * (diff_costhetajik * costhetajik_jk));
  }

  return;
}

/* Calculate three-body term phi_three(r_ij, r_ik, r_jk) and its 1st & 2nd derivatives
  dphi_three(r_ij, r_ik, r_jk), d2phi_three(r_ij, r_ik, r_jk)

  dphi has three components as derivatives of phi w.r.t. r_ij, r_ik, r_jk

  d2phi as symmetric Hessian matrix of phi has six components: [0]=(ij,ij), [3]=(ij,ik), [4]=(ij,jk)
                                                                            [1]=(ik,ik), [5]=(ik,jk)
                                                                                         [2]=(jk,jk)
*/
static void calc_phi_d2phi_three(double* a, double* lambda, double* gamma,
                                 double* sigma, double* epsilon, double* Q,
                                 double* costhetat, double rij, double rik,
                                 double rjk, double* phi, double* dphi,
                                 double* d2phi)
{
   /* local variables */
  double c1;
  double c2;
  double c3;

  double rij_cap;
  double rik_cap;

  double costhetajik;

  double diff_costhetajik;

  double costhetajik_ij;
  double costhetajik_ik;
  double costhetajik_jk;

  double costhetajik_ij_ij;
  double costhetajik_ik_ik;
  double costhetajik_jk_jk;
  double costhetajik_ij_ik;
  double costhetajik_ij_jk;
  double costhetajik_ik_jk;

  double exp_ij_ik;

  double d_ij;
  double d_ik;

  double d_ij_2;
  double d_ik_2;

  double dd_ij;
  double dd_ik;

  double g1, g2, g3, g4, g5, powQ;

  double* dg1 = (double *) malloc(3*sizeof(double));

  c1 = *lambda * *epsilon;
  c2 = c1 / *sigma;
  c3 = c2 / *sigma;

  rij_cap = rij/(*sigma);
  rik_cap = rik/(*sigma);

  costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);

  diff_costhetajik = costhetajik - *costhetat;

  /* Derivatives of cosines w.r.t. r_ij, r_ik, r_jk */
  costhetajik_ij = (*sigma)*(pow(rij,2) - pow(rik,2) +
                   pow(rjk,2))/(2*rij*rij*rik);
  costhetajik_ik = (*sigma)*(pow(rik,2) - pow(rij,2) +
                   pow(rjk,2))/(2*rij*rik*rik);
  costhetajik_jk = -(*sigma)*rjk/(rij*rik);

  /* Hessian matrix of cosine */
  costhetajik_ij_ij = (*sigma * *sigma)*(pow(rik,2) -
                      pow(rjk,2))/(rij*rij*rij*rik);
  costhetajik_ik_ik = (*sigma * *sigma)*(pow(rij,2) -
                      pow(rjk,2))/(rij*rik*rik*rik);
  costhetajik_jk_jk = -(*sigma * *sigma)/(rij*rik);
  costhetajik_ij_ik = -(*sigma * *sigma)*(pow(rij,2) +
                      pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik*rik);
  costhetajik_ij_jk = (*sigma * *sigma)*rjk/(rij*rij*rik);
  costhetajik_ik_jk = (*sigma * *sigma)*rjk/(rik*rik*rij);

  /* Variables for simplifying terms */
  exp_ij_ik = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a)));

  d_ij = -(*gamma)*pow((rij_cap - *a),-2);
  d_ik = -(*gamma)*pow((rik_cap - *a),-2);

  d_ij_2 = d_ij * d_ij;
  d_ik_2 = d_ik * d_ik;

  dd_ij = (2 * *gamma)*pow((rij_cap - *a),-3);
  dd_ik = (2 * *gamma)*pow((rik_cap - *a),-3);

  powQ = 4 * *Q * *Q;

  dg1[0]  =  diff_costhetajik * costhetajik_ij;
  dg1[1]  =  diff_costhetajik * costhetajik_ik;
  dg1[2]  =  diff_costhetajik * costhetajik_jk;

  if ((rij_cap < *a) && (rik_cap < *a))
  {
    g1  =  pow(diff_costhetajik,2);
    g2  =  exp(- *Q * g1);
    g3  =  1 - g2;
    g4  =  2 * *Q * g2;
    g5  =  powQ * g2;

    *phi  = c1 * exp_ij_ik * g3;

    /* 1st derivative */
    /* w.r.t. ij */
    dphi[0] = c2 * exp_ij_ik *
              (g3 * d_ij + g4 * (diff_costhetajik * costhetajik_ij));
    /* w.r.t. ik */
    dphi[1] = c2 * exp_ij_ik * (g3 * d_ik + g4 *
              (diff_costhetajik * costhetajik_ik));
    /* w.r.t. jk */
    dphi[2] = c2 * exp_ij_ik * (g4 * (diff_costhetajik * costhetajik_jk));

    /* Hessian */
    d2phi[0] = c3 * exp_ij_ik * ((g3 * (d_ij_2 + dd_ij)) +
               (2 * g4 * d_ij * dg1[0]) - (g5 * pow(dg1[0],2)) +
               (g4 * (diff_costhetajik * costhetajik_ij_ij +
                pow(costhetajik_ij,2))));
    d2phi[1] = c3 * exp_ij_ik * ((g3 * (d_ik_2 + dd_ik)) +
               (2 * g4 * d_ik * dg1[1]) - (g5 * pow(dg1[1],2)) +
               (g4 * (diff_costhetajik * costhetajik_ik_ik +
                pow(costhetajik_ik,2))));
    d2phi[2] = c3 * exp_ij_ik * (- (g5 * pow(dg1[2],2)) +
               (g4 * (diff_costhetajik * costhetajik_jk_jk +
                pow(costhetajik_jk,2))));
    d2phi[3] = c3 * exp_ij_ik * ((g3 * d_ij * d_ik) +
               g4 * (d_ik * dg1[0] + d_ij * dg1[1]) - (g5 * dg1[0] * dg1[1]) +
               (g4 * (diff_costhetajik * costhetajik_ij_ik +
               costhetajik_ij * costhetajik_ik)));
    d2phi[4] = c3 * exp_ij_ik * ((g4 * d_ik * dg1[2]) - (g5 * dg1[1] * dg1[2])
               + (g4 * (diff_costhetajik * costhetajik_ik_jk +
               costhetajik_ik * costhetajik_jk)));
    d2phi[5] = c3 * exp_ij_ik * ((g4 * d_ij * dg1[2]) - (g5 * dg1[0] * dg1[2])
               + (g4 * (diff_costhetajik * costhetajik_ij_jk +
               costhetajik_ij * costhetajik_jk)));
  }

  /* d2phi[0] derivative is w.r.t. rij, rij */
  /* d2phi[1] derivative is w.r.t. rik, rik */
  /* d2phi[2] derivative is w.r.t. rjk, rjk */
  /* d2phi[3] derivative is w.r.t. rij, rik */
  /* d2phi[4] derivative is w.r.t. rij, rjk */
  /* d2phi[5] derivative is w.r.t. rik, rjk */

  free(dg1);

  return;
}

/* Calculate  four-body term phi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl) */
static void calc_phi_four(double* a, double* lambda_2, double* gamma,
                          double* sigma, double* epsilon, double* Q,
                          double* costhetat, double rij, double rik, double ril,
                          double rjk, double rjl, double rkl, double* phi)
{
   /* local variables */
   double c1;

   double rij_cap;
   double rik_cap;
   double ril_cap;

   double costhetajik, costhetajil, costhetakil;

   double diff_costhetajik, diff_costhetajil, diff_costhetakil;

   double exp_ij_ik_il;

   c1 = *lambda_2 * *epsilon;

   rij_cap = rij/(*sigma);
   rik_cap = rik/(*sigma);
   ril_cap = ril/(*sigma);

   costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);
   costhetajil = (pow(rij,2) + pow(ril,2) - pow(rjl,2))/(2*rij*ril);
   costhetakil = (pow(ril,2) + pow(rik,2) - pow(rkl,2))/(2*ril*rik);

   /* Difference of two cosines */
   diff_costhetajik = costhetajik - *costhetat; diff_costhetajil = costhetajil - *costhetat; diff_costhetakil = costhetakil - *costhetat;

   /* Variables for simplifying terms */
   exp_ij_ik_il = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a) + 1/(ril_cap - *a)));

   if ((rij_cap < *a) && (rik_cap < *a) && (ril_cap < *a))
   {
     *phi  = c1 * exp_ij_ik_il * (1 - exp(- *Q *(pow(diff_costhetajik,2) + pow(diff_costhetajil,2) + pow(diff_costhetakil,2))));
   }

   return;
}

/* Calculate four-body term phi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl) and its 1st derivative
  dphi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl)

  dphi has six components as derivatives of phi w.r.t. r_ij, r_ik, r_il, r_jk, r_jl, r_kl
*/
static void calc_phi_dphi_four(double* a, double* lambda_2, double* gamma,
                               double* sigma, double* epsilon, double* Q,
                               double* costhetat, double rij, double rik, double
                               ril, double rjk, double rjl, double rkl, double*
                               phi, double* dphi)
{
  /* local variables */
  double c1;
  double c2;

  double rij_cap;
  double rik_cap;
  double ril_cap;

  double costhetajik, costhetajil, costhetakil;

  double diff_costhetajik, diff_costhetajil, diff_costhetakil;

  /* 1st derivate of cosines w.r.t. ij, ik, il, jk, jl, kl in array order*/
  double* d_costhetajik = (double *) malloc(3*sizeof(double));
  double* d_costhetajil = (double *) malloc(3*sizeof(double));
  double* d_costhetakil = (double *) malloc(3*sizeof(double));

  /* Variables for simplifying terms */
  double exp_ij_ik_il;

  double d_ij;
  double d_ik;
  double d_il;

  double g1, g2, g3, g4;

  c1 = *lambda_2 * *epsilon;
  c2 = c1 / *sigma;

  rij_cap = rij/(*sigma);
  rik_cap = rik/(*sigma);
  ril_cap = ril/(*sigma);

  costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);
  costhetajil = (pow(rij,2) + pow(ril,2) - pow(rjl,2))/(2*rij*ril);
  costhetakil = (pow(ril,2) + pow(rik,2) - pow(rkl,2))/(2*ril*rik);

  /* Difference of two cosines */
  diff_costhetajik = costhetajik - *costhetat; diff_costhetajil
                   = costhetajil - *costhetat; diff_costhetakil
                   = costhetakil - *costhetat;

  /* 1st derivative of cosines */
  d_costhetajik[0] = (*sigma)*(pow(rij,2) - pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik); /* w.r.t. ij */
  d_costhetajik[1] = (*sigma)*(pow(rik,2) - pow(rij,2) + pow(rjk,2))/(2*rij*rik*rik); /* w.r.t. ik */
  d_costhetajik[2] = -(*sigma)*rjk/(rij*rik); /* w.r.t. jk */

  d_costhetajil[0] = (*sigma)*(pow(rij,2) - pow(ril,2) + pow(rjl,2))/(2*rij*rij*ril); /* w.r.t. ij */
  d_costhetajil[1] = (*sigma)*(pow(ril,2) - pow(rij,2) + pow(rjl,2))/(2*rij*ril*ril); /* w.r.t. il */
  d_costhetajil[2] = -(*sigma)*rjl/(rij*ril); /* w.r.t. jl */

  d_costhetakil[0] = (*sigma)*(pow(rik,2) - pow(ril,2) + pow(rkl,2))/(2*rik*rik*ril); /* w.r.t. ik */
  d_costhetakil[1] = (*sigma)*(pow(ril,2) - pow(rik,2) + pow(rkl,2))/(2*rik*ril*ril); /* w.r.t. il */
  d_costhetakil[2] = -(*sigma)*rkl/(rik*ril); /* w.r.t. kl */

  /* Variables for simplifying terms */
  exp_ij_ik_il = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a) + 1/(ril_cap - *a)));

  d_ij = -(*gamma)*pow((rij_cap - *a),-2);
  d_ik = -(*gamma)*pow((rik_cap - *a),-2);
  d_il = -(*gamma)*pow((ril_cap - *a),-2);


  if ((rij_cap < *a) && (rik_cap < *a) && (ril_cap < *a))
  {
    g1  =  pow(diff_costhetajik,2) + pow(diff_costhetajil,2) + pow(diff_costhetakil,2);
    g2  =  exp(- *Q * g1);
    g3  =  1 - g2;
    g4  =  2 * *Q * g2;

    *phi  = c1 * exp_ij_ik_il * g3;

    dphi[0] = c2 * exp_ij_ik_il *  (g3 * d_ij + g4 * (diff_costhetajik * d_costhetajik[0]
                                                + diff_costhetajil * d_costhetajil[0])); /* w.r.t. ij */
    dphi[1] = c2 * exp_ij_ik_il *  (g3 * d_ik + g4 * (diff_costhetajik * d_costhetajik[1]
                                                + diff_costhetakil * d_costhetakil[0])); /* w.r.t. ik */
    dphi[2] = c2 * exp_ij_ik_il *  (g3 * d_il + g4 * (diff_costhetajil * d_costhetajil[1]
                                                + diff_costhetakil * d_costhetakil[1])); /* w.r.t. il */
    dphi[3] = c2 * exp_ij_ik_il *  (g4 * (diff_costhetajik * d_costhetajik[2]));         /* w.r.t. jk */
    dphi[4] = c2 * exp_ij_ik_il *  (g4 * (diff_costhetajil * d_costhetajil[2]));         /* w.r.t. jl */
    dphi[5] = c2 * exp_ij_ik_il *  (g4 * (diff_costhetakil * d_costhetakil[2]));         /* w.r.t. kl */
  }

  free(d_costhetajik);
  free(d_costhetajil);
  free(d_costhetakil);

  return;
}

/* Calculate four-body term phi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl) and its 1st & 2nd derivatives
  dphi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl), d2phi_four(r_ij, r_ik, r_il, r_jk, r_jl, r_kl)

  dphi has six components as derivatives of phi w.r.t. r_ij, r_ik, r_il, r_jk, r_jl, r_kl

  d2phi as symmetric Hessian matrix of phi has 21 components:
  [0]=(ij,ij), [6]=(ij,ik), [11]=(ij,il), [15]=(ij,jk), [18]=(ij,jl), [20]=(ij,kl)
               [1]=(ik,ik), [7] =(ik,il), [12]=(ik,jk), [16]=(ik,jl), [19]=(ik,kl)
                            [2] =(il,il), [8] =(il,jk), [13]=(il,jl), [17]=(il,kl)
                                          [3] =(jk,jk), [9] =(jk,jl), [14]=(jk,kl)
                                                        [4] =(jl,jl), [10]=(jl,kl)
                                                                       [5] =(kl,kl)
*/
static void calc_phi_d2phi_four(double* a, double* lambda_2, double* gamma,
                                double* sigma, double* epsilon, double* Q,
                                double* costhetat, double rij, double rik,
                                double ril, double rjk, double rjl, double rkl,
                                double* phi, double* dphi, double* d2phi)
{
    /* local variables */
   double c1;
   double c2;
   double c3;

   double rij_cap;
   double rik_cap;
   double ril_cap;

   double costhetajik, costhetajil, costhetakil;

   double diff_costhetajik, diff_costhetajil, diff_costhetakil;

   /* 1st derivate of cosines w.r.t. ij, ik, il, jk, jl, kl in array order*/
   double* d_costhetajik = (double *) malloc(3*sizeof(double));
   double* d_costhetajil = (double *) malloc(3*sizeof(double));
   double* d_costhetakil = (double *) malloc(3*sizeof(double));

   /* Hessian of cosines in the given array order*/
   double* d2_costhetajik = (double *) malloc(6*sizeof(double));
   /* [0] = (ij,ij), [1] = (ik,ik), [2] = (jk,jk), [3] = (ij,ik),
      [4] = (ij,jk), [5] = (ik,jk) */
   double* d2_costhetajil = (double *) malloc(6*sizeof(double));
   double* d2_costhetakil = (double *) malloc(6*sizeof(double));

   double exp_ij_ik_il;

   double d_ij;
   double d_ik;
   double d_il;

   double d_ij_2;
   double d_ik_2;
   double d_il_2;

   double dd_ij;
   double dd_ik;
   double dd_il;

   double g1, g2, g3, g4, g5, powQ;

   double* dg1 = (double *) malloc(6*sizeof(double));
   double* dg2 = (double *) malloc(6*sizeof(double));

   int i;

   c1 = *lambda_2 * *epsilon;
   c2 = c1 / *sigma;
   c3 = c2 / *sigma;

   rij_cap = rij/(*sigma);
   rik_cap = rik/(*sigma);
   ril_cap = ril/(*sigma);

   costhetajik = (pow(rij,2) + pow(rik,2) - pow(rjk,2))/(2*rij*rik);
   costhetajil = (pow(rij,2) + pow(ril,2) - pow(rjl,2))/(2*rij*ril);
   costhetakil = (pow(ril,2) + pow(rik,2) - pow(rkl,2))/(2*ril*rik);

   /* Difference of two cosines */
   diff_costhetajik = costhetajik - *costhetat; diff_costhetajil
                    = costhetajil - *costhetat; diff_costhetakil
                    = costhetakil - *costhetat;

   /* 1st derivative of cosines */
   /* w.r.t. ij */
   d_costhetajik[0] = (*sigma)*(pow(rij,2) - pow(rik,2) + pow(rjk,2))/(2*rij*rij*rik);
   d_costhetajik[1] = (*sigma)*(pow(rik,2) - pow(rij,2) + pow(rjk,2))/(2*rij*rik*rik); /* w.r.t. ik */
   d_costhetajik[2] = -(*sigma)*rjk/(rij*rik); /* w.r.t. jk */

   d_costhetajil[0] = (*sigma)*(pow(rij,2) - pow(ril,2) + pow(rjl,2))/(2*rij*rij*ril); /* w.r.t. ij */
   d_costhetajil[1] = (*sigma)*(pow(ril,2) - pow(rij,2) + pow(rjl,2))/(2*rij*ril*ril); /* w.r.t. il */
   d_costhetajil[2] = -(*sigma)*rjl/(rij*ril); /* w.r.t. jl */

   d_costhetakil[0] = (*sigma)*(pow(rik,2) - pow(ril,2) + pow(rkl,2))/(2*rik*rik*ril); /* w.r.t. ik */
   d_costhetakil[1] = (*sigma)*(pow(ril,2) - pow(rik,2) + pow(rkl,2))/(2*rik*ril*ril); /* w.r.t. il */
   d_costhetakil[2] = -(*sigma)*rkl/(rik*ril); /* w.r.t. kl */

   /* Hessian matrix of cosines */
   d2_costhetajik[0] = (*sigma * *sigma)*(pow(rik,2) - pow(rjk,2))/(pow(rij,3)*rik);                         /* w.r.t. ij,ij */
   d2_costhetajik[1] = -(*sigma * *sigma)*(pow(rij,2) + pow(rik,2) + pow(rjk,2))/(2*pow(rij,2)*pow(rik,2));  /* w.r.t. ij,ik */
   d2_costhetajik[2] = (*sigma * *sigma)*rjk/(pow(rij,2)*rik);                                               /* w.r.t. ij,jk */
   d2_costhetajik[3] = (*sigma * *sigma)*(pow(rij,2) - pow(rjk,2))/(rij*pow(rik,3));                         /* w.r.t. ik,ik */
   d2_costhetajik[4] = (*sigma * *sigma)*rjk/(pow(rik,2)*rij);                                               /* w.r.t. ik,jk */
   d2_costhetajik[5] = -(*sigma * *sigma)/(rij*rik);                                                         /* w.r.t. jk,jk */

   d2_costhetajil[0] = (*sigma * *sigma)*(pow(ril,2) - pow(rjl,2))/(pow(rij,3)*ril);                         /* w.r.t. ij,ij */
   d2_costhetajil[1] = -(*sigma * *sigma)*(pow(rij,2) + pow(ril,2) + pow(rjl,2))/(2*pow(rij,2)*pow(ril,2));  /* w.r.t. ij,il */
   d2_costhetajil[2] = (*sigma * *sigma)*rjl/(pow(rij,2)*ril);                                               /* w.r.t. ij,jl */
   d2_costhetajil[3] = (*sigma * *sigma)*(pow(rij,2) - pow(rjl,2))/(rij*pow(ril,3));                         /* w.r.t. il,il */
   d2_costhetajil[4] = (*sigma * *sigma)*rjl/(pow(ril,2)*rij);                                               /* w.r.t. il,jl */
   d2_costhetajil[5] = -(*sigma * *sigma)/(rij*ril);                                                         /* w.r.t. jl,jl */

   d2_costhetakil[0] = (*sigma * *sigma)*(pow(ril,2) - pow(rkl,2))/(pow(rik,3)*ril);                         /* w.r.t. ik,ik */
   d2_costhetakil[1] = -(*sigma * *sigma)*(pow(rik,2) + pow(ril,2) + pow(rkl,2))/(2*pow(rik,2)*pow(ril,2));  /* w.r.t. ik,il */
   d2_costhetakil[2] = (*sigma * *sigma)*rkl/(pow(rik,2)*ril);                                               /* w.r.t. ik,kl */
   d2_costhetakil[3] = (*sigma * *sigma)*(pow(rik,2) - pow(rkl,2))/(rik*pow(ril,3));                         /* w.r.t. il,il */
   d2_costhetakil[4] = (*sigma * *sigma)*rkl/(pow(ril,2)*rik);                                               /* w.r.t. il,kl */
   d2_costhetakil[5] = -(*sigma * *sigma)/(rik*ril);                                                         /* w.r.t. kl,kl */

   /* Variables for simplifying terms */
   exp_ij_ik_il = exp((*gamma) * (1/(rij_cap - *a) + 1/(rik_cap - *a) + 1/(ril_cap - *a)));

   d_ij = -(*gamma)*pow((rij_cap - *a),-2);
   d_ik = -(*gamma)*pow((rik_cap - *a),-2);
   d_il = -(*gamma)*pow((ril_cap - *a),-2);

   d_ij_2 = d_ij * d_ij;
   d_ik_2 = d_ik * d_ik;
   d_il_2 = d_il * d_il;

   dd_ij = (2 * *gamma)*pow((rij_cap - *a),-3);
   dd_ik = (2 * *gamma)*pow((rik_cap - *a),-3);
   dd_il = (2 * *gamma)*pow((ril_cap - *a),-3);

   *phi = 0.0;

   for(i = 0; i < 6; i++) dphi[i] = 0.0;
   for(i = 0; i < 21; i++) d2phi[i] = 0.0;

   powQ = 4 * *Q * *Q;

   if ((rij_cap < *a) && (rik_cap < *a) && (ril_cap < *a))
   {
     g1  =  pow(diff_costhetajik,2) + pow(diff_costhetajil,2) + pow(diff_costhetakil,2);
     g2  =  exp(- *Q * g1);
     g3  =  1 - g2;
     g4  =  2 * *Q * g2;
     g5  =  powQ * g2;

     dg1[0]  =  diff_costhetajik * d_costhetajik[0] + diff_costhetajil * d_costhetajil[0];
     dg1[1]  =  diff_costhetajik * d_costhetajik[1] + diff_costhetakil * d_costhetakil[0];
     dg1[2]  =  diff_costhetajil * d_costhetajil[1] + diff_costhetakil * d_costhetakil[1];
     dg1[3]  =  diff_costhetajik * d_costhetajik[2];
     dg1[4]  =  diff_costhetajil * d_costhetajil[2];
     dg1[5]  =  diff_costhetakil * d_costhetakil[2];

     *phi  += exp_ij_ik_il * g3;

     /* 1st derivative */
     dphi[0] += exp_ij_ik_il *  (g3 * d_ij + g4 * dg1[0]); /* w.r.t. ij */
     dphi[1] += exp_ij_ik_il *  (g3 * d_ik + g4 * dg1[1]); /* w.r.t. ik */
     dphi[2] += exp_ij_ik_il *  (g3 * d_il + g4 * dg1[2]); /* w.r.t. il */
     dphi[3] += exp_ij_ik_il *  g4 * dg1[3];               /* w.r.t. jk */
     dphi[4] += exp_ij_ik_il *  g4 * dg1[4];               /* w.r.t. jl */
     dphi[5] += exp_ij_ik_il *  g4 * dg1[5];               /* w.r.t. kl */

     /* Hessian */
     dg2[0]  =  diff_costhetajik * d2_costhetajik[0] + diff_costhetajil * d2_costhetajil[0] + pow(d_costhetajik[0],2) + pow(d_costhetajil[0],2);
     dg2[1]  =  diff_costhetajik * d2_costhetajik[1] + d_costhetajik[0] * d_costhetajik[1];
     dg2[2]  =  diff_costhetajil * d2_costhetajil[1] + d_costhetajil[0] * d_costhetajil[1];
     dg2[3]  =  diff_costhetajik * d2_costhetajik[2] + d_costhetajik[0] * d_costhetajik[2];
     dg2[4]  =  diff_costhetajil * d2_costhetajil[2] + d_costhetajil[0] * d_costhetajil[2];

     d2phi[0] += exp_ij_ik_il * ((g3 * (d_ij_2 + dd_ij)) + (2 * g4 * d_ij * dg1[0]) - (g5 * pow(dg1[0],2)) + (g4 * dg2[0]));          /* w.r.t. (ij,ij) */
     d2phi[6] += exp_ij_ik_il * ((g3 * d_ij * d_ik) + g4 * (d_ik * dg1[0] + d_ij * dg1[1]) - (g5 * dg1[0] * dg1[1]) + (g4 * dg2[1])); /* w.r.t. (ij,ik) */
     d2phi[11] += exp_ij_ik_il * ((g3 * d_ij * d_il) + g4 * (d_il * dg1[0] + d_ij * dg1[2]) - (g5 * dg1[0] * dg1[2]) + (g4 * dg2[2])); /* w.r.t. (ij,il) */
     d2phi[15] += exp_ij_ik_il * ((g4 * d_ij * dg1[3]) - (g5 * dg1[0] * dg1[3]) + (g4 * dg2[3]));                                      /* w.r.t. (ij,jk) */
     d2phi[18] += exp_ij_ik_il * ((g4 * d_ij * dg1[4]) - (g5 * dg1[0] * dg1[4]) + (g4 * dg2[4]));                                      /* w.r.t. (ij,jl) */
     d2phi[20] += exp_ij_ik_il * ((g4 * d_ij * dg1[5]) - (g5 * dg1[0] * dg1[5]));                                                      /* w.r.t. (ij,kl) */

     dg2[0]  =  diff_costhetajik * d2_costhetajik[3] + diff_costhetakil * d2_costhetakil[0] + pow(d_costhetajik[1],2) + pow(d_costhetakil[0],2);
     dg2[1]  =  diff_costhetakil * d2_costhetakil[1] + d_costhetakil[0] * d_costhetakil[1];
     dg2[2]  =  diff_costhetajik * d2_costhetajik[4] + d_costhetajik[1] * d_costhetajik[2];
     dg2[3]  =  diff_costhetakil * d2_costhetakil[2] + d_costhetakil[0] * d_costhetakil[2];

     d2phi[1] += exp_ij_ik_il * ((g3 * (d_ik_2 + dd_ik)) + (2 * g4 * d_ik * dg1[1]) - (g5 * pow(dg1[1],2)) + (g4 * dg2[0]));          /* w.r.t. (ik,ik) */
     d2phi[7] += exp_ij_ik_il * ((g3 * d_ik * d_il) + g4 * (d_il * dg1[1] + d_ik * dg1[2]) - (g5 * dg1[1] * dg1[2]) + (g4 * dg2[1])); /* w.r.t. (ik,il) */
     d2phi[12] += exp_ij_ik_il * ((g4 * d_ik * dg1[3]) - (g5 * dg1[1] * dg1[3]) + (g4 * dg2[2]));                                      /* w.r.t. (ik,jk) */
     d2phi[16] += exp_ij_ik_il * ((g4 * d_ik * dg1[4]) - (g5 * dg1[1] * dg1[4]));                                                      /* w.r.t. (ik,jl) */
     d2phi[19] += exp_ij_ik_il * ((g4 * d_ik * dg1[5]) - (g5 * dg1[1] * dg1[5]) + (g4 * dg2[3]));                                     /* w.r.t. (ik,kl) */

     dg2[0]  =  diff_costhetajil * d2_costhetajil[3] + diff_costhetakil * d2_costhetakil[3] + pow(d_costhetajil[1],2) + pow(d_costhetakil[1],2);
     dg2[1]  =  diff_costhetajil * d2_costhetajil[4] + d_costhetajil[1] * d_costhetajil[2];
     dg2[2]  =  diff_costhetakil * d2_costhetakil[4] + d_costhetakil[1] * d_costhetakil[2];

     d2phi[2] += exp_ij_ik_il * ((g3 * (d_il_2 + dd_il)) + (2 * g4 * d_il * dg1[2]) - (g5 * pow(dg1[2],2)) + (g4 * dg2[0]));         /* w.r.t. (il,il) */
     d2phi[8] += exp_ij_ik_il * ((g4 * d_il * dg1[3]) - (g5 * dg1[2] * dg1[3]));                                                     /* w.r.t. (il,jk) */
     d2phi[13] += exp_ij_ik_il * ((g4 * d_il * dg1[4]) - (g5 * dg1[2] * dg1[4]) + (g4 * dg2[1]));                                     /* w.r.t. (il,jl) */
     d2phi[17] += exp_ij_ik_il * ((g4 * d_il * dg1[5]) - (g5 * dg1[2] * dg1[5]) + (g4 * dg2[2]));                                     /* w.r.t. (il,kl) */

     dg2[0]  =  diff_costhetajik * d2_costhetajik[5] +  pow(d_costhetajik[2],2);

     d2phi[3] += exp_ij_ik_il * (- (g5 * pow(dg1[3],2)) + (g4 * dg2[0]));                                                            /* w.r.t. (jk,jk) */
     d2phi[9] += exp_ij_ik_il * (- g5 * dg1[3] * dg1[4]);                                                                            /* w.r.t. (jk,jl) */
     d2phi[14] += exp_ij_ik_il * (- g5 * dg1[3] * dg1[5]);                                                                            /* w.r.t. (jk,kl) */

     dg2[0]  =  diff_costhetajil * d2_costhetajil[5] +  pow(d_costhetajil[2],2);

     d2phi[4] += exp_ij_ik_il * (- (g5 * pow(dg1[4],2)) + (g4 * dg2[0]));                                                            /* w.r.t. (jl,jl) */
     d2phi[10] += exp_ij_ik_il * (- g5 * dg1[4] * dg1[5]);                                                                            /* w.r.t. (jl,kl) */

     dg2[0]  =  diff_costhetakil * d2_costhetakil[5] +  pow(d_costhetakil[2],2);

     d2phi[5] += exp_ij_ik_il * (- (g5 * pow(dg1[5],2)) + (g4 * dg2[0]));                                                            /* w.r.t. (kl,kl) */
   }

   *phi  *= c1;

   for(i = 0; i < 6; i++) dphi[i] *= c2;
   for(i = 0; i < 21; i++) d2phi[i] *= c3;

    free(d_costhetajik);
    free(d_costhetajil);
    free(d_costhetakil);

    free(d2_costhetajik);
    free(d2_costhetajil);
    free(d2_costhetakil);

    free(dg1);
    free(dg2);

   return;
}

/* compute function */
#include "KIM_ModelComputeLogMacros.h"
static int compute(KIM_ModelCompute const * const modelCompute,
                   KIM_ModelComputeArguments const * const modelComputeArguments)
{
  /* local variables */
  double R1, R2, R3, R4, R5, R6;
  double R_pairs[2];
  double *pR_pairs = &(R_pairs[0]);
  double Rsqij, Rsqik, Rsqil, Rsqjk, Rsqjl, Rsqkl;
  double phi_two;
  double dphi_two;
  double d2phi_two;
  double dEidr_two;
  double d2Eidr_two;
  double phi_three;
  double* dphi_three;
  double* d2phi_three;
  double* dEidr_three;
  double* d2Eidr_three;
  double phi_four;
  double* dphi_four;
  double* d2phi_four;
  double* dEidr_four;
  double* d2Eidr_four;
  double Rij[DIM];
  double Rik[DIM];
  double Ril[DIM];
  double Rjk[DIM];
  double Rjl[DIM];
  double Rkl[DIM];
  double *pRij = &(Rij[0]);
  double *pRik = &(Rik[0]);
  double *pRil = &(Ril[0]);
  double *pRjk = &(Rjk[0]);
  double *pRjl = &(Rjl[0]);
  double *pRkl = &(Rkl[0]);
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
  int kk;
  int l;
  int ll;
  int kdim;
  int const * neighListOfCurrentAtom;
  struct model_buffer* buffer;
  int comp_energy;
  int comp_force;
  int comp_particleEnergy;
  int comp_process_dEdr;
  int comp_process_d2Edr2;

  int* nAtoms;
  int* particleSpeciesCodes;
  int* particleContributing;
  double* cutsq;
  double* A;
  double* B;
  double* p;
  double* q;
  double* a;
  double* lambda;
  double* lambda_2;
  double* gamma;
  double* sigma;
  double* epsilon;
  double* Q;
  double* costhetat;
  double* coords;
  double* energy;
  double* force;
  double* particleEnergy;
  int numOfAtomNeigh;
  int iSpecies;
  int jSpecies;
  int kSpecies;
  int lSpecies;
  int interaction_index;

  dphi_three = (double *) malloc(3*sizeof(double));
  d2phi_three = (double *) malloc(6*sizeof(double));
  dEidr_three = (double *) malloc(3*sizeof(double));
  d2Eidr_three = (double *) malloc(6*sizeof(double));

  dphi_four = (double *) malloc(6*sizeof(double));
  d2phi_four = (double *) malloc(21*sizeof(double));
  dEidr_four = (double *) malloc(6*sizeof(double));
  d2Eidr_four = (double *) malloc(21*sizeof(double));

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

  /* get buffer from KIM object */
  KIM_ModelCompute_GetModelBufferPointer(modelCompute, (void **) &buffer);

  /* unpack the Model's parameters stored in the KIM API object */
  if(buffer->num_interactions == 1)
  {
       interaction_index = 0;
       A         = &(buffer->A)[interaction_index];
       B         = &(buffer->B)[interaction_index];
       p         = &(buffer->p)[interaction_index];
       q         = &(buffer->q)[interaction_index];
       lambda    = &(buffer->lambda)[interaction_index];
       lambda_2  = &(buffer->lambda_2)[interaction_index];
       gamma     = &(buffer->gamma)[interaction_index];
       sigma     = &(buffer->sigma)[interaction_index];
       epsilon   = &(buffer->epsilon)[interaction_index];
       Q         = &(buffer->Q)[interaction_index];
       costhetat = &(buffer->costhetat)[interaction_index];
       a         = &(buffer->a)[interaction_index];
       cutsq     = &(buffer->cutsq)[interaction_index];
  }

  /* Check to be sure that the atom types are correct */
  /**/
  for (i = 0; i < *nAtoms; ++i)
  {
    if ( particleSpeciesCodes[i] > 2)
    {
      ier = TRUE;
      LOG_ERROR("Unexpected species type detected");
      return ier;
    }
  }

  /* initialize potential energies, forces, and virial term */
  if (comp_particleEnergy)
  {
     for (i = 0; i < *nAtoms; ++i)
     {
        particleEnergy[i] = 0.0;
     }
  }
  else if (comp_energy)
  {
     *energy = 0.0;
  }

  if (comp_force)
  {
     for (i = 0; i < *nAtoms; ++i)
     {
        for (kdim = 0; kdim < DIM; ++kdim)
        {
           force[i*DIM + kdim] = 0.0;
        }
     }
  }

  /* Compute enery and forces */

  /* loop over particles and compute energy and forces */
  for(i=0; i < *nAtoms; ++i)
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
        j = neighListOfCurrentAtom[jj]; /* get neighbor ID */

        jSpecies = particleSpeciesCodes[j];

        /* get corresponding parameters if there are two atomic species */
        if(buffer->num_interactions == 16)
        {
          if (iSpecies == SPEC1 && jSpecies == SPEC1) interaction_index = 0;
          else if (iSpecies == SPEC2 && jSpecies == SPEC2) interaction_index = 15;
          else interaction_index = 7;

          A         = &(buffer->A)[interaction_index];
          B         = &(buffer->B)[interaction_index];
          p         = &(buffer->p)[interaction_index];
          q         = &(buffer->q)[interaction_index];
          lambda    = &(buffer->lambda)[interaction_index];
          lambda_2  = &(buffer->lambda_2)[interaction_index];
          gamma     = &(buffer->gamma)[interaction_index];
          sigma     = &(buffer->sigma)[interaction_index];
          epsilon   = &(buffer->epsilon)[interaction_index];
          Q         = &(buffer->Q)[interaction_index];
          costhetat = &(buffer->costhetat)[interaction_index];
          a         = &(buffer->a)[interaction_index];
          cutsq     = &(buffer->cutsq)[interaction_index];
        }

        /* compute relative position vector and squared distance */
        Rsqij = 0.0;
        for (kdim = 0; kdim < DIM; ++kdim)
        {
          Rij[kdim] = coords[j*DIM + kdim] - coords[i*DIM + kdim];

          /* compute squared distance */
          Rsqij += Rij[kdim]*Rij[kdim];
        }

        /* compute energy and force */
        if (Rsqij > *cutsq) continue; /* particles are not interacting  */
        R1 = sqrt(Rsqij);
        if (comp_process_d2Edr2)
        {
          /* compute pair potential and its derivatives */
          calc_phi_d2phi_two(A, B, p, q, a, sigma, epsilon,
                             R1, &phi_two, &dphi_two, &d2phi_two);

          /* compute dEidr */
          dEidr_two  = 0.5*dphi_two;
          d2Eidr_two = 0.5*d2phi_two;
        }
        else if (comp_force || comp_process_dEdr)
        {
          /* compute pair potential and its derivative */
          calc_phi_dphi_two(A, B, p, q, a, sigma, epsilon,
                            R1, &phi_two, &dphi_two);

          /* compute dEidr */
          dEidr_two = 0.5*dphi_two;
        }
        else
        {
          /* compute just pair potential */
          calc_phi_two(A, B, p, q, a, sigma, epsilon,
                       R1, &phi_two);
        }

        /* contribution to energy */
        if (comp_particleEnergy)
        {
          particleEnergy[i] += 0.5*phi_two;
        }
        if (comp_energy)
        {
          *energy += 0.5*phi_two;
        }

        /* contribution to process_dEdr */
        if (comp_process_dEdr)
        {
          ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                  modelComputeArguments,
                  dEidr_two, R1, pRij, i, j);
        }

        /* contribution to process_d2Edr2 */
        if (comp_process_d2Edr2)
        {
          R_pairs[0] = R_pairs[1] = R1;
          Rij_pairs[0][0] = Rij_pairs[1][0] = Rij[0];
          Rij_pairs[0][1] = Rij_pairs[1][1] = Rij[1];
          Rij_pairs[0][2] = Rij_pairs[1][2] = Rij[2];
          i_pairs[0] = i_pairs[1] = i;
          j_pairs[0] = j_pairs[1] = j;

          ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                  modelComputeArguments,
                  d2Eidr_two, pR_pairs, pRij_pairs, pi_pairs, pj_pairs);
        }

        /* contribution to forces */
        if (comp_force)
        {
          for (kdim = 0; kdim < DIM; ++kdim)
          {
            force[i*DIM + kdim] += dEidr_two*Rij[kdim]/R1; /* accumulate force on atom i */
            force[j*DIM + kdim] -= dEidr_two*Rij[kdim]/R1; /* accumulate force on atom j */
          }
        }

        /* Start adding three body terms */
        /*********************************/
        if(jj == numOfAtomNeigh-1) continue;

        for (kk = jj+1; kk < numOfAtomNeigh; ++kk)
        {
          k = neighListOfCurrentAtom[kk]; /* get neighbor ID */
          kSpecies = particleSpeciesCodes[k];

          if(buffer->num_interactions == 16)
          {
            if (iSpecies == SPEC1 && jSpecies == SPEC1 && kSpecies == SPEC1)
              interaction_index = 0;
            else if (iSpecies == SPEC1 && jSpecies == SPEC1 && kSpecies == SPEC2)
              interaction_index = 2;
            else if (iSpecies == SPEC1 && jSpecies == SPEC2 && kSpecies == SPEC1)
              interaction_index = 4;
            else if (iSpecies == SPEC1 && jSpecies == SPEC2 && kSpecies == SPEC2)
              interaction_index = 6;
            else if (iSpecies == SPEC2 && jSpecies == SPEC1 && kSpecies == SPEC1)
              interaction_index = 8;
            else if (iSpecies == SPEC2 && jSpecies == SPEC1 && kSpecies == SPEC2)
              interaction_index = 10;
            else if (iSpecies == SPEC2 && jSpecies == SPEC2 && kSpecies == SPEC1)
               interaction_index = 12;
            else
               interaction_index = 14;

            A         = &(buffer->A)[interaction_index];
            B         = &(buffer->B)[interaction_index];
            p         = &(buffer->p)[interaction_index];
            q         = &(buffer->q)[interaction_index];
            lambda    = &(buffer->lambda)[interaction_index];
            lambda_2  = &(buffer->lambda_2)[interaction_index];
            gamma     = &(buffer->gamma)[interaction_index];
            sigma     = &(buffer->sigma)[interaction_index];
            epsilon   = &(buffer->epsilon)[interaction_index];
            Q         = &(buffer->Q)[interaction_index];
            costhetat = &(buffer->costhetat)[interaction_index];
            a         = &(buffer->a)[interaction_index];
            cutsq     = &(buffer->cutsq)[interaction_index];
          }

          /* compute relative position vector and squared distance */
          Rsqik = 0.0;
          Rsqjk = 0.0;
          for (kdim = 0; kdim < DIM; ++kdim)
          {
            Rik[kdim] = coords[k*DIM + kdim] - coords[i*DIM + kdim];
            Rjk[kdim] = Rik[kdim] - Rij[kdim];

            /* compute squared distance */
            Rsqik += Rik[kdim]*Rik[kdim];
            Rsqjk += Rjk[kdim]*Rjk[kdim];
          }

          /* compute energy and force */
          if (Rsqik > *cutsq) continue; /* particles are interacting ? */

          R2 = sqrt(Rsqik);
          R4 = sqrt(Rsqjk);

          if (comp_process_d2Edr2)
          {
            /* compute three-body potential and its derivatives */
            calc_phi_d2phi_three(a, lambda, gamma, sigma, epsilon, Q, costhetat,
                                 R1, R2, R4, &phi_three, dphi_three, d2phi_three);

            /* compute dEidr */
            dEidr_three[0] = dphi_three[0];
            dEidr_three[1] = dphi_three[1];
            dEidr_three[2] = dphi_three[2];

            d2Eidr_three[0] = d2phi_three[0];
            d2Eidr_three[1] = d2phi_three[1];
            d2Eidr_three[2] = d2phi_three[2];
            d2Eidr_three[3] = d2phi_three[3];
            d2Eidr_three[4] = d2phi_three[4];
            d2Eidr_three[5] = d2phi_three[5];
          }
          else if (comp_force || comp_process_dEdr)
          {
            /* compute three-body potential and its derivative */
            calc_phi_dphi_three(a, lambda, gamma, sigma, epsilon, Q, costhetat,
                                   R1, R2, R4, &phi_three, dphi_three);

            /* compute dEidr */
            dEidr_three[0] = dphi_three[0];
            dEidr_three[1] = dphi_three[1];
            dEidr_three[2] = dphi_three[2];
          }
          else
          {
            /* compute just three-body potential */
            calc_phi_three(a, lambda, gamma, sigma, epsilon, Q, costhetat,
                           R1, R2, R4, &phi_three);
          }

          /* contribution to energy */
          if (comp_particleEnergy)
          {
            particleEnergy[i] += phi_three;
          }
          if (comp_energy)
          {
            *energy += phi_three;
          }

          /* contribution to process_dEdr */
          if (comp_process_dEdr)
          {
            ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                    modelComputeArguments,
                    dEidr_three[0], R1, pRij, i, j);
            ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                    modelComputeArguments,
                    dEidr_three[1], R2, pRik, i, k);
            ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                    modelComputeArguments,
                    dEidr_three[2], R4, pRjk, i, k);
          }

          /* contribution to process_d2Edr2 */
          if (comp_process_d2Edr2)
          {
            R_pairs[0] = R_pairs[1] = R1;
            Rij_pairs[0][0] = Rij_pairs[1][0] = Rij[0];
            Rij_pairs[0][1] = Rij_pairs[1][1] = Rij[1];
            Rij_pairs[0][2] = Rij_pairs[1][2] = Rij[2];
            i_pairs[0] = i_pairs[1] = i;
            j_pairs[0] = j_pairs[1] = j;

            ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                    modelComputeArguments,
                    d2Eidr_three[0], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

            R_pairs[0] = R_pairs[1] = R2;
            Rij_pairs[0][0] = Rij_pairs[1][0] = Rik[0];
            Rij_pairs[0][1] = Rij_pairs[1][1] = Rik[1];
            Rij_pairs[0][2] = Rij_pairs[1][2] = Rik[2];
            i_pairs[0] = i_pairs[1] = i;
            j_pairs[0] = j_pairs[1] = k;

            ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                    modelComputeArguments,
                    d2Eidr_three[1], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

            R_pairs[0] = R_pairs[1] = R4;
            Rij_pairs[0][0] = Rij_pairs[1][0] = Rjk[0];
            Rij_pairs[0][1] = Rij_pairs[1][1] = Rjk[1];
            Rij_pairs[0][2] = Rij_pairs[1][2] = Rjk[2];
            i_pairs[0] = i_pairs[1] = j;
            j_pairs[0] = j_pairs[1] = k;

            ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                    modelComputeArguments,
                    d2Eidr_three[2], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

            R_pairs[0] = R1;
            R_pairs[1] = R2;
            Rij_pairs[0][0] = Rij[0];
            Rij_pairs[0][1] = Rij[1];
            Rij_pairs[0][2] = Rij[2];
            Rij_pairs[1][0] = Rik[0];
            Rij_pairs[1][1] = Rik[1];
            Rij_pairs[1][2] = Rik[2];
            i_pairs[0] = i;
            j_pairs[0] = j;
            i_pairs[1] = i;
            j_pairs[1] = k;

            ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                    modelComputeArguments,
                    d2Eidr_three[3], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

            R_pairs[0] = R2;
            R_pairs[1] = R1;
            Rij_pairs[0][0] = Rik[0];
            Rij_pairs[0][1] = Rik[1];
            Rij_pairs[0][2] = Rik[2];
            Rij_pairs[1][0] = Rij[0];
            Rij_pairs[1][1] = Rij[1];
            Rij_pairs[1][2] = Rij[2];
            i_pairs[0] = i;
            j_pairs[0] = k;
            i_pairs[1] = i;
            j_pairs[1] = j;

            ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                    modelComputeArguments,
                    d2Eidr_three[3], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

            R_pairs[0] = R1;
            R_pairs[1] = R4;
            Rij_pairs[0][0] = Rij[0];
            Rij_pairs[0][1] = Rij[1];
            Rij_pairs[0][2] = Rij[2];
            Rij_pairs[1][0] = Rjk[0];
            Rij_pairs[1][1] = Rjk[1];
            Rij_pairs[1][2] = Rjk[2];
            i_pairs[0] = i;
            j_pairs[0] = j;
            i_pairs[1] = j;
            j_pairs[1] = k;

            ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                    modelComputeArguments,
                    d2Eidr_three[4], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

            R_pairs[0] = R4;
            R_pairs[1] = R1;
            Rij_pairs[0][0] = Rjk[0];
            Rij_pairs[0][1] = Rjk[1];
            Rij_pairs[0][2] = Rjk[2];
            Rij_pairs[1][0] = Rij[0];
            Rij_pairs[1][1] = Rij[1];
            Rij_pairs[1][2] = Rij[2];
            i_pairs[0] = j;
            j_pairs[0] = k;
            i_pairs[1] = i;
            j_pairs[1] = j;

            ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                    modelComputeArguments,
                    d2Eidr_three[4], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

            R_pairs[0] = R2;
            R_pairs[1] = R4;
            Rij_pairs[0][0] = Rik[0];
            Rij_pairs[0][1] = Rik[1];
            Rij_pairs[0][2] = Rik[2];
            Rij_pairs[1][0] = Rjk[0];
            Rij_pairs[1][1] = Rjk[1];
            Rij_pairs[1][2] = Rjk[2];
            i_pairs[0] = i;
            j_pairs[0] = k;
            i_pairs[1] = j;
            j_pairs[1] = k;

            ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                    modelComputeArguments,
                    d2Eidr_three[5], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

            R_pairs[0] = R4;
            R_pairs[1] = R2;
            Rij_pairs[0][0] = Rjk[0];
            Rij_pairs[0][1] = Rjk[1];
            Rij_pairs[0][2] = Rjk[2];
            Rij_pairs[1][0] = Rik[0];
            Rij_pairs[1][1] = Rik[1];
            Rij_pairs[1][2] = Rik[2];
            i_pairs[0] = j;
            j_pairs[0] = k;
            i_pairs[1] = i;
            j_pairs[1] = k;

            ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                    modelComputeArguments,
                    d2Eidr_three[5], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);
            }

            /* contribution to forces */
            if (comp_force)
            {
              for (kdim = 0; kdim < DIM; ++kdim)
              {
                /* accumulate force on atom i */
                force[i*DIM + kdim] += dEidr_three[0]*Rij[kdim]/R1 +
                                       dEidr_three[1]*Rik[kdim]/R2;
                /* accumulate force on atom j */
                force[j*DIM + kdim] -= dEidr_three[0]*Rij[kdim]/R1 -
                                       dEidr_three[2]*Rjk[kdim]/R4;
                /* accumulate force on atom k */
                force[k*DIM + kdim] -= dEidr_three[2]*Rjk[kdim]/R4 +
                                       dEidr_three[1]*Rik[kdim]/R2;
               }
            }

            /* Start adding four body terms */
            /*********************************/
            for (ll = kk+1; ll < numOfAtomNeigh; ++ll)
            {
              l = neighListOfCurrentAtom[ll]; /* get neighbor ID */
              lSpecies = particleSpeciesCodes[l];

              if (buffer->num_interactions == 16)
              {
                if (iSpecies == SPEC1 && jSpecies == SPEC1 && kSpecies == SPEC1 && lSpecies == SPEC1)
                  interaction_index = 0;
                else if (iSpecies == SPEC1 && jSpecies == SPEC1 && kSpecies == SPEC1 && lSpecies == SPEC2)
                  interaction_index = 1;
                else if (iSpecies == SPEC1 && jSpecies == SPEC1 && kSpecies == SPEC2 && lSpecies == SPEC1)
                  interaction_index = 2;
                else if (iSpecies == SPEC1 && jSpecies == SPEC1 && kSpecies == SPEC2 && lSpecies == SPEC2)
                  interaction_index = 3;
                else if (iSpecies == SPEC1 && jSpecies == SPEC2 && kSpecies == SPEC1 && lSpecies == SPEC1)
                  interaction_index = 4;
                else if (iSpecies == SPEC1 && jSpecies == SPEC2 && kSpecies == SPEC1 && lSpecies == SPEC2)
                  interaction_index = 5;
                else if (iSpecies == SPEC1 && jSpecies == SPEC2 && kSpecies == SPEC2 && lSpecies == SPEC1)
                  interaction_index = 6;
                else if (iSpecies == SPEC1 && jSpecies == SPEC2 && kSpecies == SPEC2 && lSpecies == SPEC2)
                  interaction_index = 7;
                else if (iSpecies == SPEC2 && jSpecies == SPEC1 && kSpecies == SPEC1 && lSpecies == SPEC1)
                  interaction_index = 8;
                else if (iSpecies == SPEC2 && jSpecies == SPEC1 && kSpecies == SPEC1 && lSpecies == SPEC2)
                  interaction_index = 9;
                else if (iSpecies == SPEC2 && jSpecies == SPEC1 && kSpecies == SPEC2 && lSpecies == SPEC1)
                  interaction_index = 10;
                else if (iSpecies == SPEC2 && jSpecies == SPEC1 && kSpecies == SPEC2 && lSpecies == SPEC2)
                  interaction_index = 11;
                else if (iSpecies == SPEC2 && jSpecies == SPEC2 && kSpecies == SPEC1 && lSpecies == SPEC1)
                  interaction_index = 12;
                else if (iSpecies == SPEC2 && jSpecies == SPEC2 && kSpecies == SPEC1 && lSpecies == SPEC2)
                  interaction_index = 13;
                else if (iSpecies == SPEC2 && jSpecies == SPEC2 && kSpecies == SPEC2 && lSpecies == SPEC1)
                  interaction_index = 14;
                else
                  interaction_index = 15;

                A         = &(buffer->A)[interaction_index];
                B         = &(buffer->B)[interaction_index];
                p         = &(buffer->p)[interaction_index];
                q         = &(buffer->q)[interaction_index];
                lambda    = &(buffer->lambda)[interaction_index];
                lambda_2  = &(buffer->lambda_2)[interaction_index];
                gamma     = &(buffer->gamma)[interaction_index];
                sigma     = &(buffer->sigma)[interaction_index];
                epsilon   = &(buffer->epsilon)[interaction_index];
                Q         = &(buffer->Q)[interaction_index];
                costhetat = &(buffer->costhetat)[interaction_index];
                a         = &(buffer->a)[interaction_index];
                cutsq     = &(buffer->cutsq)[interaction_index];
              }

              /* compute relative position vector and squared distance */
              Rsqil = 0.0;
              Rsqjl = 0.0;
              Rsqkl = 0.0;
              for (kdim = 0; kdim < DIM; ++kdim)
              {
                Ril[kdim] = coords[l*DIM + kdim] - coords[i*DIM + kdim];
                Rjl[kdim] = Ril[kdim] - Rij[kdim];
                Rkl[kdim] = Ril[kdim] - Rik[kdim];

                /* compute squared distance */
                Rsqil += Ril[kdim]*Ril[kdim];
                Rsqjl += Rjl[kdim]*Rjl[kdim];
                Rsqkl += Rkl[kdim]*Rkl[kdim];
              }

              /* compute energy and force */
              /* if (Rsqil > *cutsq || Rsqjl > *cutsq || Rsqkl > *cutsq) continue; */ /* particles are interacting ? */
              if (Rsqil > *cutsq) continue; /* particles are interacting ? */

              R3 = sqrt(Rsqil);
              R5 = sqrt(Rsqjl);
              R6 = sqrt(Rsqkl);

              if (comp_process_d2Edr2)
              {
                /* compute four-body potential and its derivatives */
                calc_phi_d2phi_four(a, lambda_2, gamma, sigma, epsilon, Q, costhetat,
                                    R1, R2, R3, R4, R5, R6, &phi_four, dphi_four, d2phi_four);

                /* compute dEidr */
                dEidr_four[0]   =  dphi_four[0];
                dEidr_four[1]   =  dphi_four[1];
                dEidr_four[2]   =  dphi_four[2];
                dEidr_four[3]   =  dphi_four[3];
                dEidr_four[4]   =  dphi_four[4];
                dEidr_four[5]   =  dphi_four[5];

                d2Eidr_four[0]  =  d2phi_four[0];
                d2Eidr_four[1]  =  d2phi_four[1];
                d2Eidr_four[2]  =  d2phi_four[2];
                d2Eidr_four[3]  =  d2phi_four[3];
                d2Eidr_four[4]  =  d2phi_four[4];
                d2Eidr_four[5]  =  d2phi_four[5];
                d2Eidr_four[6]  =  d2phi_four[6];
                d2Eidr_four[7]  =  d2phi_four[7];
                d2Eidr_four[8]  =  d2phi_four[8];
                d2Eidr_four[9]  =  d2phi_four[9];
                d2Eidr_four[10] =  d2phi_four[10];
                d2Eidr_four[11] =  d2phi_four[11];
                d2Eidr_four[12] =  d2phi_four[12];
                d2Eidr_four[13] =  d2phi_four[13];
                d2Eidr_four[14] =  d2phi_four[14];
                d2Eidr_four[15] =  d2phi_four[15];
                d2Eidr_four[16] =  d2phi_four[16];
                d2Eidr_four[17] =  d2phi_four[17];
                d2Eidr_four[18] =  d2phi_four[18];
                d2Eidr_four[19] =  d2phi_four[19];
                d2Eidr_four[20] =  d2phi_four[20];
              }
              else if (comp_force || comp_process_dEdr)
              {
                /* compute four-body potential and its derivative */
                calc_phi_dphi_four(a, lambda_2, gamma, sigma, epsilon, Q, costhetat,
                                   R1, R2, R3, R4, R5, R6, &phi_four, dphi_four);

                /* compute dEidr */
                dEidr_four[0]  =  dphi_four[0];
                dEidr_four[1]  =  dphi_four[1];
                dEidr_four[2]  =  dphi_four[2];
                dEidr_four[3]  =  dphi_four[3];
                dEidr_four[4]  =  dphi_four[4];
                dEidr_four[5]  =  dphi_four[5];
              }
              else
              {
                /* compute just four-body potential */
                calc_phi_four(a, lambda_2, gamma, sigma, epsilon, Q, costhetat,
                              R1, R2, R3, R4, R5, R6, &phi_four);
              }

              /* contribution to energy */
              if (comp_particleEnergy)
              {
                particleEnergy[i] += phi_four;
              }
              if (comp_energy)
              {
                 *energy += phi_four;
              }

              /* contribution to process_dEdr */
              if (comp_process_dEdr)
              {
                ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                        modelComputeArguments,
                        dEidr_four[0], R1, pRij, i, j);
                ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                        modelComputeArguments,
                        dEidr_four[1], R2, pRik, i, k);
                ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                        modelComputeArguments,
                        dEidr_four[2], R3, pRil, i, l);
                ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                        modelComputeArguments,
                        dEidr_four[3], R4, pRjk, j, k);
                ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                        modelComputeArguments,
                        dEidr_four[4], R5, pRjl, j, l);
                ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                        modelComputeArguments,
                        dEidr_four[5], R6, pRkl, k, l);
              }

              /* contribution to process_d2Edr2 */
              if (comp_process_d2Edr2)
              {
                R_pairs[0] = R_pairs[1] = R1;
                Rij_pairs[0][0] = Rij_pairs[1][0] = Rij[0];
                Rij_pairs[0][1] = Rij_pairs[1][1] = Rij[1];
                Rij_pairs[0][2] = Rij_pairs[1][2] = Rij[2];
                i_pairs[0] = i_pairs[1] = i;
                j_pairs[0] = j_pairs[1] = j;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[0], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R1;
                R_pairs[1] = R2;
                Rij_pairs[0][0] =  Rij[0];
                Rij_pairs[0][1] =  Rij[1];
                Rij_pairs[0][2] =  Rij[2];
                Rij_pairs[1][0] =  Rik[0];
                Rij_pairs[1][1] =  Rik[1];
                Rij_pairs[1][2] =  Rik[2];
                i_pairs[0]  = i;
                j_pairs[0]  = j;
                i_pairs[1]  = i;
                j_pairs[1]  = k;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[6], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R2;
                R_pairs[1] = R1;
                Rij_pairs[0][0] =  Rik[0];
                Rij_pairs[0][1] =  Rik[1];
                Rij_pairs[0][2] =  Rik[2];
                Rij_pairs[1][0] =  Rij[0];
                Rij_pairs[1][1] =  Rij[1];
                Rij_pairs[1][2] =  Rij[2];
                i_pairs[0]  = i;
                j_pairs[0]  = k;
                i_pairs[1]  = i;
                j_pairs[1]  = j;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[6], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R1;
                R_pairs[1] = R3;
                Rij_pairs[0][0] =  Rij[0];
                Rij_pairs[0][1] =  Rij[1];
                Rij_pairs[0][2] =  Rij[2];
                Rij_pairs[1][0] =  Ril[0];
                Rij_pairs[1][1] =  Ril[1];
                Rij_pairs[1][2] =  Ril[2];
                i_pairs[0]  = i;
                j_pairs[0]  = j;
                i_pairs[1]  = i;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[11], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R3;
                R_pairs[1] = R1;
                Rij_pairs[0][0] =  Ril[0];
                Rij_pairs[0][1] =  Ril[1];
                Rij_pairs[0][2] =  Ril[2];
                Rij_pairs[1][0] =  Rij[0];
                Rij_pairs[1][1] =  Rij[1];
                Rij_pairs[1][2] =  Rij[2];
                i_pairs[0]  = i;
                j_pairs[0]  = l;
                i_pairs[1]  = i;
                j_pairs[1]  = j;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[11], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R1;
                R_pairs[1] = R4;
                Rij_pairs[0][0] =  Rij[0];
                Rij_pairs[0][1] =  Rij[1];
                Rij_pairs[0][2] =  Rij[2];
                Rij_pairs[1][0] =  Rjk[0];
                Rij_pairs[1][1] =  Rjk[1];
                Rij_pairs[1][2] =  Rjk[2];
                i_pairs[0]  = i;
                j_pairs[0]  = j;
                i_pairs[1]  = j;
                j_pairs[1]  = k;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[15], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R4;
                R_pairs[1] = R1;
                Rij_pairs[0][0] =  Rjk[0];
                Rij_pairs[0][1] =  Rjk[1];
                Rij_pairs[0][2] =  Rjk[2];
                Rij_pairs[1][0] =  Rij[0];
                Rij_pairs[1][1] =  Rij[1];
                Rij_pairs[1][2] =  Rij[2];
                i_pairs[0]  = j;
                j_pairs[0]  = k;
                i_pairs[1]  = i;
                j_pairs[1]  = j;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[15], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R1;
                R_pairs[1] = R5;
                Rij_pairs[0][0] =  Rij[0];
                Rij_pairs[0][1] =  Rij[1];
                Rij_pairs[0][2] =  Rij[2];
                Rij_pairs[1][0] =  Rjl[0];
                Rij_pairs[1][1] =  Rjl[1];
                Rij_pairs[1][2] =  Rjl[2];
                i_pairs[0]  = i;
                j_pairs[0]  = j;
                i_pairs[1]  = j;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[18], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R5;
                R_pairs[1] = R1;
                Rij_pairs[0][0] =  Rjl[0];
                Rij_pairs[0][1] =  Rjl[1];
                Rij_pairs[0][2] =  Rjl[2];
                Rij_pairs[1][0] =  Rij[0];
                Rij_pairs[1][1] =  Rij[1];
                Rij_pairs[1][2] =  Rij[2];
                i_pairs[0]  = j;
                j_pairs[0]  = l;
                i_pairs[1]  = i;
                j_pairs[1]  = j;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[18], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R1;
                R_pairs[1] = R6;
                Rij_pairs[0][0] =  Rij[0];
                Rij_pairs[0][1] =  Rij[1];
                Rij_pairs[0][2] =  Rij[2];
                Rij_pairs[1][0] =  Rkl[0];
                Rij_pairs[1][1] =  Rkl[1];
                Rij_pairs[1][2] =  Rkl[2];
                i_pairs[0]  = i;
                j_pairs[0]  = j;
                i_pairs[1]  = k;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[20], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R6;
                R_pairs[1] = R1;
                Rij_pairs[0][0] =  Rkl[0];
                Rij_pairs[0][1] =  Rkl[1];
                Rij_pairs[0][2] =  Rkl[2];
                Rij_pairs[1][0] =  Rij[0];
                Rij_pairs[1][1] =  Rij[1];
                Rij_pairs[1][2] =  Rij[2];
                i_pairs[0]  = k;
                j_pairs[0]  = l;
                i_pairs[1]  = i;
                j_pairs[1]  = j;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[20], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R_pairs[1] = R2;
                Rij_pairs[0][0] = Rij_pairs[1][0] = Rik[0];
                Rij_pairs[0][1] = Rij_pairs[1][1] = Rik[1];
                Rij_pairs[0][2] = Rij_pairs[1][2] = Rik[2];
                i_pairs[0] = i_pairs[1] = i;
                j_pairs[0] = j_pairs[1] = k;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[1], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R2;
                R_pairs[1] = R3;
                Rij_pairs[0][0] =  Rik[0];
                Rij_pairs[0][1] =  Rik[1];
                Rij_pairs[0][2] =  Rik[2];
                Rij_pairs[1][0] =  Ril[0];
                Rij_pairs[1][1] =  Ril[1];
                Rij_pairs[1][2] =  Ril[2];
                i_pairs[0]  = i;
                j_pairs[0]  = k;
                i_pairs[1]  = i;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[7], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R3;
                R_pairs[1] = R2;
                Rij_pairs[0][0] =  Ril[0];
                Rij_pairs[0][1] =  Ril[1];
                Rij_pairs[0][2] =  Ril[2];
                Rij_pairs[1][0] =  Rik[0];
                Rij_pairs[1][1] =  Rik[1];
                Rij_pairs[1][2] =  Rik[2];
                i_pairs[0]  = i;
                j_pairs[0]  = l;
                i_pairs[1]  = i;
                j_pairs[1]  = k;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[7], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R2;
                R_pairs[1] = R4;
                Rij_pairs[0][0] =  Rik[0];
                Rij_pairs[0][1] =  Rik[1];
                Rij_pairs[0][2] =  Rik[2];
                Rij_pairs[1][0] =  Rjk[0];
                Rij_pairs[1][1] =  Rjk[1];
                Rij_pairs[1][2] =  Rjk[2];
                i_pairs[0]  = i;
                j_pairs[0]  = k;
                i_pairs[1]  = j;
                j_pairs[1]  = k;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[12], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R4;
                R_pairs[1] = R2;
                Rij_pairs[0][0] =  Rjk[0];
                Rij_pairs[0][1] =  Rjk[1];
                Rij_pairs[0][2] =  Rjk[2];
                Rij_pairs[1][0] =  Rik[0];
                Rij_pairs[1][1] =  Rik[1];
                Rij_pairs[1][2] =  Rik[2];
                i_pairs[0]  = j;
                j_pairs[0]  = k;
                i_pairs[1]  = i;
                j_pairs[1]  = k;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[12], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R2;
                R_pairs[1] = R5;
                Rij_pairs[0][0] =  Rik[0];
                Rij_pairs[0][1] =  Rik[1];
                Rij_pairs[0][2] =  Rik[2];
                Rij_pairs[1][0] =  Rjl[0];
                Rij_pairs[1][1] =  Rjl[1];
                Rij_pairs[1][2] =  Rjl[2];
                i_pairs[0]  = i;
                j_pairs[0]  = k;
                i_pairs[1]  = j;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[16], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R5;
                R_pairs[1] = R2;
                Rij_pairs[0][0] =  Rjl[0];
                Rij_pairs[0][1] =  Rjl[1];
                Rij_pairs[0][2] =  Rjl[2];
                Rij_pairs[1][0] =  Rik[0];
                Rij_pairs[1][1] =  Rik[1];
                Rij_pairs[1][2] =  Rik[2];
                i_pairs[0]  = j;
                j_pairs[0]  = l;
                i_pairs[1]  = i;
                j_pairs[1]  = k;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[16], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R2;
                R_pairs[1] = R6;
                Rij_pairs[0][0] =  Rik[0];
                Rij_pairs[0][1] =  Rik[1];
                Rij_pairs[0][2] =  Rik[2];
                Rij_pairs[1][0] =  Rkl[0];
                Rij_pairs[1][1] =  Rkl[1];
                Rij_pairs[1][2] =  Rkl[2];
                i_pairs[0]  = i;
                j_pairs[0]  = k;
                i_pairs[1]  = k;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[19], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R6;
                R_pairs[1] = R2;
                Rij_pairs[0][0] =  Rkl[0];
                Rij_pairs[0][1] =  Rkl[1];
                Rij_pairs[0][2] =  Rkl[2];
                Rij_pairs[1][0] =  Rik[0];
                Rij_pairs[1][1] =  Rik[1];
                Rij_pairs[1][2] =  Rik[2];
                i_pairs[0]  = k;
                j_pairs[0]  = l;
                i_pairs[1]  = i;
                j_pairs[1]  = k;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[19], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R_pairs[1] = R3;
                Rij_pairs[0][0] = Rij_pairs[1][0] = Ril[0];
                Rij_pairs[0][1] = Rij_pairs[1][1] = Ril[1];
                Rij_pairs[0][2] = Rij_pairs[1][2] = Ril[2];
                i_pairs[0] = i_pairs[1] = i;
                j_pairs[0] = j_pairs[1] = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[2], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R3;
                R_pairs[1] = R4;
                Rij_pairs[0][0] =  Ril[0];
                Rij_pairs[0][1] =  Ril[1];
                Rij_pairs[0][2] =  Ril[2];
                Rij_pairs[1][0] =  Rjk[0];
                Rij_pairs[1][1] =  Rjk[1];
                Rij_pairs[1][2] =  Rjk[2];
                i_pairs[0]  = i;
                j_pairs[0]  = l;
                i_pairs[1]  = j;
                j_pairs[1]  = k;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[8], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R4;
                R_pairs[1] = R3;
                Rij_pairs[0][0] =  Rjk[0];
                Rij_pairs[0][1] =  Rjk[1];
                Rij_pairs[0][2] =  Rjk[2];
                Rij_pairs[1][0] =  Ril[0];
                Rij_pairs[1][1] =  Ril[1];
                Rij_pairs[1][2] =  Ril[2];
                i_pairs[0]  = j;
                j_pairs[0]  = k;
                i_pairs[1]  = i;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[8], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R3;
                R_pairs[1] = R5;
                Rij_pairs[0][0] =  Ril[0];
                Rij_pairs[0][1] =  Ril[1];
                Rij_pairs[0][2] =  Ril[2];
                Rij_pairs[1][0] =  Rjl[0];
                Rij_pairs[1][1] =  Rjl[1];
                Rij_pairs[1][2] =  Rjl[2];
                i_pairs[0]  = i;
                j_pairs[0]  = l;
                i_pairs[1]  = j;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[13], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R5;
                R_pairs[1] = R3;
                Rij_pairs[0][0] =  Rjl[0];
                Rij_pairs[0][1] =  Rjl[1];
                Rij_pairs[0][2] =  Rjl[2];
                Rij_pairs[1][0] =  Ril[0];
                Rij_pairs[1][1] =  Ril[1];
                Rij_pairs[1][2] =  Ril[2];
                i_pairs[0]  = j;
                j_pairs[0]  = l;
                i_pairs[1]  = i;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[13], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R3;
                R_pairs[1] = R6;
                Rij_pairs[0][0] =  Ril[0];
                Rij_pairs[0][1] =  Ril[1];
                Rij_pairs[0][2] =  Ril[2];
                Rij_pairs[1][0] =  Rkl[0];
                Rij_pairs[1][1] =  Rkl[1];
                Rij_pairs[1][2] =  Rkl[2];
                i_pairs[0]  = i;
                j_pairs[0]  = l;
                i_pairs[1]  = k;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[17], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R6;
                R_pairs[1] = R3;
                Rij_pairs[0][0] =  Rkl[0];
                Rij_pairs[0][1] =  Rkl[1];
                Rij_pairs[0][2] =  Rkl[2];
                Rij_pairs[1][0] =  Ril[0];
                Rij_pairs[1][1] =  Ril[1];
                Rij_pairs[1][2] =  Ril[2];
                i_pairs[0]  = k;
                j_pairs[0]  = l;
                i_pairs[1]  = i;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[17], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R_pairs[1] = R4;
                Rij_pairs[0][0] = Rij_pairs[1][0] = Rjk[0];
                Rij_pairs[0][1] = Rij_pairs[1][1] = Rjk[1];
                Rij_pairs[0][2] = Rij_pairs[1][2] = Rjk[2];
                i_pairs[0] = i_pairs[1] = j;
                j_pairs[0] = j_pairs[1] = k;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[3], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R4;
                R_pairs[1] = R5;
                Rij_pairs[0][0] =  Rjk[0];
                Rij_pairs[0][1] =  Rjk[1];
                Rij_pairs[0][2] =  Rjk[2];
                Rij_pairs[1][0] =  Rjl[0];
                Rij_pairs[1][1] =  Rjl[1];
                Rij_pairs[1][2] =  Rjl[2];
                i_pairs[0]  = j;
                j_pairs[0]  = k;
                i_pairs[1]  = j;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[9], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R5;
                R_pairs[1] = R4;
                Rij_pairs[0][0] =  Rjl[0];
                Rij_pairs[0][1] =  Rjl[1];
                Rij_pairs[0][2] =  Rjl[2];
                Rij_pairs[1][0] =  Rjk[0];
                Rij_pairs[1][1] =  Rjk[1];
                Rij_pairs[1][2] =  Rjk[2];
                i_pairs[0]  = j;
                j_pairs[0]  = l;
                i_pairs[1]  = j;
                j_pairs[1]  = k;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[9], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R4;
                R_pairs[1] = R6;
                Rij_pairs[0][0] =  Rjk[0];
                Rij_pairs[0][1] =  Rjk[1];
                Rij_pairs[0][2] =  Rjk[2];
                Rij_pairs[1][0] =  Rkl[0];
                Rij_pairs[1][1] =  Rkl[1];
                Rij_pairs[1][2] =  Rkl[2];
                i_pairs[0]  = j;
                j_pairs[0]  = k;
                i_pairs[1]  = k;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[14], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R6;
                R_pairs[1] = R4;
                Rij_pairs[0][0] =  Rkl[0];
                Rij_pairs[0][1] =  Rkl[1];
                Rij_pairs[0][2] =  Rkl[2];
                Rij_pairs[1][0] =  Rjk[0];
                Rij_pairs[1][1] =  Rjk[1];
                Rij_pairs[1][2] =  Rjk[2];
                i_pairs[0]  = k;
                j_pairs[0]  = l;
                i_pairs[1]  = j;
                j_pairs[1]  = k;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[14], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R_pairs[1] = R5;
                Rij_pairs[0][0] = Rij_pairs[1][0] = Rjl[0];
                Rij_pairs[0][1] = Rij_pairs[1][1] = Rjl[1];
                Rij_pairs[0][2] = Rij_pairs[1][2] = Rjl[2];
                i_pairs[0] = i_pairs[1] = j;
                j_pairs[0] = j_pairs[1] = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[4], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R5;
                R_pairs[1] = R6;
                Rij_pairs[0][0] =  Rjl[0];
                Rij_pairs[0][1] =  Rjl[1];
                Rij_pairs[0][2] =  Rjl[2];
                Rij_pairs[1][0] =  Rkl[0];
                Rij_pairs[1][1] =  Rkl[1];
                Rij_pairs[1][2] =  Rkl[2];
                i_pairs[0]  = j;
                j_pairs[0]  = l;
                i_pairs[1]  = k;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[10], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R6;
                R_pairs[1] = R5;
                Rij_pairs[0][0] =  Rkl[0];
                Rij_pairs[0][1] =  Rkl[1];
                Rij_pairs[0][2] =  Rkl[2];
                Rij_pairs[1][0] =  Rjl[0];
                Rij_pairs[1][1] =  Rjl[1];
                Rij_pairs[1][2] =  Rjl[2];
                i_pairs[0]  = k;
                j_pairs[0]  = l;
                i_pairs[1]  = j;
                j_pairs[1]  = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[10], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

                R_pairs[0] = R_pairs[1] = R6;
                Rij_pairs[0][0] = Rij_pairs[1][0] = Rkl[0];
                Rij_pairs[0][1] = Rij_pairs[1][1] = Rkl[1];
                Rij_pairs[0][2] = Rij_pairs[1][2] = Rkl[2];
                i_pairs[0] = i_pairs[1] = k;
                j_pairs[0] = j_pairs[1] = l;

                ier = KIM_ModelComputeArguments_ProcessD2EDr2Term(
                        modelComputeArguments,
                        d2Eidr_four[5], pR_pairs, pRij_pairs, pi_pairs, pj_pairs);

              }
              /* contribution to forces */
              if (comp_force)
              {
                for (kdim = 0; kdim < DIM; ++kdim)
                {
                  /* accumulate force on atom i */
                  force[i*DIM + kdim] += dEidr_four[0]*Rij[kdim]/R1
                                         + dEidr_four[1]*Rik[kdim]/R2
                                         + dEidr_four[2]*Ril[kdim]/R3;
                  /* accumulate force on atom j */
                  force[j*DIM + kdim] += -dEidr_four[0]*Rij[kdim]/R1
                                         + dEidr_four[3]*Rjk[kdim]/R4
                                         + dEidr_four[4]*Rjl[kdim]/R5;
                  /* accumulate force on atom k */
                  force[k*DIM + kdim] += -dEidr_four[1]*Rik[kdim]/R2
                                         - dEidr_four[3]*Rjk[kdim]/R4
                                         + dEidr_four[5]*Rkl[kdim]/R6;
                  /* accumulate force on atom l */
                  force[l*DIM + kdim] += -dEidr_four[2]*Ril[kdim]/R3
                                         - dEidr_four[4]*Rjl[kdim]/R5
                                         - dEidr_four[5]*Rkl[kdim]/R6;
                 }
              }
            } /* loop on ll */
            /* End adding four body terms */
            /******************************/
        } /* loop on kk */

        /* End adding three body terms */
        /*******************************/
      } /* loop over neighbors jj of atom i */
    } /* Check on whether atom i is contributing */
  }

   /* perform final tasks */
   if (comp_particleEnergy && comp_energy)
   {
      *energy = 0.0;
      for (k = 0; k < *nAtoms; ++k)
      {
         *energy += particleEnergy[k];
      }
   }

   free(dphi_three);
   free(d2phi_three);
   free(dEidr_three);
   free(d2Eidr_three);

   free(dphi_four);
   free(d2phi_four);
   free(dEidr_four);
   free(d2Eidr_four);

   /* everything is great */
   ier = FALSE;

   return ier;
}

/* Reinitialization function */
/* We do not publish any parameters via "SetParameterPointer",
   so this function just trivially reregisters the existing
   influence distance in the API object */
#include "KIM_ModelRefreshLogMacros.h"
static int refresh(KIM_ModelRefresh * const modelRefresh)
{
   /* Local variables */
  int ier;
  struct model_buffer* buffer;

  /* get buffer from KIM object */
  LOG_INFORMATION("Retrieving Model buffer");
  KIM_ModelRefresh_GetModelBufferPointer(modelRefresh, (void **) &buffer);

  LOG_INFORMATION("Registering influence distance pointer");
  KIM_ModelRefresh_SetInfluenceDistancePointer(modelRefresh,
    &(buffer->influenceDistance));

  /* Register cutoff pointer */
  LOG_INFORMATION("Registering cutoff pointer");
  KIM_ModelRefresh_SetNeighborListPointers(modelRefresh, 1,
    &(buffer->influenceDistance), &(buffer->paddingNeighborHints),
    &(buffer->halfListHints));

  ier = FALSE;
  return ier;
}

/* destroy function */
#include "KIM_ModelDestroyLogMacros.h"
static int destroy(KIM_ModelDestroy * const modelDestroy)
{
   /* Local variables */
   struct model_buffer* buffer;
   int ier;

   /* get model buffer from KIM object */
   KIM_ModelDestroy_GetModelBufferPointer(modelDestroy, (void **) &buffer);

   /* free parameters */
   free(buffer->cutsq);
   free(buffer->A);
   free(buffer->B);
   free(buffer->p);
   free(buffer->q);
   free(buffer->a);
   free(buffer->lambda);
   free(buffer->lambda_2);
   free(buffer->gamma);
   free(buffer->sigma);
   free(buffer->epsilon);
   free(buffer->Q);
   free(buffer->costhetat);

   /* destroy the buffer */
   free(buffer);

   ier = FALSE;
   return ier;
}

/* compute_arguments create routine */
#include "KIM_ModelComputeArgumentsCreateLogMacros.h"
static int compute_arguments_create(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsCreate * const modelComputeArgumentsCreate)
{
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
  if (ier == TRUE)
  {
    LOG_ERROR("Unable to set argument support status");
  }

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
    LOG_ERROR("Unable to set callback support status");
  }
  return ier;
}

/* compute_arguments destroy routine */
#include "KIM_ModelComputeArgumentsDestroyLogMacros.h"
static int compute_arguments_destroy(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy)
{
  /* Nothing to be done */
  return FALSE;
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
  int num_interactions,
  double * const A,
  double * const sigma,
  double * const lambda,
  double * const lambda_2,
  double * const epsilon)
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
    for (i=0; i< num_interactions; ++i)
    {
    sigma[i] *= convertLength;
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
    for (i=0; i< num_interactions; ++i)
    {
      A[i] *= convertEnergy;
      lambda[i] *= convertEnergy;
      lambda_2[i] *= convertEnergy;
      epsilon[i] *= convertEnergy;
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


/* Initialization function */
#include "KIM_ModelDriverCreateLogMacros.h"
int model_driver_create(KIM_ModelDriverCreate *const modelDriverCreate,
                        KIM_LengthUnit const requestedLengthUnit,
                        KIM_EnergyUnit const requestedEnergyUnit,
                        KIM_ChargeUnit const requestedChargeUnit,
                        KIM_TemperatureUnit const requestedTemperatureUnit,
                        KIM_TimeUnit const requestedTimeUnit)
{
  /* KIM variables */
  const char* paramfile1name;
  char species1NameString[100], species2NameString[100];
  KIM_SpeciesName species1Name, species2Name;
  int numberOfParameterFiles;

  /* Local variables */
  FILE* fid;
  double* cutsq;
  double* A;
  double* B;
  double* p;
  double* q;
  double* a;
  double* lambda;
  double* lambda_2;
  double* gamma;
  double* sigma;
  double* epsilon;
  double* Q;
  double* costhetat;
  int ier;
  struct model_buffer* buffer;

  double max_cutoff;
  int num_species;
  int num_interactions;
  fpos_t filepos;
  char dummy[255];
  int i;

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
  KIM_ModelDriverCreate_SetRefreshPointer(modelDriverCreate,
                                          KIM_LANGUAGE_NAME_c,
                                          (func *) refresh);

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
  ier = KIM_ModelDriverCreate_GetParameterFileName(modelDriverCreate,
                                                   0,
                                                   &paramfile1name);
  if (ier == TRUE)
  {
    LOG_ERROR("Unable to get Model parameter file name.");
    return ier;
  }

  /* Read in model parameters from parameter file */
  fid = fopen(paramfile1name, "r");
  if (fid == NULL)
  {
     ier = TRUE;
     LOG_ERROR("Unable to open parameter file for "
               "Four_Body_Mistriotis_Flytzanis_Farantos parameters");
     return ier;
  }

  /* get rid of initial comments in the parameter file */
  fgetpos(fid, &filepos);
  fgets(dummy, 255, fid);
  while (dummy[0] == '#' || isspace(dummy[0]))
  {
    fgetpos(fid, &filepos);
    fgets(dummy, 255, fid);
  }
  fsetpos(fid, &filepos);

  /* read number of species */
  ier = fscanf(fid, "%d", &num_species);
  /* read number of species end */
  if (!(num_species == 1 || num_species == 2))
  {
     ier = TRUE;
     LOG_ERROR("Model parameter file must specify either 1 or 2 atomic species");
     return ier;
  }
  if (ier != 1)
  {
     ier = TRUE;
     LOG_ERROR("Error reading first line of parameter file");
     return ier;
  }

  if (num_species == 1)
  {
    num_interactions = 1;
    ier = fscanf(fid, "%s\n", species1NameString);
    if (ier != 1)
    {
      ier = TRUE;
      LOG_ERROR("Error reading species name from parameter file");
      return ier;
    }
    /* Register species */
    LOG_INFORMATION("Setting species code");
    species1Name = KIM_SpeciesName_FromString(species1NameString);
    ier = KIM_ModelDriverCreate_SetSpeciesCode(modelDriverCreate,
                                                      species1Name, SPEC1);
    if(ier == TRUE)
    {
      LOG_ERROR("Unable to set species code");
      return ier;
    }
  }
  else if (num_species == 2)
  {
    /* For two species in four body system, we have num_interactions
       = \sum_{i=0}^{i=4} 4Ci * (4-i)C(4-i) = 2^4 */
    num_interactions = 16;
    fgetpos(fid, &filepos);
    ier = fscanf(fid, "%s %s\n", &species1NameString, &species2NameString);
    if (ier != 2)
    {
      ier = TRUE;
      LOG_ERROR("Error reading species name from parameter file");
      return ier;
    }
    /* Register species */
    LOG_INFORMATION("Setting species codes");
    species1Name = KIM_SpeciesName_FromString(species1NameString);
    ier = KIM_ModelDriverCreate_SetSpeciesCode(modelDriverCreate,
                                                      species1Name, SPEC1);
    species2Name = KIM_SpeciesName_FromString(species2NameString);
    ier = KIM_ModelDriverCreate_SetSpeciesCode(modelDriverCreate,
                                                      species2Name, SPEC2);
    if(ier == TRUE)
    {
      LOG_ERROR("Unable to set species codes");
      return ier;
    }
  }

  /* Allocate buffer */
  LOG_INFORMATION("Allocating memory for Model buffer");
  buffer = (struct model_buffer*) malloc(sizeof(struct model_buffer));
  if (buffer == NULL)
  {
    ier = TRUE;
    LOG_ERROR("Could not allocate Model buffer");
    return ier;
  }

  /* Store model buffer in KIM object */
  LOG_INFORMATION("Registering Model buffer");
  KIM_ModelDriverCreate_SetModelBufferPointer(modelDriverCreate, (void*) buffer);

  /* Allocate memory for buffer parameters */
  cutsq     = (double*) malloc(num_interactions*sizeof(double));
  A         = (double*) malloc(num_interactions*sizeof(double));
  B         = (double*) malloc(num_interactions*sizeof(double));
  p         = (double*) malloc(num_interactions*sizeof(double));
  q         = (double*) malloc(num_interactions*sizeof(double));
  a         = (double*) malloc(num_interactions*sizeof(double));
  lambda    = (double*) malloc(num_interactions*sizeof(double));
  lambda_2  = (double*) malloc(num_interactions*sizeof(double));
  gamma     = (double*) malloc(num_interactions*sizeof(double));
  sigma     = (double*) malloc(num_interactions*sizeof(double));
  epsilon   = (double*) malloc(num_interactions*sizeof(double));
  Q         = (double*) malloc(num_interactions*sizeof(double));
  costhetat = (double*) malloc(num_interactions*sizeof(double));

  if(A==NULL || B==NULL        || p==NULL        || q==NULL     || a==NULL
             || lambda==NULL   || lambda_2==NULL || gamma==NULL || sigma==NULL
             || epsilon==NULL  || Q ==NULL       || costhetat==NULL)
  {
     ier = TRUE;
     LOG_ERROR("Failed to allocate local memory for parameters in model_driver_create");
     return ier;
  }

  /* read parameters */
  for (i=0; i < num_interactions; ++i)
  {
    /* get rid of comments begin */
    fgetpos(fid, &filepos);
    fgets(dummy, 255, fid);
    while (dummy[0] == '#' || isspace(dummy[0]))
    {
      fgetpos(fid, &filepos);
      fgets(dummy, 255, fid);
    }
    fsetpos(fid, &filepos);
    /* get rid of comments end */
    ier = fscanf(fid, "%lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n %lf\n",
                 &A[i],
                 &B[i],
                 &p[i],
                 &q[i],
                 &a[i],
                 &lambda[i],
                 &lambda_2[i],
                 &gamma[i],
                 &sigma[i],
                 &epsilon[i],
                 &Q[i],
                 &costhetat[i]);

    /* check that we read the right number of parameters */
    if (12 != ier)
    {
      ier = TRUE;
      LOG_ERROR("Unable to read all Four_Body_Mistriotis_Flytzanis_Farantos parameters");
      free(A);
      free(B);
      free(p);
      free(q);
      free(a);
      free(lambda);
      free(lambda_2);
      free(gamma);
      free(sigma);
      free(epsilon);
      free(Q);
      free(costhetat);
      return ier;
    }
  }
  fclose(fid);

  /* convert parameters to appropriate units (in-place) */
  ier = ConvertUnits(
          modelDriverCreate,
          requestedLengthUnit,
          requestedEnergyUnit,
          requestedChargeUnit,
          requestedTemperatureUnit,
          requestedTimeUnit,
          num_interactions,
          A, /* original units eV */
          sigma, /* original units of angstrom */
          lambda, /* original units eV */
          lambda_2, /* original units eV */
          epsilon /* original units of eV */
        );
  if (ier == TRUE)
  {
    LOG_ERROR("Failed to convert Model parameter units");
    free(A);
    free(B);
    free(p);
    free(q);
    free(a);
    free(lambda);
    free(lambda_2);
    free(gamma);
    free(sigma);
    free(epsilon);
    free(Q);
    free(costhetat);
    return ier;
  }

  /* We use only a single neighbor list here, so just set the
     influenceDistance and cutoff to the max over all interactions */
  max_cutoff = 0.0;
  for (i=0; i < num_interactions; ++i)
  {
    cutsq[i] = (a[i]*sigma[i])*(a[i]*sigma[i]);

    if (a[i]*sigma[i] > max_cutoff)
    {
      max_cutoff = a[i]*sigma[i];
    }
  }

  /* store parameters in buffer */
  buffer->influenceDistance = max_cutoff;
  buffer->cutsq             = cutsq;
  buffer->A                 = A;
  buffer->B                 = B;
  buffer->p                 = p;
  buffer->q                 = q;
  buffer->a                 = a;
  buffer->lambda            = lambda;
  buffer->lambda_2          = lambda_2;
  buffer->gamma             = gamma;
  buffer->sigma             = sigma;
  buffer->epsilon           = epsilon;
  buffer->Q                 = Q;
  buffer->costhetat         = costhetat;
  buffer->num_interactions  = num_interactions;

  /* Request that simulator omit neighbors of padding atoms */
  buffer->paddingNeighborHints = 1;

  /* Request full neighbor list from simulator */
  buffer->halfListHints = 0;

  /* Register influence distance pointer */
  LOG_INFORMATION("Registering influence distance pointer");
  KIM_ModelDriverCreate_SetInfluenceDistancePointer(modelDriverCreate,
    &(buffer->influenceDistance));

  /* Register cutoff pointer */
  LOG_INFORMATION("Registering cutoff pointer");
  KIM_ModelDriverCreate_SetNeighborListPointers(modelDriverCreate, 1,
    &(buffer->influenceDistance), &(buffer->paddingNeighborHints),
    &(buffer->halfListHints));

  ier = FALSE;
  return ier;
}


