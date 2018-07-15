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
* Copyright (c) 2013,   Institute for Theoretical and Applied Physics
*            University of Stuttgart, D-70550 Stuttgart, Germany.
*       All rights reserved.
*
* Contributors:
*    Daniel Schopf
*
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "KIM_ModelDriverHeaders.h"
#include "IMD_EAM.h"

#define TRUE 1
#define FALSE 0

/* functions implementing the IMD EAM force routines */

/****************************************************************
 *
 *  Model compute function
 *
 ****************************************************************/
#include "KIM_ModelComputeLogMacros.h"
static int compute(KIM_ModelCompute const * const modelCompute,
                   KIM_ModelComputeArguments const * const modelComputeArguments){

  /* general variables */
  int ier;
  int i, j, k, jj, l;

  /* pointer to model buffer */
  model_buffer *buffer;

  /* flags for available computes */
  int comp_energy;
  int comp_force;
  int comp_particleEnergy;
  int comp_process_dEdr;

  /* pointers to objects in the KIM API */
  int *nAtoms;
  int *particleSpeciesCodes;
  int *particleContributing;
  double *coords;
  double *energy;
  double *force;
  double *particleEnergy;

  /* variables / pointers for objects in the buffer */
  int ntypes;
  double *dF_val;
  double *rho_val;

  /* variables for dealing with the KIM API */
  int const *neighListOfCurrentAtom;
  int numOfAtomNeigh;
  double Rij[DIM];
  double *pRij = &(Rij[0]);

  /* variables used for calculations */
  int col1, col2;
  int inc;
  int is_short = 0;
  int it, jt;
  double r2;
  double rho = 0.;
  double rho_i_prime, rho_j_prime;
  double R;
  double phi;
  double dphi;

  /* get buffer from KIM object */
  KIM_ModelCompute_GetModelBufferPointer(modelCompute, (void**) &buffer);

  /* unpack info from buffer */
  rho_val = buffer->rho_val;
  dF_val = buffer->dF_val;
  ntypes = buffer->ntypes;

  inc = ntypes * ntypes;

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
  if (ier)
  {
    LOG_ERROR("GetArgumentPointer");
    return ier;
  }

  /* check to see if we have been asked to compute the energy, forces, and particleEnergies */
  comp_energy = (energy != NULL);
  comp_force = (force != NULL);
  comp_particleEnergy = (particleEnergy != NULL);

  /* check to see if process_dEdr was requested */
  KIM_ModelComputeArguments_IsCallbackPresent(
     modelComputeArguments,
     KIM_COMPUTE_CALLBACK_NAME_ProcessDEDrTerm,
     &comp_process_dEdr);

  /* Check to be sure that the species are correct */
  ier = TRUE; /* assume an error */
  for (i = 0; i < *nAtoms; ++i)
  {
    if (particleSpeciesCodes[i] < 0 || particleSpeciesCodes[i]  > ntypes)
    {
      LOG_ERROR("Unexpected species code detected");
      return ier;
    }
  }
  ier = FALSE;  /* everything is ok */

  /* initialize potential energies, forces, and virial term */
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

  /* check if the arrays for rho and dF are big enough */
  if (buffer->table_len < *nAtoms) {
    buffer->rho_val = (double *)realloc(buffer->rho_val, (*nAtoms) * sizeof(double));
    buffer->dF_val = (double *)realloc(buffer->dF_val, (*nAtoms) * sizeof(double));
    if (NULL == buffer->rho_val || NULL == buffer->dF_val) {
      LOG_ERROR("Failed to reallocate arrays 'rho_val' and 'dF_val' in model buffer");
      return ier;
    }
    buffer->table_len = *nAtoms;
  }
  rho_val = buffer->rho_val;
  dF_val = buffer->dF_val;

  /* reset arrays */
  for (i = 0; i < *nAtoms; ++i)
  {
    rho_val[i] = 0.0;
    dF_val[i] = 0.0;
  }

  /* first loop over atoms, calculates pair energies and forces as well as eam energies */
  for (i = 0; i < *nAtoms; i++)
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

      it = particleSpeciesCodes[i];

      /* loop over all atoms in the neighbor list */
      for (jj = 0; jj < numOfAtomNeigh; jj++)
      {
        j = neighListOfCurrentAtom[jj];

        /* effective half list */
        if (i < j)
        {
          /* set up particle types and cols */
          jt = particleSpeciesCodes[j];
          col1 = it * ntypes + jt;
          col2 = jt * ntypes + it;

          /* calculate distance */
          r2 = 0.0;
          for (l = 0; l < DIM; l++)
          {
            Rij[l] = coords[j * DIM + l] - coords[i * DIM + l];
            r2 += Rij[l] * Rij[l];
          }
          R = sqrt(r2);

          /* check if we are within the cutoff radius */
          if (r2 <= buffer->pair_pot.end[col1])
          {
            /* calculate pair potential and the derivative at r2 */
            PAIR_INT(phi, dphi, buffer->pair_pot, col1, inc, r2, is_short);
            if (is_short) {
              short_dist_warning(0, i, j, particleSpeciesCodes, Rij, R);
              is_short = 0;
            }

            if (!particleContributing[j])
            {
              dphi *= 0.5;
            }

            if (comp_force){
              for (l = 0; l < DIM; l++) {
                force[i * DIM + l] += Rij[l] * dphi;
                force[j * DIM + l] -= Rij[l] * dphi;
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

            if (comp_particleEnergy)
            {
              particleEnergy[i] += 0.5*phi;
              if (particleContributing[j])
              {
                particleEnergy[j] += 0.5*phi;
              }
            }

            if (comp_process_dEdr)
            {
              dphi *= R;
              ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                      modelComputeArguments, dphi, R, pRij, i, j);
            }
          }

          /* calculate contribution to density */
          if (r2 < buffer->transfer_pot.end[col1]) {
            VAL_FUNC(rho, buffer->transfer_pot, col1, inc, r2, is_short);
            if (is_short)
            {
              short_dist_warning(1, i, j, particleSpeciesCodes, Rij, R);
              is_short = 0;
            }
            rho_val[i] += rho;
          }

          /* Also add density onto atom j if it is contributing */
          if (particleContributing[j])
          {
            if (jt == it)
            {
              if (r2 < buffer->transfer_pot.end[col2])
              {
                rho_val[j] += rho;
              }
            }
            else
            {
              if (r2 < buffer->transfer_pot.end[col2])
              {
                VAL_FUNC(rho, buffer->transfer_pot, col2, inc, r2, is_short);
                if (is_short)
                {
                  short_dist_warning(1, i, j, particleSpeciesCodes, Rij, R);
                  is_short = 0;
                }
                rho_val[j] += rho;
              }
            }
          }
        } /* if (i < j) */
      } /* loop over neighbors */

      /* calculate the embedding energies */
      PAIR_INT2(phi, dF_val[i], buffer->embed_pot, particleSpeciesCodes[i], ntypes, rho_val[i]);
      if (comp_energy) *energy += phi;
      if (comp_particleEnergy) particleEnergy[i] += phi;
    } /* check on whether atom i is contributing */
  } /* outer loop over atoms */


  /* second loop over atoms, calculates eam forces */
  for (i = 0; i < *nAtoms; i++)
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

      it = particleSpeciesCodes[i];

      /* loop over all atoms in the neighbor list */
      for (jj = 0; jj < numOfAtomNeigh; jj++)
      {
        j = neighListOfCurrentAtom[jj];

        /* effective half list */
        if (i < j)
        {
          /* set up particle types and cols */
          jt = particleSpeciesCodes[j];
          col1 = jt * ntypes + it;
          col2 = it * ntypes + jt;

          /* calculate distance */
          r2 = 0.0;
          for (l = 0; l < DIM; l++)
          {
            Rij[l] = coords[j * DIM + l] - coords[i * DIM + l];
            r2 += (Rij[l] * Rij[l]);
          }
          R = sqrt(r2);

          /* check if we are within the cutoff radius */
          if ((r2 < buffer->transfer_pot.end[col1]) || (r2 < buffer->transfer_pot.end[col2]))
          {
            /* calculate derivative of the density function for both atoms */
            DERIV_FUNC(rho_i_prime, buffer->transfer_pot, col1, inc, r2, is_short);
            if (is_short) {
              short_dist_warning(2, i, j, particleSpeciesCodes, Rij, R);
              is_short = 0;
            }
            if (col1 == col2) {
              rho_j_prime = rho_i_prime;
            }
            else
            {
              DERIV_FUNC(rho_j_prime, buffer->transfer_pot, col2, inc, r2, is_short);
              if (is_short)
              {
                short_dist_warning(2, i, j, particleSpeciesCodes, Rij, R);
                is_short = 0;
              }
            }

            /* combine all contributions */
            if (particleContributing[j])
            {
              dphi = 0.5 * (dF_val[i] * rho_j_prime + dF_val[j] * rho_i_prime);
            }
            else
            {
              dphi = 0.5 * dF_val[i] * rho_j_prime;
            }

            if (comp_force)
            {
              for (l = 0; l < DIM; l++)
              {
                force[i * DIM + l] += Rij[l] * dphi;
                force[j * DIM + l] -= Rij[l] * dphi;
              }
            }

            if (comp_process_dEdr)
            {
              dphi *= R;
              ier = KIM_ModelComputeArguments_ProcessDEDrTerm(
                      modelComputeArguments, dphi, R, pRij, i, j);
            }
          } /* check on cutoff distance */
        } /* if (i < j) */
      } /* loop over neighbors */
    } /* check on whether atom i is contributing */
  } /* loop over atoms */

  /* No errors */
  ier = FALSE;
  return ier;
}


/****************************************************************
 *
 * read potential table; choose format according to header
 *
 ****************************************************************/
#include "KIM_ModelDriverCreateLogMacros.h"
void read_pot_table(KIM_ModelDriverCreate * const modelDriverCreate,
  pot_table_t *pt, char *filename, int ncols, int ntypes, int radial)
{
  FILE *infile = NULL;
  char  buffer[1024], msg[255];
  char *res;
  int   have_header = 0, have_format = 0, end_header = 0;
  int   size = ncols;
  int   format = 2;    /* 2 for EAM2, 1 otherwise */
  int   i;

  /* read header */
  /* open file */
  infile = fopen(filename, "r");
  if (NULL == infile) {
    sprintf(msg, "Could not open potential file:\n\t\t %s", filename);
    LOG_ERROR(msg);
    exit(EXIT_FAILURE);
  }

  /* read the header */
  do
  {
    /* read one line */
    res = fgets(buffer, 1024, infile);
    if (NULL == res) {
      sprintf(msg, "Unexpected end of file in %s", filename);
      LOG_ERROR(msg);
      exit(EXIT_FAILURE);
    }
    /* see if it is a header line */
    if (buffer[0] == '#')
    {
      have_header = 1;
      /* stop after last header line */
      end_header = (buffer[1] == 'E');
      /* see if it is the format line */
      if (buffer[1] == 'F')
      {
        /* format complete? */
        if (2 != sscanf((const char *)(buffer + 2), "%d%d", &format, &size))
        {
          sprintf(msg, "Corrupted format header line in file %s", filename);
          LOG_ERROR(msg);
          exit(EXIT_FAILURE);
        }
        /* right number of columns? */
        if (size != ncols)
        {
          sprintf(msg, "Wrong number of data columns in file %%s\nShould be %d, is %d", ncols, size);
          LOG_ERROR(msg);
          exit(EXIT_FAILURE);
        }
        /* recognized format? */
        if ((format != 1) && (format != 2))
        {
          sprintf(msg, "Unrecognized format specified for file %s", filename);
          LOG_ERROR(msg);
          exit(EXIT_FAILURE);
        }
        have_format = 1;
      }
    }
    else if (have_header)
    {
      /* header does not end properly */
      sprintf(msg, "Corrupted header in file %s", filename);
      LOG_ERROR(msg);
      exit(EXIT_FAILURE);
    }
    else
    {
      /* we have no header, stop reading further */
      end_header = 1;
    }
  } while (!end_header);

  /* did we have a format in the header */
  if ((have_header) && (!have_format))
  {
    sprintf(msg, "Format not specified in header of file %s", filename);
    LOG_ERROR(msg);
    exit(EXIT_FAILURE);
  }

  /* rewind if there was no header */
  if (!have_header)
    rewind(infile);

  /* warn if we have no header */
  if (!have_header)
  {
    fprintf(stderr, "WARNING: File %s has no header !\n", filename);
    fflush(stderr);
  }

  /* have read header */

  /* allocate info block of function table */
  pt->maxsteps = 0;
  pt->ncols = ncols;
  pt->begin = (double *)malloc(ncols * sizeof(double));
  pt->end = (double *)malloc(ncols * sizeof(double));
  pt->step = (double *)malloc(ncols * sizeof(double));
  pt->invstep = (double *)malloc(ncols * sizeof(double));
  pt->len = (int *)malloc(ncols * sizeof(int));
  if ((pt->begin == NULL) || (pt->end == NULL) || (pt->step == NULL) || (pt->invstep == NULL)
    || (pt->len == NULL)) {
    sprintf(msg, "Cannot allocate info block for function table %s.", filename);
    LOG_ERROR(msg);
    exit(EXIT_FAILURE);
  }

  /* catch the case where potential is identically zero */
  for (i = 0; i < ncols; ++i)
  {
    pt->end[i] = 0.0;
    pt->len[i] = 0;
  }

  /* read table */
  if (format == 1)
    read_pot_table1(modelDriverCreate, pt, ncols, ntypes, filename, infile, radial);
  if (format == 2)
    read_pot_table2(modelDriverCreate, pt, ncols, ntypes, filename, infile, radial);
  fclose(infile);

  init_threepoint(pt, ncols);

  return;
}


/****************************************************************
 *
 * read potential in first format: each line contains
 *
 * r**2 V00 V01 V02 ... V10 V11 V12 ... VNN
 *
 * N is the number of different atom types
 *
 * Note that it is assumed that the r**2 are aequidistant.
 *
 ****************************************************************/
#include "KIM_ModelDriverCreateLogMacros.h"
void read_pot_table1(KIM_ModelDriverCreate * const modelDriverCreate,
  pot_table_t *pt, int ncols, int ntypes, char *filename, FILE *infile, int radial)
{
  char  msg[255];
  int   i, k;
  int   tablesize, npot = 0;
  double val, delta;
  double r2, r2_start = 0.0, r2_step;

  /* allocate the function table */
  pt->maxsteps = PSTEP;
  tablesize = ncols * pt->maxsteps;
  pt->table = (double *)malloc(tablesize * sizeof(double));
  if (NULL == pt->table)
  {
    sprintf(msg, "Cannot allocate memory for function table %s.", filename);
    LOG_ERROR(msg);
    exit(EXIT_FAILURE);
  }

  /* input loop */
  while (!feof(infile))
  {
    /* still some space left? */
    if (((npot % PSTEP) == 0) && (npot > 0))
    {
      pt->maxsteps += PSTEP;
      tablesize = ncols * pt->maxsteps;
      pt->table = (double *)realloc(pt->table, tablesize * sizeof(double));
      if (NULL == pt->table)
      {
        sprintf(msg, "Cannot extend memory for function table %s.", filename);
        LOG_ERROR(msg);
        exit(EXIT_FAILURE);
      }
    }

    /*  read in potential */
    if (1 != fscanf(infile, "%lf", &r2))
      break;
    if (npot == 0)
      r2_start = r2;    /* catch first value */
    for (i = 0; i < ncols; ++i)
    {
      if (1 != fscanf(infile, "%lf", &val))
      {
        LOG_ERROR("Line incomplete in potential file");
        exit(EXIT_FAILURE);
      }
      *PTR_2D(pt->table, npot, i, pt->maxsteps, ncols) = val;
      if (val != 0.)
      {    /* catch last non-zero value */
        pt->end[i] = r2;
        pt->len[i] = npot + 1;
      }
    }
    ++npot;
  }

  r2_step = (r2 - r2_start) / (npot - 1);

  /* fill info block, and shift potential to zero */
  for (i = 0; i < ncols; ++i)
  {
    pt->begin[i] = r2_start;
    pt->step[i] = r2_step;
    pt->invstep[i] = 1.0 / r2_step;
    delta = *PTR_2D(pt->table, (npot - 1), i, pt->maxsteps, ncols);
    /* if function of r2, shift potential and adjust cellsz */
    if (radial)
    {
      if (delta != 0.)
      {
        printf("Potential %1d%1d shifted by %f\n", (i / ntypes), (i % ntypes), delta);
        for (k = 0; k < npot; ++k)
          *PTR_2D(pt->table, k, i, pt->table, ncols) -= delta;
      }
    }
  }

  /* increase table size for security */
  tablesize = ncols * (pt->maxsteps + 2);
  pt->table = (double *)realloc(pt->table, tablesize * sizeof(double));
  if (NULL == pt->table)
  {
    sprintf(msg, "Cannot extend memory for function table %s.", filename);
    LOG_ERROR(msg);
    exit(EXIT_FAILURE);
  }

  return;
}


/****************************************************************
 *
 *  read potential in second format: at the beginning <ncols> times
 *  a line of the form
 *
 *  r_begin r_end r_step,
 *
 *  then the values of the potential (one per line), first those
 *  for atom pair 00, then an empty line (for gnuplot), then 01 and so on.
 *  Analogously, if there is only one column per atom type.
 *
 *  Note that it is assumed that the r**2 are aequidistant.
 *
 ****************************************************************/
#include "KIM_ModelDriverCreateLogMacros.h"
void read_pot_table2(KIM_ModelDriverCreate * const modelDriverCreate,
  pot_table_t *pt, int ncols, int ntypes, char *filename, FILE *infile, int radial)
{
  char  msg[255];
  int   i, k;
  int   tablesize;
  double val, numstep, delta;

  /* read the info block of the function table */
  for (i = 0; i < ncols; i++)
  {
    if (3 != fscanf(infile, "%lf %lf %lf", &pt->begin[i], &pt->end[i], &pt->step[i]))
    {
      sprintf(msg, "Info line %d in %s corrupt.", i + 1, filename);
      LOG_ERROR(msg);
      exit(EXIT_FAILURE);
    }
    pt->invstep[i] = 1.0 / pt->step[i];
    numstep = 1 + (pt->end[i] - pt->begin[i]) / pt->step[i];
    pt->len[i] = (int)(numstep + 0.49);
    pt->maxsteps = MAX(pt->maxsteps, pt->len[i]);

    /* some security against rounding errors */
    if (fabs(pt->len[i] - numstep) >= 0.1)
    {
      fprintf(stderr, "numstep = %f rounded to %d in file %s.\n", numstep, pt->len[i], filename);
      fflush(stderr);
    }
  }

  /* allocate the function table */
  /* allow some extra values at the end for interpolation */
  tablesize = ncols * (pt->maxsteps + 2);
  pt->table = (double *)malloc(tablesize * sizeof(double));
  if (NULL == pt->table)
  {
    sprintf(msg, "Cannot allocate memory for function table %s.", filename);
    LOG_ERROR(msg);
    exit(EXIT_FAILURE);
  }

  /* input loop */
  for (i = 0; i < ncols; i++)
  {
    for (k = 0; k < pt->len[i]; k++)
    {
      if (1 != fscanf(infile, "%lf", &val))
      {
        sprintf(msg, "wrong format in file %s.", filename);
        LOG_ERROR(msg);
        exit(EXIT_FAILURE);
      }
      *PTR_2D(pt->table, k, i, pt->maxsteps, ncols) = val;
    }
  }

  /* if function of r2, shift potential if necessary */
  if (radial)
  {
    for (i = 0; i < ncols; i++)
    {
      delta = *PTR_2D(pt->table, pt->len[i] - 1, i, pt->maxsteps, ncols);
      if (delta != 0.0)
      {
        printf("Potential %1d%1d shifted by %f\n", (i / ntypes), (i % ntypes), delta);
        for (k = 0; k < pt->len[i]; k++)
          *PTR_2D(pt->table, k, i, pt->table, ncols) -= delta;
      }
    }
  }

  return;
}


/****************************************************************
 *
 *  init_threepoint -- initialize for 3point interpolation
 *
 ****************************************************************/

void init_threepoint(pot_table_t *pt, int ncols)
{
  int   col, n;
  double *y;

  /* loop over columns */
  for (col = 0; col < ncols; col++)
  {

    y = pt->table + col;
    n = pt->len[col];

    /* for security, we continue the last interpolation polynomial */
    y[n * ncols] = 3 * y[(n - 1) * ncols] - 3 * y[(n - 2) * ncols] + y[(n - 3) * ncols];
    y[(n + 1) * ncols] = 6 * y[(n - 1) * ncols] - 8 * y[(n - 2) * ncols] + 3 * y[(n - 3) * ncols];

  }

  return;
}


/****************************************************************
 *
 * warning function for short distances
 * this is called every time a short distance is encountered
 *
 ****************************************************************/

void short_dist_warning(int type, int i, int j, int *types, double *Rij, double R)
{
  int   l;

  if (0 == type)
    fprintf(stderr, "Short distance in the pair potential!\n");
  else if (1 == type)
    fprintf(stderr, "Short distance in the transfer function!\n");
  else
    fprintf(stderr, "Short distance in the embedding function!\n");

  fprintf(stderr, "Involved particles are %d (type %d) and %d (type %d).\n", i, types[i], j, types[j]);
  fprintf(stderr, "Relative position vector is (");

  for (l = 0; l < (DIM - 1); l++)
    fprintf(stderr, "%f,", Rij[l]);

  fprintf(stderr, "%f), distance is %f\n\n", Rij[DIM - 1], R);
  fflush(stderr);

}


/****************************************************************
 *
 * model driver refresh function
 * this function is called by the KIM API after the potentials
 * have changed
 *
 * for IMD potentials this does not make sense, because we have no
 * published parameters which could be changed, so we provide a
 * dummy function which doesn't do anything (other than reregister
 * the cutoff and influence distance pointers, as required by the
 * KIM API)
 *
 ****************************************************************/
#include "KIM_ModelRefreshLogMacros.h"
static int refresh(KIM_ModelRefresh * const modelRefresh)
{
  /* Local variables */
  model_buffer *buffer;

  /* get model buffer from KIM object */
  LOG_INFORMATION("Getting model buffer");
  KIM_ModelRefresh_GetModelBufferPointer(modelRefresh,
                                         (void **) &buffer);

  LOG_INFORMATION("Resetting influence distance and cutoffs");
  KIM_ModelRefresh_SetInfluenceDistancePointer(
      modelRefresh, &(buffer->influenceDistance));
  KIM_ModelRefresh_SetNeighborListPointers(
      modelRefresh, 1,
      &(buffer->cutoff),
      &(buffer->paddingNeighborHints),
      &(buffer->halfListHints));

  /* No errors */
  return FALSE;
}


/****************************************************************
 *
 * model driver destroy function
 * this is called by the KIM API after the calculation is done to
 * allow us to free any memory we allocated
 *
 ****************************************************************/
int destroy(KIM_ModelDestroy * const modelDestroy)
{
  model_buffer *buffer;

  /* get model buffer from KIM object */
  KIM_ModelDestroy_GetModelBufferPointer(modelDestroy,
    (void**) &buffer);

  /* destroy all variables we allocated */
  free(buffer->pair_pot.begin);
  free(buffer->pair_pot.end);
  free(buffer->pair_pot.step);
  free(buffer->pair_pot.invstep);
  free(buffer->pair_pot.len);
  free(buffer->pair_pot.table);
  free(buffer->transfer_pot.begin);
  free(buffer->transfer_pot.end);
  free(buffer->transfer_pot.step);
  free(buffer->transfer_pot.invstep);
  free(buffer->transfer_pot.len);
  free(buffer->transfer_pot.table);
  free(buffer->embed_pot.begin);
  free(buffer->embed_pot.end);
  free(buffer->embed_pot.step);
  free(buffer->embed_pot.invstep);
  free(buffer->embed_pot.len);
  free(buffer->embed_pot.table);

  free(buffer->dF_val);
  free(buffer->rho_val);

  /* destroy the buffer */
  free(buffer);

  /* No errors */
  return FALSE;
}


/****************************************************************
 *
 * compute arguments create routine
 *
 ****************************************************************/
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

  /* register call backs */
  LOG_INFORMATION("Register call back supportStatus");
  error = error ||
      KIM_ModelComputeArgumentsCreate_SetCallbackSupportStatus(
          modelComputeArgumentsCreate,
          KIM_COMPUTE_CALLBACK_NAME_ProcessDEDrTerm,
          KIM_SUPPORT_STATUS_optional);

  if (error)
  {
    LOG_ERROR("Unable to successfully initialize compute arguments");
    return TRUE;
  }
  else
    return FALSE;
}

/****************************************************************
 *
 * compute arguments destroy routine
 *
 ****************************************************************/
#include "KIM_ModelComputeArgumentsDestroyLogMacros.h"
static int compute_arguments_destroy(
    KIM_ModelCompute const * const modelCompute,
    KIM_ModelComputeArgumentsDestroy * const modelComputeArgumentsDestroy)
{
  /* nothing to be done */

  return FALSE;
}


/****************************************************************
 *
 * model driver create function
 * this is called by the KIM API and reads the parameter file
 * of the model implemented by this model driver
 *
 ****************************************************************/
#include "KIM_ModelDriverCreateLogMacros.h"
int create(
    KIM_ModelDriverCreate * const modelDriverCreate,
    KIM_LengthUnit const requestedLengthUnit,
    KIM_EnergyUnit const requestedEnergyUnit,
    KIM_ChargeUnit const requestedChargeUnit,
    KIM_TemperatureUnit const requestedTemperatureUnit,
    KIM_TimeUnit const requestedTimeUnit)
{
  /* KIM variables */
  int numberOfParameterFiles;
  FILE *parameterFilePointers[3];

  /* local variables */
  char  msg[255], speciesNameString[100];
  KIM_SpeciesName speciesName;
  char const *species_paramfile_name, *pairpot_paramfile_name,
    *transfer_paramfile_name, *embed_paramfile_name;
  int   i, ier = FALSE;
  int   type, ntypes;
  double max_cutoff=0.0;
  model_buffer *buffer;
  FILE *infile;

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

  /* Check to make sure we have exactly three parameter files: one that
     simply contains the number of species and the order in which they
     appear in the remaining parameter files, one for the pair potential
     tabulation, one for the density function tabulation, and one for the
     embedding function tabulation. */
  if (numberOfParameterFiles != 4)
  {
    ier = TRUE;
    LOG_ERROR("Incorrect number of parameter files.");
    return ier;
  }
  /* get parameter file names */
  ier = KIM_ModelDriverCreate_GetParameterFileName(
          modelDriverCreate,
          0,
          &species_paramfile_name);
  if (ier)
  {
    LOG_ERROR("Unable to get species parameter file name.");
    return ier;
  }
  ier = KIM_ModelDriverCreate_GetParameterFileName(
          modelDriverCreate,
          1,
          &pairpot_paramfile_name);
  if (ier)
  {
    LOG_ERROR("Unable to get pair potential parameter file name.");
    return ier;
  }
  ier = KIM_ModelDriverCreate_GetParameterFileName(
          modelDriverCreate,
          2,
          &transfer_paramfile_name);
  if (ier)
  {
    LOG_ERROR("Unable to get density function parameter file name.");
    return ier;
  }
  ier = KIM_ModelDriverCreate_GetParameterFileName(
          modelDriverCreate,
          3,
          &embed_paramfile_name);
  if (ier)
  {
    LOG_ERROR("Unable to get embedding function parameter file name.");
    return ier;
  }

  /* read the species parameter file and set particleSpeciesCodes */
  infile = fopen(species_paramfile_name, "r");
  if (NULL == infile)
  {
    ier = TRUE;
    sprintf(msg,"Unable to open species parameter file:\n\t\t %s",
      species_paramfile_name);
    LOG_ERROR(msg);
    return ier;
  }
  ier = fscanf(infile, "%d\n", &ntypes);
  if (ier != 1)
  {
    ier = TRUE;
    sprintf(msg,"Could not read number of species types from parameter file:\n\t\t %s",
      species_paramfile_name);
    LOG_ERROR(msg);
    return ier;
  }
  for (type=0; type<ntypes; ++type)
  {
    ier = fscanf(infile, "%s\n", speciesNameString);
    if (ier != 1)
    {
      ier = TRUE;
      sprintf(msg,"Could not read all species types from parameter file:\n\t\t %s",
        species_paramfile_name);
      LOG_ERROR(msg);
      return ier;
    }
    /* Get the KIM API's species code using this string */
    speciesName = KIM_SpeciesName_FromString(speciesNameString);
    ier = KIM_ModelDriverCreate_SetSpeciesCode(
            modelDriverCreate,
            speciesName,
            type);
    if (ier == TRUE)
    {
      LOG_ERROR("Unable to set species code.");
      return ier;
    }
  }

  /* allocate buffer */
  buffer = (model_buffer*) malloc(sizeof(model_buffer));
  if (NULL == buffer)
  {
    ier = TRUE;
    LOG_ERROR("malloc");
    return ier;
  }

  /* Set ntypes in model buffer */
  buffer->ntypes = ntypes;

  /* Register model buffer pointer in KIM API */
  KIM_ModelDriverCreate_SetModelBufferPointer(modelDriverCreate,
    (void*) buffer);

  /* read the tabulated pair potential file */
  read_pot_table(modelDriverCreate, &(buffer->pair_pot), pairpot_paramfile_name,
    ntypes * ntypes, ntypes, 1);
  /* read the tabulated electron density function */
  read_pot_table(modelDriverCreate, &(buffer->transfer_pot), transfer_paramfile_name,
    ntypes * ntypes, ntypes, 1);
  /* read the tabulated embedding energy function */
  read_pot_table(modelDriverCreate, &(buffer->embed_pot), embed_paramfile_name,
    ntypes, ntypes, 0);

  /* calculate the cutoff */
  /* the cutoff is the maximum cutoff of all potentials
   * this is needed so the model can give us all needed neighbors
   * the compute routines check for the cutoff of the different potentials
   */
  for (i = 0; i < ntypes * ntypes; i++)
    max_cutoff = MAX(max_cutoff, buffer->pair_pot.end[i]);
  for (i = 0; i < ntypes * ntypes; i++)
    max_cutoff = MAX(max_cutoff, buffer->transfer_pot.end[i]);
  for (i = 0; i < ntypes; i++)
    max_cutoff = MAX(max_cutoff, buffer->embed_pot.end[i]);
  /* the distance in IMD is r^2, we need the cutoff in angstrom */
  max_cutoff = sqrt(max_cutoff);

  buffer->influenceDistance = max_cutoff;
  buffer->cutoff = max_cutoff;

  /* store model cutoffs in KIM object */
  KIM_ModelDriverCreate_SetInfluenceDistancePointer(
      modelDriverCreate,
      &(buffer->influenceDistance));
  KIM_ModelDriverCreate_SetNeighborListPointers(
      modelDriverCreate, 1,
      &(buffer->influenceDistance),
      &(buffer->paddingNeighborHints),
      &(buffer->halfListHints));

  /* allocate the arrays for the density and embedding values */
  buffer->rho_val = (double *)malloc(1 * sizeof(double));
  buffer->dF_val = (double *)malloc(1 * sizeof(double));
  if (NULL == buffer->rho_val || NULL == buffer->dF_val) {
    ier = TRUE;
    LOG_ERROR("Failed to allocate memory for 'rho_val' and 'dF_val' arrays in model buffer");
    return ier;
  }
  /* Initial value for table_len (set in compute function) */
  buffer->table_len = 1;

  /* Request omission of neighbors of padding atoms if possible */
  buffer->paddingNeighborHints = 1;

  /* Request half list if possible */
  buffer->halfListHints = 1;

  /* No errors */
  ier = FALSE;
  return ier;
}
