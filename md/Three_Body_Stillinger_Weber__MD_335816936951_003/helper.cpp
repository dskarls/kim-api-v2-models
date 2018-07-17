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

#include "helper.hpp"

void ProcessVirialTerm(double const dEidr, double const rij,
    double const* const r_ij, int const i, int const j, VectorOfSizeSix virial)
{
  double const v = dEidr / rij;

  virial[0] += v * r_ij[0] * r_ij[0];
  virial[1] += v * r_ij[1] * r_ij[1];
  virial[2] += v * r_ij[2] * r_ij[2];
  virial[3] += v * r_ij[1] * r_ij[2];
  virial[4] += v * r_ij[0] * r_ij[2];
  virial[5] += v * r_ij[0] * r_ij[1];
}


void ProcessParticleVirialTerm(double const dEidr, double const rij,
    double const* const r_ij, int const i, int const j,
    VectorOfSizeSix* const particleVirial)
{
  double const v = dEidr / rij;
  VectorOfSizeSix vir;

  vir[0] = 0.5 * v * r_ij[0] * r_ij[0];
  vir[1] = 0.5 * v * r_ij[1] * r_ij[1];
  vir[2] = 0.5 * v * r_ij[2] * r_ij[2];
  vir[3] = 0.5 * v * r_ij[1] * r_ij[2];
  vir[4] = 0.5 * v * r_ij[0] * r_ij[2];
  vir[5] = 0.5 * v * r_ij[0] * r_ij[1];

  for (int k = 0; k < 6; ++k) {
    particleVirial[i][k] += vir[k];
    particleVirial[j][k] += vir[k];
  }
}
