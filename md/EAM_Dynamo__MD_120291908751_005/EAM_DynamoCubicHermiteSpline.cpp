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
// Copyright (c) 2017, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Ryan S. Elliott
//


#include "EAM_Implementation.hpp"
#include "EAM_DynamoCubicHermiteSpline.hpp"

//******************************************************************************
void EAM_Implementation::SplineInterpolate(double const* const dat,
                                           double const delta,
                                           int const n,
                                           double* const coe)
{ // setup convenient pointers (spline) into the coefficients (coe) array
  double** const spline = new double*[n];  // deleted at end of function
  for (int i = 0; i < n; ++i)
  {
    spline[i] = &coe[i * NUMBER_SPLINE_COEFF];
  }

  for (int m = 0; m < n; ++m)
  {
    spline[m][F_CONSTANT] = dat[m];
  }

  // Parts of this function originally obtained under CDDL from Steve Plimpton
  spline[0][F_LINEAR] = spline[1][F_CONSTANT] - spline[0][F_CONSTANT];
  spline[1][F_LINEAR] = HALF * (spline[2][F_CONSTANT] - spline[0][F_CONSTANT]);
  spline[n - 2][F_LINEAR] = HALF * (spline[n - 1][F_CONSTANT] -
                                    spline[n - 3][F_CONSTANT]);
  spline[n - 1][F_LINEAR] = spline[n - 1][F_CONSTANT]
      - spline[n - 2][F_CONSTANT];

  for (int m = 2; m <= n - 3; ++m)
  {
    spline[m][F_LINEAR] = ((spline[m - 2][F_CONSTANT] -
                            spline[m + 2][F_CONSTANT]) +
                           8.0 * (spline[m + 1][F_CONSTANT] -
                                  spline[m - 1][F_CONSTANT])) / 12.0;
  }

  for (int m = 0; m <= n - 2; ++m)
  {
    spline[m][F_QUADRATIC] = 3.0 * (spline[m + 1][F_CONSTANT] -
                                    spline[m][F_CONSTANT]) -
        2.0 * spline[m][F_LINEAR] - spline[m + 1][F_LINEAR];
    spline[m][F_CUBIC] = spline[m][F_LINEAR] + spline[m + 1][F_LINEAR] -
        2.0 * (spline[m + 1][F_CONSTANT] - spline[m][F_CONSTANT]);
  }

  spline[n - 1][F_QUADRATIC] = 0.0;
  spline[n - 1][F_CUBIC] = 0.0;

  for (int m = 0; m < n; ++m)
  {
    spline[m][DF_CONSTANT] = spline[m][F_LINEAR] / delta;
    spline[m][DF_LINEAR] = 2.0 * spline[m][F_QUADRATIC] / delta;
    spline[m][DF_QUADRATIC] = 3.0 * spline[m][F_CUBIC] / delta;
  }

  for (int m = 0; m < n; ++m)
  {
    spline[m][D2F_CONSTANT] = spline[m][DF_LINEAR] / delta;
    spline[m][D2F_LINEAR] = 2.0 * spline[m][DF_QUADRATIC] / delta;
  }

  delete [] spline;
}
