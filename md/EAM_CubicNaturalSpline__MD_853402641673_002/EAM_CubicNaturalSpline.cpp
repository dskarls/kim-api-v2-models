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
// Copyright (c) 2014, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Ryan S. Elliott
//    Mingjian Wen
//

#include "EAM_Implementation.hpp"
#include "EAM_CubicNaturalSpline.hpp"

//******************************************************************************
/* QC natural cubic spline interpolation. */
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

  // local variables
  double* const f7  = new double[n];
  double* const CSu = new double[n];
  double CSp;
  double CSqn;
  double CSun;

  // calculate the second derivatives
  f7[0] = 0.0;
  CSu[0] = 0.0;
  for (int i = 1; i < n-1; i++)
  {
    CSp = 0.5*f7[i-1]+2.0;
    f7[i] = -0.5/CSp;
    CSu[i] = (dat[i+1]-2.0*dat[i]+dat[i-1])/delta;
    CSu[i] = (3.0*CSu[i]/delta-0.5*CSu[i-1])/CSp;
  }
  CSqn = 0.0;
  CSun = 0.0;
  f7[n-1] = (CSun-CSqn*CSu[n-2])/(CSqn*f7[n-2]+1.0);
  for (int i = n-2; i >= 0; i--)
  {
    f7[i] = f7[i]*f7[i+1]+CSu[i];
  }

  // calculate the interpolation coefficients
  for (int m = 0; m < n-1; ++m)
  {
    spline[m][F_CONSTANT] = dat[m];
    spline[m][F_LINEAR] = dat[m+1]-dat[m]-delta*delta*(f7[m+1]+2.0*f7[m])/6.0;
    spline[m][F_QUADRATIC] = delta*delta*f7[m]/2.0;
    spline[m][F_CUBIC] = delta*delta*(f7[m+1]-f7[m])/6.0;
  }
  spline[n-1][F_LINEAR] = 0.0;
  spline[n-1][F_QUADRATIC] = 0.0;
  spline[n-1][F_CUBIC] = 0.0;

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

  delete [] f7;
  delete [] CSu;
  delete [] spline;
}
