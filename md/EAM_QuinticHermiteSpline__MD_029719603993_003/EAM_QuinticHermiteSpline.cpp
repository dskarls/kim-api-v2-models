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
#include "EAM_QuinticHermiteSpline.hpp"

//************************************************************************************
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

    // linear term and quadratic term are estimated using 7 points finite difference
    //linear coeff
    spline[0][F_LINEAR] =(-11*dat[0]+18*dat[1]-9*dat[2]+2*dat[3])/6;
    spline[1][F_LINEAR] =(-3*dat[0]-10*dat[1]+18*dat[2]-6*dat[3]+dat[4])/12;
    spline[2][F_LINEAR] =dat[0]/20-dat[1]/2-dat[2]/3+dat[3]-dat[4]/4+dat[5]/30;
    spline[n-3][F_LINEAR] =-dat[n-6]/30+dat[n-5]/4-dat[n-4]+dat[n-3]/3
   						   +dat[n-2]/2-dat[n-1]/20;
    spline[n-2][F_LINEAR] =(-dat[n-5]+6*dat[n-4]-18*dat[n-3]+10*dat[n-2]+3*dat[n-1])/12;
    spline[n-1][F_LINEAR] =(-2*dat[n-4]+9*dat[n-3]-18*dat[n-2]+11*dat[n-1])/6;

    for (int m = 3; m < n-3; ++m)
    {
        spline[m][F_LINEAR] = -dat[m-3]/60+3*dat[m-2]/20-3*dat[m-1]/4+3*dat[m+1]/4
        -3*dat[m+2]/20+dat[m+3]/60;
    }

    // quadratic coeff
    spline[0][F_QUADRATIC] =(2*dat[0]-5*dat[1]+4*dat[2]-dat[3])/2;
    spline[1][F_QUADRATIC] =((11*dat[0]-20*dat[1]+6*dat[2]+4*dat[3]-dat[4])/12)/2;
    spline[2][F_QUADRATIC] =(-dat[0]/12+4*dat[1]/3-5*dat[2]/2+4*dat[3]/3-dat[4]/12)/2;
    spline[n-3][F_QUADRATIC] =(-dat[n-5]/12+4*dat[n-4]/3-5*dat[n-3]/2+4*dat[n-2]/3
                               -dat[n-1]/12)/2;
    spline[n-2][F_QUADRATIC] =((-dat[n-5]+4*dat[n-4]+6*dat[n-3]-20*dat[n-2]
                                +11*dat[n-1])/12)/2;
    spline[n-1][F_QUADRATIC] =(-dat[n-4]+4*dat[n-3]-5*dat[n-2]+2*dat[n-1])/2;

    for (int m = 3; m < n-3; ++m)
    {
        spline[m][F_QUADRATIC] = (dat[m-3]/90-3*dat[m-2]/20+3*dat[m-1]/2-49*dat[m]/18
                                  +3*dat[m+1]/2-3*dat[m+2]/20+dat[m+3]/90)/2;
    }

    // cubic, quartic, and quintic coeff
    for (int m = 0; m < n-1; ++m)
    {
        spline[m][F_CUBIC]   = 10*spline[m+1][F_CONSTANT]-4*spline[m+1][F_LINEAR]
       						   +spline[m+1][F_QUADRATIC]-10*spline[m][F_CONSTANT]
      						   -6*spline[m][F_LINEAR]-3*spline[m][F_QUADRATIC];

        spline[m][F_QUARTIC] = -15*spline[m+1][F_CONSTANT]+7*spline[m+1][F_LINEAR]
     						   -2*spline[m+1][F_QUADRATIC]+15*spline[m][F_CONSTANT]
      						   +8*spline[m][F_LINEAR]+3*spline[m][F_QUADRATIC];

        spline[m][F_QUINTIC] = 6*spline[m+1][F_CONSTANT]-3*spline[m+1][F_LINEAR]
       						   +spline[m+1][F_QUADRATIC]-6*spline[m][F_CONSTANT]
     						   -3*spline[m][F_LINEAR]-spline[m][F_QUADRATIC];
    }
    spline[n-1][F_CUBIC]=0;
    spline[n-1][F_QUARTIC]=0;
    spline[n-1][F_QUINTIC]=0;

    //derivatives
    for (int m = 0; m < n-1; ++m)
    {
        spline[m][DF_CONSTANT]  = spline[m][F_LINEAR] / delta;
        spline[m][DF_LINEAR]    = 2.0 * spline[m][F_QUADRATIC] / delta;
        spline[m][DF_QUADRATIC] = 3.0 * spline[m][F_CUBIC] / delta;
        spline[m][DF_CUBIC]     = 4.0 * spline[m][F_QUARTIC] / delta;
        spline[m][DF_QUARTIC]   = 5.0 * spline[m][F_QUINTIC] / delta;
    }

    for (int m = 0; m < n-1; ++m)
    {
        spline[m][D2F_CONSTANT]  = spline[m][DF_LINEAR] / delta;
        spline[m][D2F_LINEAR]    = 2.0 * spline[m][DF_QUADRATIC] / delta;
        spline[m][D2F_QUADRATIC] = 3.0 * spline[m][DF_CUBIC] / delta;
        spline[m][D2F_CUBIC]     = 4.0 * spline[m][DF_QUARTIC] / delta;

    }

    delete [] spline;
}


