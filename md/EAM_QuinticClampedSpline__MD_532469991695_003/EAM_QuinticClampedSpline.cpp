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
#include "EAM_QuinticClampedSpline.hpp"

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
        
    // B.C.s end first and second derivatives  four point finite difference    
    spline[0][F_LINEAR] =(-11*dat[0]+18*dat[1]-9*dat[2]+2*dat[3])/6;
    spline[0][F_QUADRATIC] =(2*dat[0]-5*dat[1]+4*dat[2]-dat[3])/2;
    spline[n-1][F_LINEAR] =(-2*dat[n-4]+9*dat[n-3]-18*dat[n-2]+11*dat[n-1])/6;
    spline[n-1][F_QUADRATIC] =(-dat[n-4]+4*dat[n-3]-5*dat[n-2]+2*dat[n-1])/2;
        
    //local variables    
    double * a =new double [2*(n-2)];
    double * b =new double [2*(n-2)];
    double * c =new double [2*(n-2)];
    double * d =new double [2*(n-2)];
    double * e =new double [2*(n-2)];
    double * f =new double [2*(n-2)];
    double * g =new double [2*(n-2)];
    double * h =new double [2*(n-2)];
    double * xn=new double [2*(n-2)];
    double xmult;    
    
    // put the matrix values into the variables 
    for (int i=0; i<n-2;++i)
    {
        d[2*i]=6;
        d[2*i+1]=16;
        
        c[2*i]=0;
        c[2*i+1]=-4;
        e[2*i]=0;
        e[2*i+1]=-2;
        
        b[2*i]=-1;
        b[2*i+1]=7;
        f[2*i]=-1;
        f[2*i+1]=7;
        
        a[2*i]=2;
        a[2*i+1]=0;
        g[2*i]=4;
        g[2*i+1]=0;
        
        h[2*i]=10*(dat[i+2]-2*dat[i+1]+dat[i]);
        h[2*i+1]=15*(dat[i+2]-dat[i]);
    }
    
    h[0]=h[0]+spline[0][F_QUADRATIC]+4*spline[0][F_LINEAR];
    h[1]=h[1]-2*spline[0][F_QUADRATIC]-7*spline[0][F_LINEAR];
    h[2*(n-2)-2]=h[2*(n-2)-2]+spline[n-1][F_QUADRATIC]-4*spline[n-1][F_LINEAR];
    h[2*(n-2)-1]=h[2*(n-2)-1]+2*spline[n-1][F_QUADRATIC]-7*spline[n-1][F_LINEAR];
       
    // gaussian forward elimination 
    for (int i=1;i<2*(n-2)-2; ++i)
    {
        xmult = c[i-1]/d[i-1];
        d[i]=d[i]-xmult*e[i-1];
        e[i]=e[i]-xmult*f[i-1];
        f[i]=f[i]-xmult*g[i-1];
        h[i]=h[i]-xmult*h[i-1]; 
        
        xmult = b[i-1]/d[i-1];
        c[i]=c[i]-xmult*e[i-1];
        d[i+1]=d[i+1]-xmult*f[i-1];
        e[i+1]=e[i+1]-xmult*g[i-1];
        h[i+1]=h[i+1]-xmult*h[i-1];
                
        xmult = a[i-1]/d[i-1];
        b[i]=b[i]-xmult*e[i-1];
        c[i+1]=c[i+1]-xmult*f[i-1];
        d[i+2]=d[i+2]-xmult*g[i-1];
        h[i+2]=h[i+2]-xmult*h[i-1];
    }
    
    // elimination of the last two lines
    xmult = c[2*(n-2)-3]/d[2*(n-2)-3];
    d[2*(n-2)-2]=d[2*(n-2)-2]-xmult*e[2*(n-2)-3];
    e[2*(n-2)-2]=e[2*(n-2)-2]-xmult*f[2*(n-2)-3];
    h[2*(n-2)-2]=h[2*(n-2)-2]-xmult*h[2*(n-2)-3];    
    
    xmult = b[2*(n-2)-3]/d[2*(n-2)-3];
    c[2*(n-2)-2]=c[2*(n-2)-2]-xmult*e[2*(n-2)-3];
    d[2*(n-2)-1]=d[2*(n-2)-1]-xmult*f[2*(n-2)-3];
    h[2*(n-2)-1]=h[2*(n-2)-1]-xmult*h[2*(n-2)-3];
    
    //eliminaiton of the last line
    xmult = c[2*(n-2)-2]/d[2*(n-2)-2];
    d[2*(n-2)-1]=d[2*(n-2)-1]-xmult*e[2*(n-2)-2];
    h[2*(n-2)-1]=h[2*(n-2)-1]-xmult*h[2*(n-2)-2];    
       
    // gaussian backward substitution   
    xn[2*(n-2)-1]=h[2*(n-2)-1]/d[2*(n-2)-1];
    xn[2*(n-2)-2]=(h[2*(n-2)-2]-e[2*(n-2)-2]*xn[2*(n-2)-1])/d[2*(n-2)-2];
    xn[2*(n-2)-3]=(h[2*(n-2)-3]-e[2*(n-2)-3]*xn[2*(n-2)-2]
                   -f[2*(n-2)-3]*xn[2*(n-2)-1])/d[2*(n-2)-3];
    
    for (int i= 2*(n-2)-4;i>=0;--i)
    {
        xn[i]=(h[i]-e[i]*xn[i+1]-f[i]*xn[i+2]-g[i]*xn[i+3])/d[i];
    }
    
    // put the coefficinets into coe
    for (int i=1; i<n-1;i++)
    {
        spline[i][F_QUADRATIC]=xn[2*i-2];
        spline[i][F_LINEAR]=xn[2*i-1];
    }
        
    for (int m = 0; m < n-1; ++m)
    {
        spline[m][F_CONSTANT] = dat[m];
        
        spline[m][F_CUBIC] =10*(dat[m+1]-dat[m])
        -6*spline[m][F_LINEAR]-4*spline[m+1][F_LINEAR]
        -3*spline[m][F_QUADRATIC]+spline[m+1][F_QUADRATIC];
        
        spline[m][F_QUARTIC]=-15*(dat[m+1]-dat[m])
        +8*spline[m][F_LINEAR]+7*spline[m+1][F_LINEAR]
        +3*spline[m][F_QUADRATIC]-2*spline[m+1][F_QUADRATIC];
        
        
        spline[m][F_QUINTIC]=6*(dat[m+1]-dat[m])
        -3*spline[m][F_LINEAR]-3*spline[m+1][F_LINEAR]
        -spline[m][F_QUADRATIC]+spline[m+1][F_QUADRATIC];
    }    
    
    for (int m = 0; m < n-1; ++m)
    {
        spline[m][DF_CONSTANT] = spline[m][F_LINEAR] / delta;
        spline[m][DF_LINEAR] = 2.0 * spline[m][F_QUADRATIC] / delta;
        spline[m][DF_QUADRATIC] = 3.0 * spline[m][F_CUBIC] / delta;
        spline[m][DF_CUBIC] = 4.0 * spline[m][F_QUARTIC] / delta;
        spline[m][DF_QUARTIC] = 5.0 * spline[m][F_QUINTIC] / delta;
    }
    
    for (int m = 0; m < n-1; ++m)
    {
        spline[m][D2F_CONSTANT] = spline[m][DF_LINEAR] / delta;
        spline[m][D2F_LINEAR] = 2.0 * spline[m][DF_QUADRATIC] / delta;
        spline[m][D2F_QUADRATIC] = 3.0 * spline[m][DF_CUBIC] / delta;
        spline[m][D2F_CUBIC] = 4.0 * spline[m][DF_QUARTIC] / delta;
    }
   
    delete []a;
    delete []b;
    delete []c;
    delete []d;
    delete []e;
    delete []f;
    delete []h;
    delete []xn;
    delete [] spline;   
}

