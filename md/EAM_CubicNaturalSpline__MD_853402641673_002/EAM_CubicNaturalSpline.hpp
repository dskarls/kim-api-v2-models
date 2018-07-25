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
// Copyright (c) 2014--2018, Regents of the University of Minnesota.
// All rights reserved.
//
// Contributors:
//    Ryan S. Elliott
//    Mingjian Wen
//


#ifndef EAM_CUBIC_NATURAL_SPLINE_HPP_
#define EAM_CUBIC_NATURAL_SPLINE_HPP_

#define NUMBER_SPLINE_COEFF 9
#define D2F_LINEAR 0
#define D2F_CONSTANT 1
#define DF_QUADRATIC 2
#define DF_LINEAR 3
#define DF_CONSTANT 4
#define F_CUBIC 5
#define F_QUADRATIC 6
#define F_LINEAR 7
#define F_CONSTANT 8


//==============================================================================
//
// Definition of MACROs for improved efficiency
//
//==============================================================================

//******************************************************************************
// MACRO to compute parameters for cubic spline interpolation
// (used for efficiency)
//
// X - function argument
// H - 1/dX where dX is the spline knot spacing
// N - number of knots in the spline
// INDX = int(X/dX) * NUMBER_SPLINE_COEFF
// DELTAX = X/dX - int(X/dX)
#define GET_DELTAX_AND_INDEX(X, H, N, DELTAX, INDX)             \
  DELTAX  = std::max(X, 0.0) * H;                               \
  INDX    = static_cast<int>(DELTAX);                           \
  INDX    = std::min(INDX, N - 1);                              \
  DELTAX -= static_cast<double>(INDX);                          \
  INDX   *= NUMBER_SPLINE_COEFF;

//******************************************************************************
// MACRO to interpolate F(X) (used for efficiency)
//
// DX - DELTAX as computed in GET_DELTAX_AND_INDEX
// I  - INDX   as computed in GET_DELTAX_AND_INDEX
// F  - F(X)
#define INTERPOLATE_F(COEFF, DX, I, F)                    \
  F = COEFF[I + F_CUBIC] * DX + COEFF[I + F_QUADRATIC];   \
  F = F * DX + COEFF[I + F_LINEAR];                       \
  F = F * DX + COEFF[I + F_CONSTANT];

//******************************************************************************
// MACRO to interpolate dF(X)/dX (used for efficiency)
//
// DX - DELTAX as computed in GET_DELTAX_AND_INDEX
// I  - INDX   as computed in GET_DELTAX_AND_INDEX
// DF - dF(X)/dX
#define INTERPOLATE_DF(COEFF, DX, I, DF)                        \
  DF = COEFF[I + DF_QUADRATIC] * DX + COEFF[I + DF_LINEAR];     \
  DF = DF * DX + COEFF[I + DF_CONSTANT];

//******************************************************************************
// MACRO to interpolate d^2F(X)/dX^2 (used for efficiency)
//
// DX - DELTAX as computed in GET_DELTAX_AND_INDEX
// I  - INDX   as computed in GET_DELTAX_AND_INDEX
// D2F- d^2F(X)/dX^2
#define INTERPOLATE_D2F(COEFF, DX, I, D2F)                      \
  D2F = COEFF[I + D2F_LINEAR] * DX + COEFF[I + D2F_CONSTANT];


#endif  // EAM_CUBIC_NATURAL_SPLINE_HPP_
