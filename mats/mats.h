//
// NOTE(michiel): This code was started based on the newlib code in mingw's source.
//

// NOTE(michiel): We don't raise errors, most signaling things are removed.


// TODO(michiel): Add function list
// NOTE(michiel): Sqrt
/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

// TODO(michiel): Add function list (complete it)
// NOTE(michiel): Sin/Cos
/*
   Copyright (c) 2018 Arm Ltd.  All rights reserved.

   SPDX-License-Identifier: BSD-3-Clause

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
   3. The name of the company may not be used to endorse or promote
      products derived from this software without specific prior written
      permission.

   THIS SOFTWARE IS PROVIDED BY ARM LTD ``AS IS'' AND ANY EXPRESS OR IMPLIED
   WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
   IN NO EVENT SHALL ARM LTD BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
   TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */


#include "mats_defines.h"
#include "mats_constants.h"

// NOTE(michiel): Abbreviations:
// inf         | infinity       (all exponent bits set and mantissa == 0)
// NaN         | Not a Number   (all exponent bits set and mantissa != 0)
// qNaN        | quiet NaN      (highest mantissa bit is 0)
// sNaN        | signaling NaN  (highest mantissa bit is 1)
// IVLD        | Invalid input
// DBZ         | Divide by zero

// NOTE(michiel): IEEE754 standard domain/error handling
//  Function   | Domain        | Exceptions             | Special (no exception, unless noted otherwise)
//  exp        | [-inf, inf]   | under/overflow         | exp(-inf) = 0, exp(inf) = inf
//  exp2       | [-inf, inf]   | under/overflow         | exp2(-inf) = 0, exp2(inf) = inf
//  exp10      | [-inf, inf]   | under/overflow         | exp10(-inf) = 0, exp10(inf) = inf
//  expm1      | [-inf, inf]   | under/overflow         | expm1(+0) = +0, expm1(-0) = -0, exp(-inf) = -1, exp(inf) = inf
//  log        | [0, inf]      | x = 0 => DBZ           | log(inf) = inf, log(+-0) = -inf with divideByZero,
//             |               | x < 0 => IVLD          | log(1) = 0
//  log2       | [0, inf]      | x = 0 => DBZ           | log2(inf) = inf, log2(+-0) = -inf with divideByZero,
//             |               | x < 0 => IVLD          | log2(1) = 0
//  log10      | [0, inf]      | x = 0 => DBZ           | log10(inf) = inf, log10(+-0) = -inf with divideByZero,
//             |               | x < 0 => IVLD          | log10(1) = 0
//  log1p      | [-1, inf]     | x = -1 => DBZ          | log1p(inf) = inf, log1p(-1) = -inf with divideByZero,
//             |               | x < -1 => IVLD         | log1p(+0) = +0, log1p(-0) = -0

//  pow        | [-inf, inf]x[-inf, inf] | pow(x, +-0)      = 1       for any x (even for zero/qNaN/inf)
//                                       | pow(+-0, y)      = +-inf   signals divideByZero for y an odd integer < 0
//                                       | pow(+-0, -inf)   = inf     no exceptions
//                                       | pow(+-0,  inf)   = +0      no exceptions
//                                       | pow(+-0, y)      = inf     signals divideByZero for finite y < 0 and not an odd integer
//                                       | pow(+-0, y)      = +-0     for finite y > 0 and odd integer
//                                       | pow(+-0, y)      = +0      for finite y > 0 and not odd integer
//                                       | pow(-1, +-inf)   = 1       no exceptions
//                                       | pow(+1, y)       = 1       for any y (even a qNaN)
//                                       | pow(x, y)                  signals invalid for finite x < 0 and finite non-integer y

// cos         | (-inf, inf)   | |x| = inf => IVLD      |
// sin         | (-inf, inf)   | |x| = inf => IVLD      | sin(+0) = +0, sin(-0) = -0
//             |               | underflow              |
// tan         | (-inf, inf)   | |x| = inf => IVLD      | tan(+0) = +0, tan(-0) = -0
//             |               | underflow              |
// acos        | [-1, 1]       | |x| > 1 => IVLD        | acos(1)  = +0
// asin        | [-1, 1]       | |x| > 1 => IVLD        | asin(+0) = +0, asin(-0) = -0
//             |               | underflow              |
// atan        | [-inf, inf]   | underflow              | atan(+0) = +0, atan(-0) = -0, atan(+-inf) = +-pi/2 and signals inexact

// atan2       | [-inf, inf]x[-inf, inf] | atan2(+y, x) for finite x > 0 is atan(|y/x|), which can signal inexact or underflow
//                                       | atan2(+y, x) for finite x < 0 is pi - atan(|y/x|), which can signal inexact
//                                       | atan2(+-0, -0)     = +-pi
//                                       | atan2(+-0, +0)     = +-0
//                                       | atan2(+-0, x)      = +-pi     for x < 0
//                                       | atan2(+-0, x)      = +-0      for x > 0
//                                       | atan2(y, +-0)      = -pi/2    for y < 0
//                                       | atan2(y, +-0)      = +pi/2    for y > 0
//                                       | atan2(+-y, -inf)   = +-pi     for finite y > 0
//                                       | atan2(+-y, +inf)   = +-0      for finite y > 0
//                                       | atan2(+-inf, x)    = +-pi/2   for finite x
//                                       | atan2(+-inf, -inf) = +-3pi/4
//                                       | atan2(+-inf, +inf) = +-pi/4

// cosh        | [-inf, inf]   | overflow               | cosh(+-inf) = +inf
// sinh        | [-inf, inf]   | under/overflow         | sinh(+0) = +0, sinh(-0) = -0, sinh(+-inf) = +-inf
// tanh        | [-inf, inf]   | underflow              | tanh(+0) = +0, tanh(-0) = -0, tanh(+-inf) = +-1
// acosh       | [1, inf]      | x < 1 => IVLD          | acosh(1) = +0, acosh(+inf) = +inf
// asinh       | [-inf, inf]   | underflow              | asinh(+0) = +0, asinh(-0) = -0
// atanh       | [-1, 1]       | underflow              | atanh(+0) = +0, atanh(-0) = -0,
//             |               | |x| = 1 => DBZ         | atanh(+-1) = +-inf and signals divideByZero
//             |               | |x| > 1 => IVLD        |

// TODO(michiel)
// - clean up all error returns
//   - replace all underflows with 0, maybe even treat subnormals as being 0
//   - replace all overflows with infinity
//   - replace divide by zero by infinity (hard coded)
//   - replace invalid by NaN (maybe set lowest bits with function id)
//   - need to handle NaN/inf and zero, especially on sse
//   - could ignore 0 signage
// - remove all libberdip stuff and add a optional file for type defines
// - 64bit

// NOTE(michiel): XX = 32 or 64
// NOTE(michiel): absoluteXX
#include "mats_common.h"
// NOTE(michiel): floorXX/ceilXX/roundXX/truncXX/modulusXX/remainderXX
#include "mats_rounding.h"
// NOTE(michiel): sqrtXX/hypotXX/expXX/exp2_XX/pow2_32/logXX/log2_XX/log10_XX/expm1_XX/log1pXX/log1p_fast32
#include "mats_elem.h"
// NOTE(michiel): pow32/pow10_32/exp10_32
#include "mats_elem_ext.h"
// NOTE(michiel): cos32/sin32/sincos32/tan32/acos32/asin32/atan32/atan2_32
//                cosh32/sinh32/sinhcosh32/tanh32/acosh32/asinh32/atanh32
#include "mats_trig.h"
