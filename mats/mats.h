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

// TODO(michiel)
// - clean up all error returns
// - remove all libberdip stuff and add a optional file for type defines
// - 64bit

// NOTE(michiel): absolute32
#include "mats_common.h"
// NOTE(michiel): floor32/ceil32/round32/trunc32/modulus32/remainder32
#include "mats_rounding.h"
// NOTE(michiel): sqrt32/hypot32/exp32/exp2_32/pow2_32/log32/log2_32/log10_32/expm1_32/log1p32/log1p_fast32
#include "mats_elem.h"
// NOTE(michiel): pow32/pow10_32/exp10_32
#include "mats_elem_ext.h"
// NOTE(michiel): cos32/sin32/sincos32/tan32/acos32/asin32/atan32/atan2_32
//                cosh32/sinh32/sinhcosh32/tanh32/acosh32/asinh32/atanh32
#include "mats_trig.h"
