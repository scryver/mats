/* Data definitions for sinf, cosf and sincosf.
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

#define TOINT_INTRINSICS 0

/* The constants and polynomials for sine and cosine.  */
struct ArmSinCos
{
    f64 sign[4];		/* Sign of sine in quadrants 0..3.  */
    f64 hpi_inv;		/* 2 / PI ( * 2^24 if !TOINT_INTRINSICS).  */
    f64 hpi;			/* PI / 2.  */
    f64 c0, c1, c2, c3, c4;	/* Cosine polynomial.  */
    f64 s1, s2, s3;		/* Sine polynomial.  */
};

/* The constants and polynomials for sine and cosine.  The 2nd entry
   computes -cos (x) rather than cos (x) to get negation for free.  */
global ArmSinCos gArmSinCosF32Table[2] =
{
    {
        { 1.0, -1.0, -1.0, 1.0 },
#if TOINT_INTRINSICS
        0x1.45F306DC9C883p-1,
#else
        0x1.45F306DC9C883p+23,
#endif
        0x1.921FB54442D18p0,
        0x1p0,
        -0x1.ffffffd0c621cp-2,
        0x1.55553e1068f19p-5,
        -0x1.6c087e89a359dp-10,
        0x1.99343027bf8c3p-16,
        -0x1.555545995a603p-3,
        0x1.1107605230bc4p-7,
        -0x1.994eb3774cf24p-13
    },
    {
        { 1.0, -1.0, -1.0, 1.0 },
#if TOINT_INTRINSICS
        0x1.45F306DC9C883p-1,
#else
        0x1.45F306DC9C883p+23,
#endif
        0x1.921FB54442D18p0,
        -0x1p0,
        0x1.ffffffd0c621cp-2,
        -0x1.55553e1068f19p-5,
        0x1.6c087e89a359dp-10,
        -0x1.99343027bf8c3p-16,
        -0x1.555545995a603p-3,
        0x1.1107605230bc4p-7,
        -0x1.994eb3774cf24p-13
    }
};

/* Table with 4/PI to 192 bit precision.  To avoid unaligned accesses
   only 8 new bits are added per entry, making the table 4 times larger.  */
global u32 gArm4OverPi[24] =
{
    0xa2,       0xa2f9,	  0xa2f983,   0xa2f9836e,
    0xf9836e4e, 0x836e4e44, 0x6e4e4415, 0x4e441529,
    0x441529fc, 0x1529fc27, 0x29fc2757, 0xfc2757d1,
    0x2757d1f5, 0x57d1f534, 0xd1f534dd, 0xf534ddc0,
    0x34ddc0db, 0xddc0db62, 0xc0db6295, 0xdb629599,
    0x6295993c, 0x95993c43, 0x993c4390, 0x3c439041
};
/* 2PI * 2^-64.  */
global const f64 g2PiPowMin64 = 0x1.921FB54442D18p-62;
/* PI / 4.  */
global const f64 gPiOver4 = 0x1.921FB54442D18p-1;

/* Top 12 bits of the float representation with the sign bit cleared.  */
func u32
abstop12(f32 x)
{
    //return (u32f32(x).u >> 20) & 0x7ff;
    return u32f32(x).u & 0x7ff00000;
}

func void
force_eval_float(f32 f)
{
    volatile f32 x = f;
    unused(x);
}

/* Compute the sine and cosine of inputs X and X2 (X squared), using the
   polynomial P and store the results in SINP and COSP.  N is the quadrant,
   if odd the cosine and sine polynomials are swapped.  */
func void
sincosf_poly(f64 x, f64 x2, ArmSinCos *p, int n, f32 *sinp, f32 *cosp)
{
    f64 x3, x4, x5, x6, s, c, c1, c2, s1;

    x4 = x2 * x2;
    x3 = x2 * x;
    c2 = p->c3 + x2 * p->c4;
    s1 = p->s2 + x2 * p->s3;

    /* Swap sin/cos result based on quadrant.  */
    f32 *tmp = (n & 1 ? cosp : sinp);
    cosp = (n & 1 ? sinp : cosp);
    sinp = tmp;

    c1 = p->c0 + x2 * p->c1;
    x5 = x3 * x2;
    x6 = x4 * x2;

    s = x + x3 * p->s1;
    c = c1 + x4 * p->c2;

    *sinp = s + x5 * s1;
    *cosp = c + x6 * c2;
}

/* Return the sine of inputs X and X2 (X squared) using the polynomial P.
   N is the quadrant, and if odd the cosine polynomial is used.  */
func f32
sinf_poly(f64 x, f64 x2, ArmSinCos *p, int n)
{
    f64 x3, x4, x6, x7, s, c, c1, c2, s1;

    if ((n & 1) == 0)
    {
        x3 = x * x2;
        s1 = p->s2 + x2 * p->s3;

        x7 = x3 * x2;
        s = x + x3 * p->s1;

        return s + x7 * s1;
    }
    else
    {
        x4 = x2 * x2;
        c2 = p->c3 + x2 * p->c4;
        c1 = p->c0 + x2 * p->c1;

        x6 = x4 * x2;
        c = c1 + x4 * p->c2;

        return c + x6 * c2;
    }
}

/* Fast range reduction using single multiply-subtract.  Return the modulo of
   X as a value between -PI/4 and PI/4 and store the quadrant in NP.
   The values for PI/2 and 2/PI are accessed via P.  Since PI/2 as a double
   is accurate to 55 bits and the worst-case cancellation happens at 6 * PI/4,
   the result is accurate for |X| <= 120.0.  */
func f64
reduce_fast(f64 x, ArmSinCos *p, int *np)
{
    f64 r;
#if TOINT_INTRINSICS
    /* Use fast round and lround instructions when available.  */
    r = x * p->hpi_inv;
    *np = converttoint(r);
    return x - roundtoint(r) * p->hpi;
#else
    /* Use scaled float to int conversion with explicit rounding.
       hpi_inv is prescaled by 2^24 so the quadrant ends up in bits 24..31.
       This avoids inaccuracies introduced by truncating negative values.  */
    r = x * p->hpi_inv;
    s32 n = ((s32)r + 0x800000) >> 24;
    *np = n;
    return x - n * p->hpi;
#endif
}

/* Reduce the range of XI to a multiple of PI/2 using fast integer arithmetic.
   XI is a reinterpreted float and must be >= 2.0f (the sign bit is ignored).
   Return the modulo between -PI/4 and PI/4 and store the quadrant in NP.
   Reduction uses a table of 4/PI with 192 bits of precision.  A 32x96->128 bit
   multiply computes the exact 2.62-bit fixed-point modulo.  Since the result
   can have at most 29 leading zeros after the binary point, the double
   precision result is accurate to 33 bits.  */
func f64
reduce_large(u32 xi, int *np)
{
    const u32 *arr = &gArm4OverPi[(xi >> 26) & 15];
    s32 shift = (xi >> 23) & 7;
    u64 n, res0, res1, res2;

    xi = (xi & 0xffffff) | 0x800000;
    xi <<= shift;

    res0 = xi * arr[0];
    res1 = (u64)xi * arr[4];
    res2 = (u64)xi * arr[8];
    res0 = (res2 >> 32) | (res0 << 32);
    res0 += res1;

    n = (res0 + (1ULL << 61)) >> 62;
    res0 -= n << 62;
    f64 x = (s64)res0;
    *np = n;
    return x * g2PiPowMin64;
}

/* Fast sincosf implementation.  Worst-case ULP is 0.5607, maximum relative
   error is 0.5303 * 2^-23.  A single-step range reduction is used for
   small values.  Large inputs have their range reduced using fast integer
   arithmetic.  */
func void
arm_sincosf(f32 y, f32 *sinp, f32 *cosp)
{
    f64 x = y;
    f64 s;
    int n;
    ArmSinCos *p = &gArmSinCosF32Table[0];

    if (abstop12(y) < abstop12(gPiOver4))
    {
        f64 x2 = x * x;

        if (abstop12(y) < abstop12(0x1p-12f))
        {
            if (abstop12(y) < abstop12(0x1p-126f))
            /* Force underflow for tiny y.  */
                force_eval_float(x2);
            *sinp = y;
            *cosp = 1.0f;
            return;
        }

        sincosf_poly(x, x2, p, 0, sinp, cosp);
    }
    else if (abstop12(y) < abstop12(120.0f))
    {
        x = reduce_fast(x, p, &n);

        /* Setup the signs for sin and cos.  */
        s = p->sign[n & 3];

        if (n & 2)
            p = &gArmSinCosF32Table[1];

        sincosf_poly(x * s, x * x, p, n, sinp, cosp);
    }
    else if (abstop12(y) < abstop12(F32_INF))
    {
        u32 xi = u32f32(y).u;
        s32 sign = xi >> 31;

        x = reduce_large(xi, &n);

        /* Setup signs for sin and cos - include original sign.  */
        s = p->sign[(n + sign) & 3];

        if ((n + sign) & 2)
            p = &gArmSinCosF32Table[1];

        sincosf_poly(x * s, x * x, p, n, sinp, cosp);
    }
    else
    {
        /* Return NaN if Inf or NaN for both sin and cos.  */
        *sinp = *cosp = y - y;
    }
}

/* Fast sinf implementation.  Worst-case ULP is 0.5607, maximum relative
   error is 0.5303 * 2^-23.  A single-step range reduction is used for
   small values.  Large inputs have their range reduced using fast integer
   arithmetic.  */
func f32
arm_sinf(f32 y)
{
    f64 x = y;
    f64 s;
    int n;
    ArmSinCos *p = &gArmSinCosF32Table[0];

    if (abstop12(y) < abstop12(gPiOver4))
    {
        s = x * x;

        if (abstop12(y) < abstop12(0x1p-12f))
        {
            if (abstop12(y) < abstop12(0x1p-126f))
            /* Force underflow for tiny y.  */
                force_eval_float (s);
            return y;
        }

        return sinf_poly(x, s, p, 0);
    }
    else if (abstop12(y) < abstop12(120.0f))
    {
        x = reduce_fast(x, p, &n);

        /* Setup the signs for sin and cos.  */
        s = p->sign[n & 3];

        if (n & 2)
            p = &gArmSinCosF32Table[1];

        return sinf_poly(x * s, x * x, p, n);
    }
    else if (abstop12(y) < abstop12(F32_INF))
    {
        u32 xi = u32f32(y).u;
        s32 sign = xi >> 31;

        x = reduce_large(xi, &n);

        /* Setup signs for sin and cos - include original sign.  */
        s = p->sign[(n + sign) & 3];

        if ((n + sign) & 2)
            p = &gArmSinCosF32Table[1];

        return sinf_poly(x * s, x * x, p, n);
    }
    else
        return y -y;
}

/* Fast cosf implementation.  Worst-case ULP is 0.5607, maximum relative
   error is 0.5303 * 2^-23.  A single-step range reduction is used for
   small values.  Large inputs have their range reduced using fast integer
   arithmetic.  */
func f32
arm_cosf(f32 y)
{
    f64 x = y;
    f64 s;
    int n;
    ArmSinCos *p = &gArmSinCosF32Table[0];

    if (abstop12(y) < abstop12(gPiOver4))
    {
        f64 x2 = x * x;

        if (abstop12(y) < abstop12(0x1p-12f))
            return 1.0f;

        return sinf_poly(x, x2, p, 1);
    }
    else if (abstop12 (y) < abstop12 (120.0f))
    {
        x = reduce_fast(x, p, &n);

        /* Setup the signs for sin and cos.  */
        s = p->sign[n & 3];

        if (n & 2)
            p = &gArmSinCosF32Table[1];

        return sinf_poly(x * s, x * x, p, n ^ 1);
    }
    else if (abstop12(y) < abstop12(F32_INF))
    {
        u32 xi = u32f32(y).u;
        s32 sign = xi >> 31;

        x = reduce_large(xi, &n);

        /* Setup signs for sin and cos - include original sign.  */
        s = p->sign[(n + sign) & 3];

        if ((n + sign) & 2)
            p = &gArmSinCosF32Table[1];

        return sinf_poly(x * s, x * x, p, n ^ 1);
    }
    else
        return y - y;
}
