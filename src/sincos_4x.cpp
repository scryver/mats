// NOTE(michiel): This code is based on the arm sin/cos functions in mingw's newlib.

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

union f32_4x
{
    __m128  m;
    __m128i mi;
    f32 e[4];
    u32 u[4];
};

struct SinCos4x
{
    f32_4x cos;
    f32_4x sin;
};

internal __m128
sinf_exp_poly_q0_4x(__m128 x, __m128 x2)
{
    // NOTE(michiel): x - s1*x^3 + s2*x^5 - s3*x^7
    __m128 x3 = _mm_mul_ps(x, x2);
    __m128 s1 = _mm_sub_ps(_mm_set1_ps(0x1.110760p-7f),
                           _mm_mul_ps(x2, _mm_set1_ps(0x1.994eb2p-13f)));
    __m128 x7 = _mm_mul_ps(x2, x3);
    __m128 s  = _mm_sub_ps(x, _mm_mul_ps(x3, _mm_set1_ps(0x1.555544p-3f)));
    return _mm_add_ps(s, _mm_mul_ps(x7, s1));
}

internal __m128
sinf_exp_poly_q1_4x(__m128 x2)
{
    // NOTE(michiel): c0 - c1*x^2 + c2*x^4 - c3*x^6 + c4*x^8
    __m128 x4 = _mm_mul_ps(x2, x2);
    __m128 c2 = _mm_add_ps(_mm_set1_ps(-0x1.6c087ep-10f),
                           _mm_mul_ps(x2, _mm_set1_ps(0x1.993430p-16f)));
    __m128 c1 = _mm_sub_ps(_mm_set1_ps(1.0f),
                           _mm_mul_ps(x2, _mm_set1_ps(0x1.fffffep-2f)));
    __m128 x6 = _mm_mul_ps(x4, x2);
    __m128 c  = _mm_add_ps(c1, _mm_mul_ps(x4, _mm_set1_ps(0x1.55553ep-5f)));
    return _mm_add_ps(c, _mm_mul_ps(x6, c2));
}

internal f32_4x
arm_cosf_4x(f32_4x y)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    __m128 x = y.m;

    __m128 ones = _mm_set1_ps(1.0f);
    __m128 piOver4 = _mm_set1_ps(gPiOver4);
    __m128 minVal = _mm_set1_ps(0x1p-12f);

    __m128i absTopY = _mm_and_si128(y.mi, _mm_set1_epi32(0x7FF00000));
    __m128i absTopPiOver4 = _mm_and_si128(_mm_castps_si128(piOver4), _mm_set1_epi32(0x7FF00000));
    __m128i absMinVal = _mm_and_si128(_mm_castps_si128(minVal), _mm_set1_epi32(0x7FF00000));

    __m128 smallestMask = _mm_castsi128_ps(_mm_cmplt_epi32(absTopY, absMinVal));
    __m128 smallMask = _mm_castsi128_ps(_mm_cmplt_epi32(absTopY, absTopPiOver4));

    // NOTE(michiel): Reduce input argument
    __m128 hpiInv = _mm_set1_ps(0x1.45F306p-1f);
    __m128 r = _mm_mul_ps(x, hpiInv);
    __m128i n4xMod = _mm_cvtps_epi32(r);
    __m128  xMod   = _mm_sub_ps(x, _mm_mul_ps(_mm_round_ps(r, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC),
                                              _mm_set1_ps(0x1.921FB4p0f)));

    __m128 xm = _mm_blendv_ps(xMod, x, smallMask);
    __m128i n4x = _mm_andnot_si128(_mm_castps_si128(smallMask), n4xMod);

    __m128 x2 = _mm_mul_ps(xm, xm);

    __m128 sinV = sinf_exp_poly_q0_4x(xm, x2);
    __m128 cosV = sinf_exp_poly_q1_4x(x2);

    __m128i sel0 = _mm_and_si128(n4x, _mm_set1_epi32(1));
    __m128i sel1 = _mm_and_si128(n4x, _mm_set1_epi32(2));
    sel0 = _mm_cmpeq_epi32(sel0, _mm_set1_epi32(1));
    sel1 = _mm_cmpeq_epi32(sel1, _mm_set1_epi32(2));

    __m128 final = _mm_blendv_ps(cosV, sinV, _mm_castsi128_ps(sel0));

    __m128 signMask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    __m128 selSign = _mm_castsi128_ps(_mm_xor_si128(sel0, sel1));
    __m128 minResult = _mm_xor_ps(final, signMask);
    final = _mm_blendv_ps(final, minResult, selSign);
    final = _mm_blendv_ps(final, ones, smallestMask);

    f32_4x result;
    result.m = final;
    return result;
}

internal f32_4x
arm_sinf_4x(f32_4x y)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    __m128 x = y.m;

    __m128 piOver4 = _mm_set1_ps(gPiOver4);
    __m128 minVal  = _mm_set1_ps(0x1p-12f);

    __m128i absTopY = _mm_and_si128(y.mi, _mm_set1_epi32(0x7FF00000));
    __m128i absTopPiOver4 = _mm_and_si128(_mm_castps_si128(piOver4), _mm_set1_epi32(0x7FF00000));
    __m128i absMinVal = _mm_and_si128(_mm_castps_si128(minVal), _mm_set1_epi32(0x7FF00000));

    __m128 smallestMask = _mm_castsi128_ps(_mm_cmplt_epi32(absTopY, absMinVal));
    __m128 smallMask = _mm_castsi128_ps(_mm_cmplt_epi32(absTopY, absTopPiOver4));

    // NOTE(michiel): Reduce input argument
    __m128 hpiInv = _mm_set1_ps(0x1.45F306p-1f);
    __m128 r = _mm_mul_ps(x, hpiInv);
    __m128i n4xMod = _mm_cvtps_epi32(r);
    __m128 xMod    = _mm_sub_ps(x, _mm_mul_ps(_mm_round_ps(r, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC),
                                              _mm_set1_ps(0x1.921FB4p0f)));

    __m128 xm = _mm_blendv_ps(xMod, x, smallMask);
    __m128i n4x = _mm_andnot_si128(_mm_castps_si128(smallMask), n4xMod);

    __m128 x2 = _mm_mul_ps(xm, xm);

    __m128 sinV = sinf_exp_poly_q0_4x(xm, x2);
    __m128 cosV = sinf_exp_poly_q1_4x(x2);

    __m128i sel0 = _mm_and_si128(n4x, _mm_set1_epi32(1));
    __m128i sel1 = _mm_and_si128(n4x, _mm_set1_epi32(2));
    sel0 = _mm_cmpeq_epi32(sel0, _mm_set1_epi32(1));
    sel1 = _mm_cmpeq_epi32(sel1, _mm_set1_epi32(2));

    __m128 final = _mm_blendv_ps(sinV, cosV, _mm_castsi128_ps(sel0));

    __m128 signMask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    __m128 selSign = _mm_castsi128_ps(sel1);
    __m128 minResult = _mm_xor_ps(final, signMask);
    final = _mm_blendv_ps(final, minResult, selSign);
    final = _mm_blendv_ps(final, y.m, smallestMask);

    f32_4x result;
    result.m = final;
    return result;
}

internal SinCos4x
arm_sincosf_4x(f32_4x y)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    __m128 x = y.m;

    __m128 ones = _mm_set1_ps(1.0f);
    __m128 piOver4 = _mm_set1_ps(gPiOver4);
    __m128 minVal = _mm_set1_ps(0x1p-12f);

    __m128i absTopY = _mm_and_si128(y.mi, _mm_set1_epi32(0x7FF00000));
    __m128i absTopPiOver4 = _mm_and_si128(_mm_castps_si128(piOver4), _mm_set1_epi32(0x7FF00000));
    __m128i absMinVal = _mm_and_si128(_mm_castps_si128(minVal), _mm_set1_epi32(0x7FF00000));

    __m128 smallestMask = _mm_castsi128_ps(_mm_cmplt_epi32(absTopY, absMinVal));
    __m128 smallMask = _mm_castsi128_ps(_mm_cmplt_epi32(absTopY, absTopPiOver4));

    // NOTE(michiel): Reduce input argument
    __m128 hpiInv = _mm_set1_ps(0x1.45F306p-1f);
    __m128 r = _mm_mul_ps(x, hpiInv);
    __m128i n4xMod = _mm_cvtps_epi32(r);
    __m128  xMod   = _mm_sub_ps(x, _mm_mul_ps(_mm_round_ps(r, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC),
                                              _mm_set1_ps(0x1.921FB4p0f)));

    __m128 xm = _mm_blendv_ps(xMod, x, smallMask);
    __m128i n4x = _mm_andnot_si128(_mm_castps_si128(smallMask), n4xMod);

    __m128 x2 = _mm_mul_ps(xm, xm);

    __m128 sinV = sinf_exp_poly_q0_4x(xm, x2);
    __m128 cosV = sinf_exp_poly_q1_4x(x2);

    __m128i sel0 = _mm_and_si128(n4x, _mm_set1_epi32(1));
    __m128i sel1 = _mm_and_si128(n4x, _mm_set1_epi32(2));
    sel0 = _mm_cmpeq_epi32(sel0, _mm_set1_epi32(1));
    sel1 = _mm_cmpeq_epi32(sel1, _mm_set1_epi32(2));

    __m128 sin = _mm_blendv_ps(sinV, cosV, _mm_castsi128_ps(sel0));
    __m128 cos = _mm_blendv_ps(cosV, sinV, _mm_castsi128_ps(sel0));

    __m128 signMask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
    __m128 minSin = _mm_xor_ps(sin, signMask);
    sin = _mm_blendv_ps(sin, minSin, _mm_castsi128_ps(sel1));
    sin = _mm_blendv_ps(sin, y.m, smallestMask);

    __m128 selSign = _mm_castsi128_ps(_mm_xor_si128(sel0, sel1));
    __m128 minCos = _mm_xor_ps(cos, signMask);
    cos = _mm_blendv_ps(cos, minCos, selSign);
    cos = _mm_blendv_ps(cos, ones, smallestMask);

    SinCos4x result;
    result.cos.m = cos;
    result.sin.m = sin;
    return result;
}

#if 0
// NOTE(michiel): Above algorithms are based on these slimmed down versions
internal f32
arm_exp_cosf(f32 y)
{
    i_expect(absolute(y) < 120.0f);
    f32 x = y;

    f32 result;
    if (abstop12(y) < abstop12(0x1p-12f))
    {
        result = 1.0f;
    }
    else
    {
        int n = 0;
        f32 xm = (abstop12(y) < abstop12(gPiOver4)) ? x : reduce_fast_exp(x, &n);
        f32 x2 = xm * xm;
        switch (n & 3)
        {
            case 0: { result =  sinf_exp_poly_q1(x2); } break;
            case 1: { result = -sinf_exp_poly_q0(xm, x2); } break;
            case 2: { result = -sinf_exp_poly_q1(x2); } break;
            case 3: { result =  sinf_exp_poly_q0(xm, x2); } break;
        }
    }

    return result;
}

internal f32
arm_exp_sinf(f32 y)
{
    i_expect(absolute(y) < 120.0f);
    f32 x = y;

    f32 result;
    if (abstop12(y) < abstop12(0x1p-12f))
    {
        result = y;
    }
    else
    {
        int n = 0;
        f32 xm = (abstop12(y) < abstop12(gPiOver4)) ? x : reduce_fast_exp(x, &n);
        f32 x2 = xm * xm;
        switch (n & 3)
        {
            case 0: { result =  sinf_exp_poly_q0(xm, x2); } break;
            case 1: { result =  sinf_exp_poly_q1(x2); } break;
            case 2: { result = -sinf_exp_poly_q0(xm, x2); } break;
            case 3: { result = -sinf_exp_poly_q1(x2); } break;
        }
    }

    return result;
}
#endif
