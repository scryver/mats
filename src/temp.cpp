union f32_4x
{
    __m128  m;
    __m128i mi;
    f32 e[4];
    u32 u[4];
};

func f32_4x
F32_4x(f32 value)
{
    f32_4x result;
    result.m = _mm_set1_ps(value);
    return result;
}

func f32_4x
S32_4x(s32 value)
{
    f32_4x result;
    result.m = _mm_set1_epi32(value);
    return result;
}

func f32_4x
F32_4x(f32 value0, f32 value1, f32 value2, f32 value3)
{
    f32_4x result;
    result.m = _mm_setr_ps(value0, value1, value2, value3);
    return result;
}

func f32_4x
zero_4x(void)
{
    f32_4x result;
    result.m = _mm_setzero_ps();
    return result;
}

func f32_4x
operator +(f32_4x a, f32_4x b)
{
    f32_4x result;
    result.m = _mm_add_ps(a.m, b.m);
    return result;
}

func f32_4x
operator -(f32_4x a, f32_4x b)
{
    f32_4x result;
    result.m = _mm_sub_ps(a.m, b.m);
    return result;
}

func f32_4x
operator *(f32_4x a, f32_4x b)
{
    f32_4x result;
    result.m = _mm_mul_ps(a.m, b.m);
    return result;
}

func f32_4x
operator &(f32_4x a, f32_4x b)
{
    f32_4x result;
    result.m = _mm_and_ps(a.m, b.m);
    return result;
}

func f32_4x
operator |(f32_4x a, f32_4x b)
{
    f32_4x result;
    result.m = _mm_or_ps(a.m, b.m);
    return result;
}

func f32_4x
andnot(f32_4x a, f32_4x b)
{
    // NOTE(michiel): A & ~B
    f32_4x result;
    result.m = _mm_andnot_ps(b.m, a.m);
    return result;
}

func f32_4x
select(f32_4x a, f32_4x mask, f32_4x b)
{
    // NOTE(michiel): If mask then b else a, ( (A & ~M) | (B & M) )
    f32_4x result;
    result = andnot(a, mask) | (mask & b);
    return result;
}

func f32_4x
mix_in(f32_4x a, f32_4x mask, f32_4x b)
{
    // NOTE(michiel): If mask then b with a, ( A | (B & M) )
    f32_4x result;
    result = a | (mask & b);
    return result;
}

func f32_4x
cmplt_s32_4x(f32_4x a, f32_4x b)
{
    f32_4x result;
    result.mi = _mm_cmplt_epi32(a.mi, b.mi);
    return result;
}

func f32_4x
cmpeq_s32_4x(f32_4x a, f32_4x b)
{
    f32_4x result;
    result.mi = _mm_cmpeq_epi32(a.mi, b.mi);
    return result;
}

/* Top 12 bits of the float representation with the sign bit cleared.  */
func f32_4x
abstop12_4x(f32_4x x)
{
    //return (u32f32(x).u >> 20) & 0x7ff;
    f32_4x result = x;
    result.mi = _mm_and_si128(result.mi, _mm_set1_epi32(0x7FF00000));
    return result;
}

func __m128
reduce_fast_exp_4x(__m128 x, __m128i *np)
{
    __m128 hpiInv = _mm_set1_ps(0x1.45F306p-1f);
    __m128 r = _mm_mul_ps(x, hpiInv);
    _mm_storeu_si128(np, _mm_cvtps_epi32(r));
    return _mm_sub_ps(x, _mm_mul_ps(_mm_round_ps(r, _MM_FROUND_TO_NEAREST_INT |_MM_FROUND_NO_EXC),
                                    _mm_set1_ps(0x1.921FB4p0f)));
}

func __m128
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

func __m128
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

func f32_4x
operator -(f32_4x f4)
{
    u32 signMask = (u32)(1 << 31);
    __m128 flipMask = _mm_set1_ps(*(f32 *)&signMask);
    f32_4x result;
    result.m = _mm_xor_ps(f4.m, flipMask);
    return result;
}

func f32_4x
arm_cosf_4x(f32_4x y)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    f32_4x x = y;

    f32_4x ones = F32_4x(1.0f);
    f32_4x piOver4 = F32_4x(gPiOver4);
    f32_4x minVal = F32_4x(0x1p-12f);

    f32_4x absTopY = abstop12_4x(y);
    f32_4x absTopPiOver4 = abstop12_4x(piOver4);
    f32_4x absMinVal = abstop12_4x(minVal);

    f32_4x smallestMask = cmplt_s32_4x(absTopY, absMinVal);
    f32_4x smallMask    = cmplt_s32_4x(absTopY, absTopPiOver4);

    f32_4x n4xMod;
    f32_4x xMod;
    xMod.m = reduce_fast_exp_4x(x.m, &n4xMod.mi);

    f32_4x xm = select(xMod, smallMask, x);
    f32_4x n4x = andnot(n4xMod, smallMask);
    //n4x.mi = _mm_andnot_si128(smallMask.mi, n4xMod.mi);

    f32_4x x2 = xm * xm;

    f32_4x sinV;
    f32_4x cosV;
    sinV.m = sinf_exp_poly_q0_4x(xm.m, x2.m);
    cosV.m = sinf_exp_poly_q1_4x(x2.m);

    // TODO(michiel): Maybe pick out bit 0 and bit 1 and do some or-ing instead of cmp-ing
    f32_4x selection = n4x & S32_4x(0x3);
    f32_4x sel0Mask  = cmpeq_s32_4x(selection, zero_4x());
    f32_4x sel1Mask  = cmpeq_s32_4x(selection, S32_4x(1));
    f32_4x sel2Mask  = cmpeq_s32_4x(selection, S32_4x(2));;
    f32_4x sel3Mask  = cmpeq_s32_4x(selection, S32_4x(3));

    f32_4x result = sel0Mask & cosV;
    result = mix_in(result, sel1Mask, -sinV);
    result = mix_in(result, sel2Mask, -cosV);
    result = mix_in(result, sel3Mask, sinV);
    result = select(result, smallestMask, ones);
    return result;
}

func f32_4x
arm_sinf_4x(f32_4x y)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    f32_4x x = y;

    f32_4x piOver4;
    piOver4.m = _mm_set1_ps(gPiOver4);

    f32_4x minVal;
    minVal.m = _mm_set1_ps(0x1p-12f);

    f32_4x absTopY = abstop12_4x(y);
    f32_4x absTopPiOver4 = abstop12_4x(piOver4);
    f32_4x absMinVal = abstop12_4x(minVal);

    f32_4x smallestMask;
    smallestMask.m = _mm_cmplt_epi32(absTopY.mi, absMinVal.mi);
    f32_4x smallMask;
    smallMask.m = _mm_cmplt_epi32(absTopY.mi, absTopPiOver4.mi);

    f32_4x n4xMod;
    f32_4x xMod;
    xMod.m = reduce_fast_exp_4x(x.m, &n4xMod.mi);

    f32_4x xm;
    f32_4x n4x;
    xm.m = _mm_or_ps(_mm_andnot_ps(smallMask.m, xMod.m),
                     _mm_and_ps(smallMask.m, x.m));
    n4x.mi = _mm_andnot_si128(smallMask.mi, n4xMod.mi);

    f32_4x x2;
    x2.m = _mm_mul_ps(xm.m, xm.m);

    f32_4x sinV;
    f32_4x cosV;
    sinV.m = sinf_exp_poly_q0_4x(xm.m, x2.m);
    cosV.m = sinf_exp_poly_q1_4x(x2.m);

    f32_4x selection;
    f32_4x sel0Mask;
    f32_4x sel1Mask;
    f32_4x sel2Mask;
    f32_4x sel3Mask;
    selection.mi = _mm_and_si128(n4x.mi, _mm_set1_epi32(0x3));
    sel0Mask.mi = _mm_cmpeq_epi32(selection.mi, _mm_setzero_si128());
    sel1Mask.mi = _mm_cmpeq_epi32(selection.mi, _mm_set1_epi32(1));
    sel2Mask.mi = _mm_cmpeq_epi32(selection.mi, _mm_set1_epi32(2));
    sel3Mask.mi = _mm_cmpeq_epi32(selection.mi, _mm_set1_epi32(3));

    f32_4x result;
    result.m = _mm_and_ps(sel0Mask.m, sinV.m);
    result.m = _mm_or_ps(_mm_and_ps(sel1Mask.m, cosV.m), result.m);
    result.m = _mm_or_ps(_mm_and_ps(sel2Mask.m, (-sinV).m), result.m);
    result.m = _mm_or_ps(_mm_and_ps(sel3Mask.m, (-cosV).m), result.m);
    result.m = _mm_or_ps(_mm_and_ps(smallestMask.m, y.m),
                         _mm_andnot_ps(smallestMask.m, result.m));
    return result;
}

struct SinCos4x
{
    f32_4x cos;
    f32_4x sin;
};

func SinCos4x
arm_sincosf_4x(f32_4x y)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    f32_4x x = y;

    f32_4x ones;
    ones.m = _mm_set1_ps(1.0f);

    f32_4x piOver4;
    piOver4.m = _mm_set1_ps(gPiOver4);

    f32_4x minVal;
    minVal.m = _mm_set1_ps(0x1p-12f);

    f32_4x absTopY = abstop12_4x(y);
    f32_4x absTopPiOver4 = abstop12_4x(piOver4);
    f32_4x absMinVal = abstop12_4x(minVal);

    f32_4x smallestMask;
    smallestMask.m = _mm_cmplt_ps(absTopY.m, absMinVal.m);
    f32_4x smallMask;
    smallMask.m = _mm_cmplt_ps(absTopY.m, absTopPiOver4.m);

    f32_4x n4xMod;
    f32_4x xMod;
    xMod.m = reduce_fast_exp_4x(x.m, &n4xMod.mi);

    f32_4x xm;
    f32_4x n4x;
    xm.m = _mm_or_ps(_mm_andnot_ps(smallMask.m, xMod.m),
                     _mm_and_ps(smallMask.m, x.m));
    n4x.mi = _mm_andnot_si128(smallMask.mi, n4xMod.mi);

    f32_4x x2;
    x2.m = _mm_mul_ps(xm.m, xm.m);

    f32_4x sinV;
    f32_4x cosV;
    sinV.m = sinf_exp_poly_q0_4x(xm.m, x2.m);
    cosV.m = sinf_exp_poly_q1_4x(x2.m);

    f32_4x selection;
    f32_4x sel0Mask;
    f32_4x sel1Mask;
    f32_4x sel2Mask;
    f32_4x sel3Mask;
    selection.mi = _mm_and_si128(n4x.mi, _mm_set1_epi32(0x3));
    sel0Mask.mi = _mm_cmpeq_epi32(selection.mi, _mm_setzero_si128());
    sel1Mask.mi = _mm_cmpeq_epi32(selection.mi, _mm_set1_epi32(1));
    sel2Mask.mi = _mm_cmpeq_epi32(selection.mi, _mm_set1_epi32(2));
    sel3Mask.mi = _mm_cmpeq_epi32(selection.mi, _mm_set1_epi32(3));

    SinCos4x result;
    result.cos.m = _mm_and_ps(sel0Mask.m, cosV.m);
    result.cos.m = _mm_or_ps(_mm_and_ps(sel1Mask.m, (-sinV).m), result.cos.m);
    result.cos.m = _mm_or_ps(_mm_and_ps(sel2Mask.m, (-cosV).m), result.cos.m);
    result.cos.m = _mm_or_ps(_mm_and_ps(sel3Mask.m, sinV.m), result.cos.m);
    result.cos.m = _mm_or_ps(_mm_and_ps(smallestMask.m, ones.m),
                             _mm_andnot_ps(smallestMask.m, result.cos.m));
    result.sin.m = _mm_and_ps(sel0Mask.m, sinV.m);
    result.sin.m = _mm_or_ps(_mm_and_ps(sel1Mask.m, cosV.m), result.sin.m);
    result.sin.m = _mm_or_ps(_mm_and_ps(sel2Mask.m, (-sinV).m), result.sin.m);
    result.sin.m = _mm_or_ps(_mm_and_ps(sel3Mask.m, (-cosV).m), result.sin.m);
    result.sin.m = _mm_or_ps(_mm_and_ps(smallestMask.m, y.m),
                             _mm_andnot_ps(smallestMask.m, result.sin.m));
    return result;
}

#if 0
// NOTE(michiel): Above algorithms are based on these slimmed down versions
func f32
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

func f32
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
