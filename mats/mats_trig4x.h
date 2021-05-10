internal f32_4x
sinf_poly_q0_4x(f32_4x x, f32_4x x2)
{
    // NOTE(michiel): x - s1*x^3 + s2*x^5 - s3*x^7
    f32_4x x3 = x * x2;
    f32_4x s1 = F32_4x(0x1.110760p-7f) - x2 * F32_4x(0x1.994eb2p-13f);
    f32_4x x7 = x3 * x2;
    f32_4x s = x - x3 * F32_4x(0x1.555544p-3f);
    return s + x7 * s1;
}

internal f32_4x
sinf_poly_q1_4x(f32_4x x2)
{
    // NOTE(michiel): c0 - c1*x^2 + c2*x^4 - c3*x^6 + c4*x^8
    f32_4x x4 = x2 * x2;
    f32_4x c2 = F32_4x(-0x1.6c087ep-10f) + x2 * F32_4x(0x1.993430p-16f);
    f32_4x c1 = F32_4x(1.0f) - x2 * F32_4x(0x1.fffffep-2f);
    f32_4x x6 = x4 * x2;
    f32_4x c = c1 + x4 * F32_4x(0x1.55553ep-5f);
    return c + x6 * c2;
}

internal f64_2x
sinf_poly_q0_prec_4x(f64_2x x, f64_2x x2)
{
    // NOTE(michiel): x - s1*x^3 + s2*x^5 - s3*x^7
    f64_2x x3 = x * x2;
    f64_2x s1 = F64_2x(0x1.1107605230bc4p-7) - x2 * F64_2x(0x1.994eb3774cf24p-13);
    f64_2x x7 = x3 * x2;
    f64_2x s = x - x3 * F64_2x(0x1.555545995a603p-3);
    return s + x7 * s1;
}

internal f64_2x
sinf_poly_q1_prec_4x(f64_2x x2)
{
    // NOTE(michiel): c0 - c1*x^2 + c2*x^4 - c3*x^6 + c4*x^8
    f64_2x x4 = x2 * x2;
    f64_2x c2 = F64_2x(-0x1.6c087e89a359dp-10) + x2 * F64_2x(0x1.99343027bf8c3p-16);
    f64_2x c1 = F64_2x(1.0) - x2 * F64_2x(0x1.ffffffd0c621cp-2);
    f64_2x x6 = x4 * x2;
    f64_2x c = c1 + x4 * F64_2x(0x1.55553e1068f19p-5);
    return c + x6 * c2;
}

internal f32_4x
cos32_4x(f32_4x x)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    f32_4x ones = F32_4x(1.0f);
    f32_4x piOver4 = F32_4x(0.25f * F32_PI);
    f32_4x minVal = F32_4x(0x1p-12f);

    f32_4x absMask = S32_4x(0x7FF00000);
    f32_4x absTopX = s32_4x_and(x, absMask);
    f32_4x absTopPiOver4 = s32_4x_and(piOver4, absMask);
    f32_4x absTopMinVal = s32_4x_and(minVal, absMask);

    f32_4x smallestMask = s32_4x_less(absTopX, absTopMinVal);
    f32_4x smallMask = s32_4x_less(absTopX, absTopPiOver4);

    // NOTE(michiel): Reduce input argument
    f32_4x hpiInv = F32_4x(0x1.45F306p-1f);
    f32_4x r = x * hpiInv;
    f32_4x n4xMod = s32_4x_from_f32(r);
    f32_4x xMod   = x - round(r) * F32_4x(0x1.921FB4p0f);

#if MATS_USE_SSE4
    f32_4x xm; xm.m = _mm_blendv_ps(xMod.m, x.m, smallMask.m);
#else
    f32_4x xm = select(xMod, smallMask, x);
#endif
    f32_4x n4x = and_not(n4xMod, smallMask);

    f32_4x x2 = xm * xm;

    f32_4x sinV = sinf_poly_q0_4x(xm, x2);
    f32_4x cosV = sinf_poly_q1_4x(x2);

    f32_4x sel0 = s32_4x_equal(s32_4x_and(n4x, S32_4x(1)), S32_4x(1));
    f32_4x sel1 = s32_4x_equal(s32_4x_and(n4x, S32_4x(2)), S32_4x(2));

#if MATS_USE_SSE4
    f32_4x result; result.m = _mm_blendv_ps(cosV.m, sinV.m, sel0.m);
#else
    f32_4x result = select(cosV, sel0, sinV);
#endif

    f32_4x signMask = S32_4x(0x80000000);
    f32_4x selSign = s32_4x_xor(sel0, sel1);
    f32_4x signMod = s32_4x_and(selSign, signMask);
    result = result ^ signMod;
#if MATS_USE_SSE4
    result.m = _mm_blendv_ps(result.m, ones.m, smallestMask.m);
#else
    result = select(result, smallestMask, ones);
#endif

    return result;
}

internal f32_4x
sin32_4x(f32_4x x)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    f32_4x piOver4 = F32_4x(0.25f * F32_PI);
    f32_4x minVal = F32_4x(0x1p-12f);

    f32_4x absMask = S32_4x(0x7FF00000);
    f32_4x absTopX = s32_4x_and(x, absMask);
    f32_4x absTopPiOver4 = s32_4x_and(piOver4, absMask);
    f32_4x absTopMinVal = s32_4x_and(minVal, absMask);

    f32_4x smallestMask = s32_4x_less(absTopX, absTopMinVal);
    f32_4x smallMask = s32_4x_less(absTopX, absTopPiOver4);

    // NOTE(michiel): Reduce input argument
    f32_4x hpiInv = F32_4x(0x1.45F306p-1f);
    f32_4x r = x * hpiInv;
    f32_4x n4xMod = s32_4x_from_f32(r);
    f32_4x xMod   = x - round(r) * F32_4x(0x1.921FB4p0f);

#if MATS_USE_SSE4
    f32_4x xm; xm.m = _mm_blendv_ps(xMod.m, x.m, smallMask.m);
#else
    f32_4x xm = select(xMod, smallMask, x);
#endif
    f32_4x n4x = and_not(n4xMod, smallMask);

    f32_4x x2 = xm * xm;

    f32_4x sinV = sinf_poly_q0_4x(xm, x2);
    f32_4x cosV = sinf_poly_q1_4x(x2);

    f32_4x sel0 = s32_4x_equal(s32_4x_and(n4x, S32_4x(1)), S32_4x(1));
    f32_4x sel1 = s32_4x_equal(s32_4x_and(n4x, S32_4x(2)), S32_4x(2));

#if MATS_USE_SSE4
    f32_4x result; result.m = _mm_blendv_ps(sinV.m, cosV.m, sel0.m);
#else
    f32_4x result = select(sinV, sel0, cosV);
#endif

    f32_4x signMask = S32_4x(0x80000000);
    f32_4x signMod = s32_4x_and(sel1, signMask);
    result = result ^ signMod;
#if MATS_USE_SSE4
    result.m = _mm_blendv_ps(result.m, x.m, smallestMask.m);
#else
    result = select(result, smallestMask, x);
#endif

    return result;
}

internal f32_4x
cos32_prec_4x(f32_4x y)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    f32_4x ones = F32_4x(1.0f);
    f32_4x piOver4 = F32_4x(0.25f * F32_PI);
    f32_4x minVal = F32_4x(0x1p-12f);

    f32_4x absMask = S32_4x(0x7FF00000);
    f32_4x absTopY = s32_4x_and(y, absMask);
    f32_4x absTopPiOver4 = s32_4x_and(piOver4, absMask);
    f32_4x absTopMinVal = s32_4x_and(minVal, absMask);

    f32_4x smallestMask = s32_4x_less(absTopY, absTopMinVal);
    f32_4x smallMask = s32_4x_less(absTopY, absTopPiOver4);

    // NOTE(michiel): Reduce input argument
    f64_2x hpiInv = F64_2x(0x1.45F306DC9C883p-1);
    f64_2x xdLo = F64_2x(y);
    f64_2x xdHi = F64_2x(y, true);
    f64_2x rLo = xdLo * hpiInv;
    f64_2x rHi = xdHi * hpiInv;
    f32_4x n4xMod;
    n4xMod.mi = _mm_cvtpd_epi32(rLo.md);
    __m128 tempMod = _mm_castsi128_ps(_mm_cvtpd_epi32(rHi.md));
    n4xMod.m = _mm_shuffle_ps(n4xMod.m, tempMod, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
    f64_2x xModLo = xdLo - round(rLo) * F64_2x(0x1.921FB54442D18p0);
    f64_2x xModHi = xdHi - round(rHi) * F64_2x(0x1.921FB54442D18p0);

    f64_2x smallMaskLo, smallMaskHi;
    smallMaskLo.m = _mm_shuffle_ps(smallMask.m, smallMask.m, MULTILANE_SHUFFLE_MASK(0, 0, 1, 1));
    smallMaskHi.m = _mm_shuffle_ps(smallMask.m, smallMask.m, MULTILANE_SHUFFLE_MASK(2, 2, 3, 3));
#if MATS_USE_SSE4
    f64_2x xmLo, xmHi;
    xmLo.md = _mm_blendv_pd(xModLo.md, xdLo.md, smallMaskLo.md);
    xmHi.md = _mm_blendv_pd(xModHi.md, xdHi.md, smallMaskHi.md);
#else
    f64_2x xmLo = select(xModLo, smallMaskLo, xdLo);
    f64_2x xmHi = select(xModHi, smallMaskHi, xdHi);
#endif
    f32_4x n4x = and_not(n4xMod, smallMask);

    f64_2x x2Lo = xmLo * xmLo;
    f64_2x x2Hi = xmHi * xmHi;

    f64_2x sinLo = sinf_poly_q0_prec_4x(xmLo, x2Lo);
    f64_2x sinHi = sinf_poly_q0_prec_4x(xmHi, x2Hi);
    f64_2x cosLo = sinf_poly_q1_prec_4x(x2Lo);
    f64_2x cosHi = sinf_poly_q1_prec_4x(x2Hi);

    f32_4x sins = F32_4x(sinLo, sinHi);
    f32_4x coss = F32_4x(cosLo, cosHi);

    f32_4x sel0 = s32_4x_equal(s32_4x_and(n4x, S32_4x(1)), S32_4x(1));
    f32_4x sel1 = s32_4x_equal(s32_4x_and(n4x, S32_4x(2)), S32_4x(2));

#if MATS_USE_SSE4
    f32_4x result; result.m = _mm_blendv_ps(coss.m, sins.m, sel0.m);
#else
    f32_4x result = select(coss, sel0, sins);
#endif

    f32_4x signMask = S32_4x(0x80000000);
    f32_4x selSign = s32_4x_xor(sel0, sel1);
    f32_4x signMod = s32_4x_and(selSign, signMask);
    result = result ^ signMod;
#if MATS_USE_SSE4
    result.m = _mm_blendv_ps(result.m, ones.m, smallestMask.m);
#else
    result = select(result, smallestMask, ones);
#endif

    return result;
}

internal f32_4x
sin32_prec_4x(f32_4x y)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    f32_4x piOver4 = F32_4x(0.25f * F32_PI);
    f32_4x minVal = F32_4x(0x1p-12f);

    f32_4x absMask = S32_4x(0x7FF00000);
    f32_4x absTopY = s32_4x_and(y, absMask);
    f32_4x absTopPiOver4 = s32_4x_and(piOver4, absMask);
    f32_4x absTopMinVal = s32_4x_and(minVal, absMask);

    f32_4x smallestMask = s32_4x_less(absTopY, absTopMinVal);
    f32_4x smallMask = s32_4x_less(absTopY, absTopPiOver4);

    // NOTE(michiel): Reduce input argument
    f64_2x hpiInv = F64_2x(0x1.45F306DC9C883p-1);
    f64_2x xdLo = F64_2x(y);
    f64_2x xdHi = F64_2x(y, true);
    f64_2x rLo = xdLo * hpiInv;
    f64_2x rHi = xdHi * hpiInv;
    f32_4x n4xMod;
    n4xMod.mi = _mm_cvtpd_epi32(rLo.md);
    __m128 tempMod = _mm_castsi128_ps(_mm_cvtpd_epi32(rHi.md));
    n4xMod.m = _mm_shuffle_ps(n4xMod.m, tempMod, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
    f64_2x xModLo = xdLo - round(rLo) * F64_2x(0x1.921FB54442D18p0);
    f64_2x xModHi = xdHi - round(rHi) * F64_2x(0x1.921FB54442D18p0);

    f64_2x smallMaskLo, smallMaskHi;
    smallMaskLo.m = _mm_shuffle_ps(smallMask.m, smallMask.m, MULTILANE_SHUFFLE_MASK(0, 0, 1, 1));
    smallMaskHi.m = _mm_shuffle_ps(smallMask.m, smallMask.m, MULTILANE_SHUFFLE_MASK(2, 2, 3, 3));
#if MATS_USE_SSE4
    f64_2x xmLo, xmHi;
    xmLo.md = _mm_blendv_pd(xModLo.md, xdLo.md, smallMaskLo.md);
    xmHi.md = _mm_blendv_pd(xModHi.md, xdHi.md, smallMaskHi.md);
#else
    f64_2x xmLo = select(xModLo, smallMaskLo, xdLo);
    f64_2x xmHi = select(xModHi, smallMaskHi, xdHi);
#endif
    f32_4x n4x = and_not(n4xMod, smallMask);

    f64_2x x2Lo = xmLo * xmLo;
    f64_2x x2Hi = xmHi * xmHi;

    f64_2x sinLo = sinf_poly_q0_prec_4x(xmLo, x2Lo);
    f64_2x sinHi = sinf_poly_q0_prec_4x(xmHi, x2Hi);
    f64_2x cosLo = sinf_poly_q1_prec_4x(x2Lo);
    f64_2x cosHi = sinf_poly_q1_prec_4x(x2Hi);

    f32_4x sins = F32_4x(sinLo, sinHi);
    f32_4x coss = F32_4x(cosLo, cosHi);

    f32_4x sel0 = s32_4x_equal(s32_4x_and(n4x, S32_4x(1)), S32_4x(1));
    f32_4x sel1 = s32_4x_equal(s32_4x_and(n4x, S32_4x(2)), S32_4x(2));

#if MATS_USE_SSE4
    f32_4x result; result.m = _mm_blendv_ps(sins.m, coss.m, sel0.m);
#else
    f32_4x result = select(sins, sel0, coss);
#endif

    f32_4x signMask = S32_4x(0x80000000);
    f32_4x signMod = s32_4x_and(sel1, signMask);
    result = result ^ signMod;
#if MATS_USE_SSE4
    result.m = _mm_blendv_ps(result.m, y.m, smallestMask.m);
#else
    result = select(result, smallestMask, y);
#endif

    return result;
}

internal v2_4x
sincos32_4x(f32_4x y)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    f32_4x ones = F32_4x(1.0f);
    f32_4x piOver4 = F32_4x(0.25f * F32_PI);
    f32_4x minVal = F32_4x(0x1p-12f);

    f32_4x absMask = S32_4x(0x7FF00000);
    f32_4x absTopY = s32_4x_and(y, absMask);
    f32_4x absTopPiOver4 = s32_4x_and(piOver4, absMask);
    f32_4x absTopMinVal = s32_4x_and(minVal, absMask);

    f32_4x smallestMask = s32_4x_less(absTopY, absTopMinVal);
    f32_4x smallMask = s32_4x_less(absTopY, absTopPiOver4);

    // NOTE(michiel): Reduce input argument
    f32_4x hpiInv = F32_4x(0x1.45F306p-1f);
    f32_4x r = y * hpiInv;
    f32_4x n4xMod = s32_4x_from_f32(r);
    f32_4x xMod = y - round(r) * F32_4x(0x1.921FB4p0f);

#if MATS_USE_SSE4
    f32_4x xm; xm.m = _mm_blendv_ps(xMod.m, y.m, smallMask.m);
#else
    f32_4x xm = select(xMod, smallMask, y);
#endif
    f32_4x n4x = and_not(smallMask, n4xMod);

    f32_4x x2 = xm * xm;

    f32_4x sinV = sinf_poly_q0_4x(xm, x2);
    f32_4x cosV = sinf_poly_q1_4x(x2);

    f32_4x sel0 = s32_4x_equal(s32_4x_and(n4x, S32_4x(1)), S32_4x(1));
    f32_4x sel1 = s32_4x_equal(s32_4x_and(n4x, S32_4x(2)), S32_4x(2));

#if MATS_USE_SSE4
    f32_4x sins; sins.m = _mm_blendv_ps(sinV.m, cosV.m, sel0.m);
    f32_4x coss; coss.m = _mm_blendv_ps(cosV.m, sinV.m, sel0.m);
#else
    f32_4x sins = select(sinV, sel0, cosV);
    f32_4x coss = select(cosV, sel0, sinV);
#endif

    f32_4x signMask = S32_4x(0x80000000);
    f32_4x minSin = sins ^ signMask;
    f32_4x minCos = coss ^ signMask;

#if MATS_USE_SSE4
    sins.m = _mm_blendv_ps(sins.m, minSin.m, sel1.m);
    sins.m = _mm_blendv_ps(sins.m, y.m, smallestMask.m);
    coss.m = _mm_blendv_ps(coss.m, minCos.m, (sel0 ^ sel1).m);
    coss.m = _mm_blendv_ps(coss.m, ones.m, smallestMask.m);
#else
    sins = select(sins, sel1, minSin);
    sins = select(sins, smallestMask, y);
    coss = select(coss, sel0 ^ sel1, minCos);
    coss = select(coss, smallestMask, ones);
#endif

    v2_4x result;
    result.x = coss;
    result.y = sins;
    return result;
}

internal v2_4x
sincos32_prec_4x(f32_4x y)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    f32_4x ones = F32_4x(1.0f);
    f32_4x piOver4 = F32_4x(0.25f * F32_PI);
    f32_4x minVal = F32_4x(0x1p-12f);

    f32_4x absMask = S32_4x(0x7FF00000);
    f32_4x absTopY = s32_4x_and(y, absMask);
    f32_4x absTopPiOver4 = s32_4x_and(piOver4, absMask);
    f32_4x absTopMinVal = s32_4x_and(minVal, absMask);

    f32_4x smallestMask = s32_4x_less(absTopY, absTopMinVal);
    f32_4x smallMask = s32_4x_less(absTopY, absTopPiOver4);

    // NOTE(michiel): Reduce input argument
    f64_2x hpiInv = F64_2x(0x1.45F306DC9C883p-1);
    f64_2x xdLo = F64_2x(y);
    f64_2x xdHi = F64_2x(y, true);
    f64_2x rLo = xdLo * hpiInv;
    f64_2x rHi = xdHi * hpiInv;
    f32_4x n4xMod;
    n4xMod.mi = _mm_cvtpd_epi32(rLo.md);
    __m128 tempMod = _mm_castsi128_ps(_mm_cvtpd_epi32(rHi.md));
    n4xMod.m = _mm_shuffle_ps(n4xMod.m, tempMod, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
    f64_2x xModLo = xdLo - round(rLo) * F64_2x(0x1.921FB54442D18p0);
    f64_2x xModHi = xdHi - round(rHi) * F64_2x(0x1.921FB54442D18p0);

    f64_2x smallMaskLo, smallMaskHi;
    smallMaskLo.m = _mm_shuffle_ps(smallMask.m, smallMask.m, MULTILANE_SHUFFLE_MASK(0, 0, 1, 1));
    smallMaskHi.m = _mm_shuffle_ps(smallMask.m, smallMask.m, MULTILANE_SHUFFLE_MASK(2, 2, 3, 3));
#if MATS_USE_SSE4
    f64_2x xmLo, xmHi;
    xmLo.md = _mm_blendv_pd(xModLo.md, xdLo.md, smallMaskLo.md);
    xmHi.md = _mm_blendv_pd(xModHi.md, xdHi.md, smallMaskHi.md);
#else
    f64_2x xmLo = select(xModLo, smallMaskLo, xdLo);
    f64_2x xmHi = select(xModHi, smallMaskHi, xdHi);
#endif
    f32_4x n4x = and_not(n4xMod, smallMask);

    f64_2x x2Lo = xmLo * xmLo;
    f64_2x x2Hi = xmHi * xmHi;

    f64_2x sinLo = sinf_poly_q0_prec_4x(xmLo, x2Lo);
    f64_2x sinHi = sinf_poly_q0_prec_4x(xmHi, x2Hi);
    f64_2x cosLo = sinf_poly_q1_prec_4x(x2Lo);
    f64_2x cosHi = sinf_poly_q1_prec_4x(x2Hi);

    f32_4x sins = F32_4x(sinLo, sinHi);
    f32_4x coss = F32_4x(cosLo, cosHi);

    f32_4x sel0 = s32_4x_equal(s32_4x_and(n4x, S32_4x(1)), S32_4x(1));
    f32_4x sel1 = s32_4x_equal(s32_4x_and(n4x, S32_4x(2)), S32_4x(2));

    f32_4x cosResult;
    f32_4x sinResult;
#if MATS_USE_SSE4
    cosResult.m = _mm_blendv_ps(coss.m, sins.m, sel0.m);
    sinResult.m = _mm_blendv_ps(sins.m, coss.m, sel0.m);
#else
    cosResult = select(coss, sel0, sins);
    sinResult = select(sins, sel0, coss);
#endif

    f32_4x signMask = S32_4x(0x80000000);
    f32_4x selCosSign = s32_4x_xor(sel0, sel1);
    f32_4x selSinSign = sel1;
    f32_4x minCos = cosResult ^ signMask;
    f32_4x minSin = sinResult ^ signMask;
#if MATS_USE_SSE4
    cosResult.m = _mm_blendv_ps(cosResult.m, minCos.m, selCosSign.m);
    cosResult.m = _mm_blendv_ps(cosResult.m, ones.m, smallestMask.m);
    sinResult.m = _mm_blendv_ps(sinResult.m, minSin.m, selSinSign.m);
    sinResult.m = _mm_blendv_ps(sinResult.m, y.m, smallestMask.m);
#else
    cosResult = select(cosResult, selCosSign, minCos);
    cosResult = select(cosResult, smallestMask, ones);
    sinResult = select(sinResult, selSinSign, minSin);
    sinResult = select(sinResult, smallestMask, y);
#endif

    v2_4x result;
    result.x = cosResult;
    result.y = sinResult;
    return result;
}

internal f32_4x
tan32_kernel_4x(f32_4x x, f32_4x mod)
{
    f32_4x xOrig = x;
    f32_4x hx = x & S32_4x(0x7FFFFFFF);

    f32_4x specialMask = s32_4x_less(hx, S32_4x(0x31800000));
    specialMask = s32_4x_and(s32_4x_equal(s32_4x_from_f32_trunc(x), zero_f32_4x()), specialMask);
    f32_4x special = zero_f32_4x();
    {
        f32_4x hxModP1Mask = s32_4x_equal(s32_4x_or(hx, s32_4x_add(mod, S32_4x(1))), zero_f32_4x());
        hxModP1Mask = s32_4x_and(hxModP1Mask, specialMask);
        special = (F32_4x(1.0f) / absolute(x)) & hxModP1Mask;
        f32_4x workMask = and_not(specialMask, hxModP1Mask);
        f32_4x modIsOneMask = s32_4x_equal(mod, S32_4x(1));
#if MATS_USE_SSE4
        f32_4x modIsOne; modIsOne.m = _mm_blendv_ps((F32_4x(-1.0f) / x).m, x.m, modIsOneMask.m);
#else
        f32_4x modIsOne = select(F32_4x(-1.0f) / x, modIsOneMask, x);
#endif
        special = special | (modIsOne & workMask);
    }

    f32_4x hxGreat = s32_4x_less(S32_4x(0x3F2CA140), hx);
    f32_4x xLessZero = s32_4x_less(x, zero_f32_4x());
    xLessZero = s32_4x_and(hxGreat, s32_4x_and(xLessZero, S32_4x(0x80000000)));
    x = x ^ xLessZero;
    f32_4x zPre = F32_4x(gPiOver4F32_hi) - x;
    f32_4x wPre = F32_4x(gPiOver4F32_lo) - (F32_4x(gPiOver4F32_hi) - zPre - x);
    f32_4x xPre = zPre + wPre;
#if MATS_USE_SSE4
    x.m = _mm_blendv_ps(x.m, xPre.m, hxGreat.m);
#else
    x = select(x, hxGreat, xPre);
#endif
    f32_4x x2 = x * x;
    f32_4x x4 = x2 * x2;

    f32_4x r = F32_4x(-1.8558637748e-05f) * x4;
    f32_4x v = F32_4x(2.5907305826e-05f) * x4;
    r = (r + F32_4x(7.8179444245e-05f)) * x4;
    v = (v + F32_4x(7.1407252108e-05f)) * x4;
    r = (r + F32_4x(5.8804126456e-04f)) * x4;
    v = (v + F32_4x(2.4646313977e-04f)) * x4;
    r = (r + F32_4x(3.5920790397e-03f)) * x4;
    v = (v + F32_4x(1.4562094584e-03f)) * x4;
    r = (r + F32_4x(2.1869488060e-02f)) * x4;
    v = (v + F32_4x(8.8632395491e-03f)) * x4;
    r = (r + F32_4x(1.3333334029e-01f));
    v = (v + F32_4x(5.3968254477e-02f)) * x2;

    f32_4x x3 = x2 * x;
    r = x2 * (x3 * (r + v));
    r = r + F32_4x(3.3333334327e-01f) * x3;

    f32_4x result = x + r;

    f32_4x oneOrMinus = F32_4x(1.0f) ^ (xOrig & F32_4x(0x80000000));
    v = f32_4x_from_s32(mod);
    f32_4x greatRes = oneOrMinus * (v - F32_4x(2.0f) * (x - (result * result / (result + v) - r)));

    f32_4x nonModdedMask = s32_4x_equal(mod, S32_4x(1));
    f32_4x modded = result & S32_4x(0xFFFFF000);
    v = r - (modded - x);
    f32_4x a = F32_4x(-1.0f) / result;
    f32_4x t = a & S32_4x(0xFFFFF000);
    f32_4x s = F32_4x(1.0f) + t * modded;
    modded = t + a * (s + t * v);

#if MATS_USE_SSE4
    result.m = _mm_blendv_ps(modded.m, result.m, nonModdedMask.m);
    result.m = _mm_blendv_ps(result.m, greatRes.m, hxGreat.m);
    result.m = _mm_blendv_ps(result.m, special.m, specialMask.m);
#else
    result = select(modded, nonModdedMask, result);
    result = select(result, hxGreat, greatRes);
    result = select(result, specialMask, special);
#endif

    return result;
}

internal f32_4x
tan32_4x(f32_4x x)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    f32_4x absX = x & S32_4x(0x7FFFFFFF);
    f32_4x inputSmall = s32_4x_less(absX, S32_4x(0x3F490FDA));

    // NOTE(michiel): Reduce input argument
    f64_2x hpiInv = F64_2x(0x1.45F306DC9C883p-1);
    f64_2x xdLo = F64_2x(x);
    f64_2x xdHi = F64_2x(x, true);
    f64_2x rLo = xdLo * hpiInv;
    f64_2x rHi = xdHi * hpiInv;
    f32_4x n4xMod;
    n4xMod.mi = _mm_cvtpd_epi32(rLo.md);
    __m128 tempMod = _mm_castsi128_ps(_mm_cvtpd_epi32(rHi.md));
    n4xMod.m = _mm_shuffle_ps(n4xMod.m, tempMod, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
    f64_2x xModLo = xdLo - round(rLo) * F64_2x(0x1.921FB54442D18p0);
    f64_2x xModHi = xdHi - round(rHi) * F64_2x(0x1.921FB54442D18p0);
    f32_4x xMod = F32_4x(xModLo, xModHi);

    f32_4x nMod = s32_4x_sub(S32_4x(1), s32_4x_sll(s32_4x_and(n4xMod, S32_4x(1)), 1));;
#if MATS_USE_SSE4
    f32_4x xIn; xIn.m = _mm_blendv_ps(xMod.m, x.m, inputSmall.m);
    f32_4x nIn; nIn.m = _mm_blendv_ps(nMod.m, _mm_castsi128_ps(_mm_set1_epi32(1)), inputSmall.m);
#else
    f32_4x xIn = select(xMod, inputSmall, x);
    f32_4x nIn = select(nMod, inputSmall, S32_4x(1));
#endif

    f32_4x result = tan32_kernel_4x(xIn, nIn);

    return result;
}

internal f32
acos32_temp(f32 x)
{
#define pi        3.1415925026e+00f
#define pio2_hi   1.5707962513e+00f
#define pio2_lo   7.5497894159e-08f

	s32 ix = (s32)u32f32(x).u;
    u32 hx = ix & 0x7FFFFFFF;

    b32 greaterThanOne = hx > 0x3F800000;
    b32 lessThanHalf = hx < 0x3F000000;
    b32 lessThanZero = ix < 0;
    b32 isOne = hx == 0x3F800000;
    b32 greaterThanMin = hx > 0x23000000;

    f32 polyIn = 0.0f;
    f32 p = 0.0f;;
    f32 q = 1.0f;
    if (greaterThanOne) {
        p = x - x;
        q = x - x;
    } else {
        f32 multIn = lessThanHalf ? x : 0.5f;
        f32 addIn  = (lessThanHalf || lessThanZero) ? x : -x;
        polyIn = lessThanHalf ? 0.0f : 1.0f;
        polyIn = (polyIn + addIn) * multIn;

#if MATS_USE_SSE
        WideMath polyW; polyW.m = _mm_setr_ps(polyIn, polyIn, 0.0f, 0.0f);
        WideMath polyE; polyE.m = _mm_setr_ps(polyIn, 1.0f, 0.0f, 0.0f);

        WideMath coef0; coef0.m = _mm_setr_ps( 1.6666667163e-01f,  0.0f, 0.0f, 0.0f);
        WideMath coef1; coef1.m = _mm_setr_ps(-3.2556581497e-01f,  1.0f, 0, 0);
        WideMath coef2; coef2.m = _mm_setr_ps( 2.0121252537e-01f, -2.4033949375e+00f, 0, 0);
        WideMath coef3; coef3.m = _mm_setr_ps(-4.0055535734e-02f,  2.0209457874e+00f, 0, 0);
        WideMath coef4; coef4.m = _mm_setr_ps( 7.9153501429e-04f, -6.8828397989e-01f, 0, 0);
        WideMath coef5; coef5.m = _mm_setr_ps( 3.4793309169e-05f,  7.7038154006e-02f, 0, 0);
        WideMath s1s2; s1s2.m = _mm_mul_ps(coef5.m, polyW.m);
        s1s2.m = _mm_mul_ps(_mm_add_ps(s1s2.m, coef4.m), polyW.m);
        s1s2.m = _mm_mul_ps(_mm_add_ps(s1s2.m, coef3.m), polyW.m);
        s1s2.m = _mm_mul_ps(_mm_add_ps(s1s2.m, coef2.m), polyW.m);
        s1s2.m = _mm_mul_ps(_mm_add_ps(s1s2.m, coef1.m), polyE.m);
        s1s2.m = _mm_mul_ps(_mm_add_ps(s1s2.m, coef0.m), polyE.m);

        p = s1s2.e[0];
        q = s1s2.e[1];
#else
        p = 3.4793309169e-05f * polyIn;
        q = 7.7038154006e-02f * polyIn;
        p = (p + 7.9153501429e-04f) * polyIn;
        q = (q - 6.8828397989e-01f) * polyIn;
        p = (p - 4.0055535734e-02f) * polyIn;
        q = (q + 2.0209457874e+00f) * polyIn;
        p = (p + 2.0121252537e-01f) * polyIn;
        q = (q - 2.4033949375e+00f) * polyIn;
        p = (p - 3.2556581497e-01f) * polyIn;
        q += 1.0f;
        p = (p + 1.6666667163e-01f) * polyIn;
#endif
    }

    f32 result = p / q;

    if (lessThanHalf)
    {
        f32 a = greaterThanMin ? x * result : 0.0f;
        f32 b = pio2_lo - a;
        f32 c = x - b;
        f32 d = pio2_hi - c;
        result = d;
    }
    else if (!greaterThanOne)
    {
        f32 s = sqrt32(polyIn);
        f32 c;
        f32 m;
        if (ix < 0) {
            m = s;
            c = -pio2_lo;
        } else {
            m = u32f32(u32f32(s).u & 0xFFFFF000).f;
            c = (polyIn - m * m) / (s + m);
        }
        f32 w = result * s + c;
        f32 mult = m + w;
        mult = isOne ? (lessThanZero ? -pio2_lo : 0.0f) : mult;
        f32 z = 2.0f * mult;
        result = ix < 0 ? pi - z : z;
    }

    return result;
#undef pio2_hi
#undef pio2_lo
#undef pi
}

internal f32_4x
acos32_4x(f32_4x x)
{
#define pi        3.1415925026e+00f
#define pio2_hi   1.5707962513e+00f
#define pio2_lo   7.5497894159e-08f
    f32_4x zero4x = F32_4x(0.0f);
    f32_4x one4x  = F32_4x(1.0f);
    f32_4x half4x = F32_4x(0.5f);

    f32_4x hx = x & F32_4x(0x7FFFFFFFU);

    f32_4x greaterThanOne = s32_4x_greater(hx, one4x);
    f32_4x lessThanHalf   = s32_4x_less(hx, half4x);
    f32_4x lessThanZero   = s32_4x_less(x, zero4x);
    f32_4x isOne          = s32_4x_equal(hx, one4x);
    f32_4x greaterThanMin = s32_4x_greater(hx, S32_4x(0x23000000));

    f32_4x multIn = select(half4x, lessThanHalf, x);
    f32_4x addIn  = and_not(S32_4x(0x80000000), s32_4x_or(lessThanHalf, lessThanZero));
    addIn = s32_4x_xor(addIn, x);
    f32_4x polyIn = and_not(one4x, lessThanHalf);
    polyIn = (polyIn + addIn) * multIn;;

    f32_4x p = F32_4x(3.4793309169e-05f) * polyIn;
    f32_4x q = F32_4x(7.7038154006e-02f) * polyIn;
    p = (p + F32_4x(7.9153501429e-04f)) * polyIn;
    q = (q - F32_4x(6.8828397989e-01f)) * polyIn;
    p = (p - F32_4x(4.0055535734e-02f)) * polyIn;
    q = (q + F32_4x(2.0209457874e+00f)) * polyIn;
    p = (p + F32_4x(2.0121252537e-01f)) * polyIn;
    q = (q - F32_4x(2.4033949375e+00f)) * polyIn;
    p = (p - F32_4x(3.2556581497e-01f)) * polyIn;
    q = (q + F32_4x(1.0f));
    p = (p + F32_4x(1.6666667163e-01f)) * polyIn;

    f32_4x xMinX = x - x;
    p = select(p, greaterThanOne, xMinX);
    q = select(q, greaterThanOne, xMinX);

    // NOTE(michiel): Result for |x| > 1.0
    f32_4x r = p / q;

    // NOTE(michiel): Result for |x| < 0.5
    f32_4x smallResult = x * r;
    smallResult = smallResult & greaterThanMin;
    smallResult = F32_4x(pio2_lo) - smallResult;
    smallResult = x - smallResult;
    smallResult = F32_4x(pio2_hi) - smallResult;

    // NOTE(michiel): Result for 0.5 <= |x| <= 1.0
    f32_4x s = sqrt32_4x(polyIn);
    f32_4x sMask = s32_4x_or(S32_4x(0xFFFFF000), s32_4x_and(S32_4x(0x00000FFF), lessThanZero));
    f32_4x m = s & sMask;
    f32_4x c = (polyIn - m * m) / (s + m);
    c = select(c, lessThanZero, F32_4x(-pio2_lo));
    f32_4x w = r * s + c;
    f32_4x mult = m + w;
    mult = select(mult, isOne, lessThanZero & F32_4x(-pio2_lo));
    f32_4x z = F32_4x(2.0f) * mult;

    f32_4x normalResult = select(z, lessThanZero, F32_4x(pi) - z);

    f32_4x result = select(normalResult, lessThanHalf, smallResult);
    result = select(result, greaterThanOne, r);

    return result;
#undef pio2_hi
#undef pio2_lo
#undef pi
}

internal f32
asin32_temp(f32 x)
{
#define pio2_hi 1.57079637050628662109375f
#define pio2_lo -4.37113900018624283e-8f
#define pio4_hi 0.785398185253143310546875f

    s32 ix = (s32)u32f32(x).u;
    u32 hx = ix & 0x7FFFFFFF;

    b32 greaterThanOne  = hx >  0x3F800000;
    b32 lessThanHalf    = hx <  0x3F000000;
    b32 greaterThanZero = ix >  0;
    b32 isOne           = hx == 0x3F800000;
    b32 lessThanMin     = hx <  0x32000000;
    b32 lessThanMax     = hx <  0x3F79999A;

    f32 polyIn = 0.0f;
    f32 p = 0.0f;
    f32 q = 1.0f;
    if (greaterThanOne)
    {
        p = x - x;
        q = x - x;
    }
    else
    {
        polyIn = lessThanHalf ? x * x : (1.0f - u32f32(hx).f) * 0.5f;

        p = 3.4793309169e-05f * polyIn;
        p = (p + 7.9153501429e-04f) * polyIn;
        p = (p - 4.0055535734e-02f) * polyIn;
        p = (p + 2.0121252537e-01f) * polyIn;
        p = (p - 3.2556581497e-01f) * polyIn;
        p = (p + 1.6666667163e-01f) * polyIn;
        q = 7.7038154006e-02f * polyIn;
        q = (q - 6.8828397989e-01f) * polyIn;
        q = (q + 2.0209457874e+00f) * polyIn;
        q = (q - 2.4033949375e+00f) * polyIn;
        q = q + 1.0f;
    }

    f32 s = sqrt32(polyIn);
    f32 r = p / q;

    f32 result = polyIn;
    if (isOne)
    {
        return x * pio2_hi + x * pio2_lo;
    }
    else if (lessThanHalf) {
        if (lessThanMin) {
            return x;
        } else {
            return x + x * r;
        }
    }
    else if (lessThanMax)
    {
        f32 w = s;
        u32 iw = u32f32(w).u;
        w = u32f32(iw & 0xFFFFF000).f;
        f32 c = (result - w * w) / (s + w);
        p = 2.0f * s * r - (pio2_lo - 2.0f * c);
        q = pio4_hi - 2.0f * w;
        result = pio4_hi - (p - q);
    }
    else
    {
        result = pio2_hi - (2.0f * (s + s * r) - pio2_lo);
    }

    return (greaterThanZero ? result : -result);

#undef pio4_hi
#undef pio2_hi
#undef pio2_lo
}

internal f32_4x
asin32_4x(f32_4x x)
{
    f32_4x zero4x = F32_4x(0.0f);
    f32_4x one4x  = F32_4x(1.0f);
    f32_4x two4x  = F32_4x(2.0f);
    f32_4x half4x = F32_4x(0.5f);
    f32_4x piOver2Hi = F32_4x(1.57079637050628662109375f);
    f32_4x piOver2Lo = F32_4x(-4.37113900018624283e-8f);
    f32_4x piOver4Hi = F32_4x(0.785398185253143310546875f);

    f32_4x hx = x & F32_4x(0x7FFFFFFFU);

    f32_4x greaterThanOne  = s32_4x_greater(hx, one4x);
    f32_4x lessThanHalf    = s32_4x_less(hx, half4x);
    f32_4x greaterThanZero = s32_4x_greater(x, zero4x);
    f32_4x isOne           = s32_4x_equal(hx, one4x);
    f32_4x lessThanMin     = s32_4x_less(hx, S32_4x(0x32000000));
    f32_4x lessThanMax     = s32_4x_less(hx, S32_4x(0x3F79999A));

    f32_4x polyIn = select(half4x * (one4x - hx), lessThanHalf, x * x);

    f32_4x p = F32_4x(3.4793309169e-05f) * polyIn;
    p = (p + F32_4x(7.9153501429e-04f)) * polyIn;
    p = (p - F32_4x(4.0055535734e-02f)) * polyIn;
    p = (p + F32_4x(2.0121252537e-01f)) * polyIn;
    p = (p - F32_4x(3.2556581497e-01f)) * polyIn;
    p = (p + F32_4x(1.6666667163e-01f)) * polyIn;
    f32_4x q = F32_4x(7.7038154006e-02f) * polyIn;
    q = (q - F32_4x(6.8828397989e-01f)) * polyIn;
    q = (q + F32_4x(2.0209457874e+00f)) * polyIn;
    q = (q - F32_4x(2.4033949375e+00f)) * polyIn;
    q = (q + F32_4x(1.0f));

    f32_4x xMinX = x - x;
    p = select(p, greaterThanOne, xMinX);
    q = select(q, greaterThanOne, xMinX);

    f32_4x r = p / q;

    // NOTE(michiel): Result for |x| == 1.0
    f32_4x oneResult = x * piOver2Hi + x * piOver2Lo;

    // NOTE(michiel): Result for |x| < 0.5
    f32_4x smallResult = x * r;
    smallResult = x + smallResult;
    smallResult = select(smallResult, lessThanMin, x);

    // NOTE(michiel): Result for 0.5 <= |x| < 0.9999 something
    f32_4x s = sqrt32_4x(polyIn);
    f32_4x w = s & S32_4x(0xFFFFF000);
    f32_4x c = (polyIn - w * w) / (s + w);
    f32_4x p2 = two4x * s * r - (piOver2Lo - two4x * c);
    f32_4x q2 = piOver4Hi - two4x * w;
    f32_4x normalResult = piOver4Hi - (p2 - q2);

    f32_4x result = piOver2Hi - (two4x * (s + s * r) - piOver2Lo);
    result = select(result, lessThanMax, normalResult);
    result = result ^ and_not(S32_4x(0x80000000), greaterThanZero);
    result = select(result, lessThanHalf, smallResult);
    result = select(result, isOne, oneResult);
    return result;
}

internal f32
atan32_temp(f32 x)
{
    s32 ix = (s32)u32f32(x).u;
    u32 hx = ix & 0x7FFFFFFF;

    b32 lessThanZero = ix < 0;

    f32 result;
    if (hx >= 0x50800000)
    {
        // NOTE(michiel): |x| >= 2^34
        f32 a = 1.5707962513e+00f;
        f32 b = 7.5497894159e-08f;
        if (FLT_UWORD_IS_NAN(hx)) {
            a = x;
            b = x;
        } else if (lessThanZero) {
            a = -a;
            b = -b;
        }
        result = a + b;
    }
    else
    {
        b32 smallMask    = hx < 0x3EE00000; // NOTE(michiel): |x| < 0.4375
        b32 smallestMask = hx < 0x31000000; // NOTE(michiel): |x| < 2^-29
        b32 middleMask   = hx < 0x3F980000; // NOTE(michiel): |x| < 1.1875
        b32 middlestMask = hx < 0x3F300000; // NOTE(michiel): 7/16 <= |x| < 11/16
        b32 notBigMask   = hx < 0x401C0000;
        u32 id = 0;
        if (!smallMask)
        {
            x = u32f32(hx).f;
            f32 a = 0.0f;
            f32 b = 1.0f;
            if (middlestMask) {
                a = 2.0f;
                id = 0;
            } else if (middleMask) {
                a = 1.0f;
                id = 1;
            } else if (notBigMask) {
                a = 1.0f;
                b = 1.5f;
                id = 2;
            } else {
                id = 3;
            }
            f32 nom = a * x - b;
            f32 den = a + b * x;
            x = nom / den;
        }

        f32 x2 = x * x;
        f32 x4 = x2 * x2;

        f32 s1 =   1.6285819933e-02f * x4;
        f32 s2 = - 3.6531571299e-02f * x4;
        s1 = (s1 + 4.9768779427e-02f) * x4;
        s2 = (s2 - 5.8335702866e-02f) * x4;
        s1 = (s1 + 6.6610731184e-02f) * x4;
        s2 = (s2 - 7.6918758452e-02f) * x4;
        s1 = (s1 + 9.0908870101e-02f) * x4;
        s2 = (s2 - 1.1111110449e-01f) * x4;
        s1 = (s1 + 1.4285714924e-01f) * x4;
        s2 = (s2 - 2.0000000298e-01f) * x4;
        s1 = (s1 + 3.3333334327e-01f) * x2;

        result = x * (s1 + s2);
        if (smallMask) {
            result = smallestMask ? x : x - result;
        } else {
            result = result - gAtanLoF32[id];
            result = result - x;
            result = gAtanHiF32[id] - result;
            result = lessThanZero ? -result : result;
        }
    }

    return result;
}

alignas(16) global const f32_4x gAtanHiF32_4x = {
    4.6364760399e-01f, 7.8539812565e-01f, 9.8279368877e-01f, 1.5707962513e+00f, /* atan(inf)hi 0x3fc90fda */
};

alignas(16) global const f32_4x gAtanLoF32_4x = {
    5.0121582440e-09f, 3.7748947079e-08f, 3.4473217170e-08f, 7.5497894159e-08f, /* atan(inf)lo 0x33a22168 */
};

internal f32_4x
atan32_4x(f32_4x x)
{
    f32_4x hx = x & F32_4x(0x7FFFFFFFU);

    f32_4x lessThanZero   = s32_4x_less(x, zero_f32_4x());
    f32_4x lessThanMax    = s32_4x_less(hx, S32_4x(0x50800000));
    f32_4x isNan          = s32_4x_greater(hx, S32_4x(0x7f800000));

    // NOTE(michiel): x >= 0x50800000
    f32_4x aBig = F32_4x(1.5707962513e+00f);
    f32_4x bBig = F32_4x(7.5497894159e-08f);
    aBig = select(aBig, lessThanZero, -aBig);
    bBig = select(bBig, lessThanZero, -bBig);
    aBig = select(aBig, isNan, x);
    bBig = select(bBig, isNan, x);
    f32_4x bigResult = aBig + bBig;

    // NOTE(michiel): Small result
    f32_4x smallMask    = s32_4x_less(hx, S32_4x(0x3EE00000));
    f32_4x smallestMask = s32_4x_less(hx, S32_4x(0x31000000));
    f32_4x middleMask   = s32_4x_less(hx, S32_4x(0x3F980000));
    f32_4x middlestMask = s32_4x_less(hx, S32_4x(0x3F300000));
    f32_4x notBigMask   = s32_4x_less(hx, S32_4x(0x401C0000));

    f32_4x ids = S32_4x(0x0F0E0D0C);
    f32_4x aSmall = zero_f32_4x();
    f32_4x bSmall = F32_4x(1.0f);
    aSmall = select(aSmall, notBigMask, F32_4x(1.0f));
    bSmall = select(bSmall, notBigMask, F32_4x(1.5f));
    ids    = select(ids, notBigMask, S32_4x(0x0B0A0908));
    aSmall = select(aSmall, middleMask, F32_4x(1.0f));
    bSmall = select(bSmall, middleMask, F32_4x(1.0f));
    ids    = select(ids, middleMask, S32_4x(0x07060504));
    aSmall = select(aSmall, middlestMask, F32_4x(2.0f));
    bSmall = select(bSmall, middlestMask, F32_4x(1.0f));
    ids    = select(ids, middlestMask, S32_4x(0x03020100));
    f32_4x nom = aSmall * hx - bSmall;
    f32_4x den = aSmall + bSmall * hx;
    f32_4x xMod = nom / den;

    x = select(xMod, smallMask, x);

    f32_4x x2 = x * x;
    f32_4x x4 = x2 * x2;

    f32_4x s1 = F32_4x(1.6285819933e-02f) * x4;
    f32_4x s2 = F32_4x(-3.6531571299e-02f) * x4;
    s1 = (s1 + F32_4x(4.9768779427e-02f)) * x4;
    s2 = (s2 - F32_4x(5.8335702866e-02f)) * x4;
    s1 = (s1 + F32_4x(6.6610731184e-02f)) * x4;
    s2 = (s2 - F32_4x(7.6918758452e-02f)) * x4;
    s1 = (s1 + F32_4x(9.0908870101e-02f)) * x4;
    s2 = (s2 - F32_4x(1.1111110449e-01f)) * x4;
    s1 = (s1 + F32_4x(1.4285714924e-01f)) * x4;
    s2 = (s2 - F32_4x(2.0000000298e-01f)) * x4;
    s1 = (s1 + F32_4x(3.3333334327e-01f)) * x2;

    f32_4x smallResult = x * (s1 + s2);

    f32_4x atanLo = byte_shuffle(gAtanLoF32_4x, ids);
    f32_4x atanHi = byte_shuffle(gAtanHiF32_4x, ids);
    f32_4x result = smallResult - atanLo;
    result = result - x;
    result = atanHi - result;
    result = select(result, lessThanZero, -result);

    result = select(result, smallMask, x - smallResult);
    result = select(result, smallestMask, x);

    result = select(bigResult, lessThanMax, result);

    return result;
}

#if 0
internal f32_4x
atan2_32_4x(f32_4x y, f32_4x x)
{

}
#endif
