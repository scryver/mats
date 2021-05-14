internal f64_2x
sinf_poly_q0_4x(f64_2x x, f64_2x x2)
{
    // NOTE(michiel): x - s1*x^3 + s2*x^5 - s3*x^7
    f64_2x x3 = x * x2;
    f64_2x s1 = F64_2x(0x1.1107605230bc4p-7) - x2 * F64_2x(0x1.994eb3774cf24p-13);
    f64_2x x7 = x3 * x2;
    f64_2x s = x - x3 * F64_2x(0x1.555545995a603p-3);
    return s + x7 * s1;
}

internal f64_2x
sinf_poly_q1_4x(f64_2x x2)
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
cos32_4x(f32_4x y)
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

    f64_2x sinLo = sinf_poly_q0_4x(xmLo, x2Lo);
    f64_2x sinHi = sinf_poly_q0_4x(xmHi, x2Hi);
    f64_2x cosLo = sinf_poly_q1_4x(x2Lo);
    f64_2x cosHi = sinf_poly_q1_4x(x2Hi);

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
sin32_4x(f32_4x y)
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

    f64_2x sinLo = sinf_poly_q0_4x(xmLo, x2Lo);
    f64_2x sinHi = sinf_poly_q0_4x(xmHi, x2Hi);
    f64_2x cosLo = sinf_poly_q1_4x(x2Lo);
    f64_2x cosHi = sinf_poly_q1_4x(x2Hi);

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

    f64_2x sinLo = sinf_poly_q0_4x(xmLo, x2Lo);
    f64_2x sinHi = sinf_poly_q0_4x(xmHi, x2Hi);
    f64_2x cosLo = sinf_poly_q1_4x(x2Lo);
    f64_2x cosHi = sinf_poly_q1_4x(x2Hi);

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

// NOTE(michiel): 4 calls to acos32/asin32/atan32 is faster at the moment

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

#if MATS_USE_SSE4
    f32_4x multIn; multIn.m = _mm_blendv_ps(half4x.m, x.m, lessThanHalf.m);
#else
    f32_4x multIn = select(half4x, lessThanHalf, x);
#endif
    f32_4x addIn  = and_not(S32_4x(0x80000000), s32_4x_or(lessThanHalf, lessThanZero));
    addIn = s32_4x_xor(addIn, x);
    f32_4x polyIn = and_not(one4x, lessThanHalf);
    polyIn = (polyIn + addIn) * multIn;;

#if MATS_ASINCOS_USE_SMALL_POLY
    f32_4x p = F32_4x(-8.6563630030e-03f) * polyIn;
    f32_4x q = F32_4x(-7.0662963390e-01f) * polyIn;
    p = (p - F32_4x(4.2743422091e-02f)) * polyIn;
    q = (q + F32_4x(1.0f));
    p = (p + F32_4x(1.6666586697e-01f)) * polyIn;
#else
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
#endif

    f32_4x xMinX = x - x;
#if MATS_USE_SSE4
    p.m = _mm_blendv_ps(p.m, xMinX.m, greaterThanOne.m);
    q.m = _mm_blendv_ps(q.m, xMinX.m, greaterThanOne.m);
#else
    p = select(p, greaterThanOne, xMinX);
    q = select(q, greaterThanOne, xMinX);
#endif

    // NOTE(michiel): Result for |x| > 1.0
    f32_4x r = p / q;

    // NOTE(michiel): Result for 0.5 <= |x| <= 1.0
    f32_4x s = sqrt32_4x(polyIn);
    f32_4x sMask = s32_4x_or(S32_4x(0xFFFFF000), s32_4x_and(S32_4x(0x00000FFF), lessThanZero));
    f32_4x m = s & sMask;
    f32_4x c = (polyIn - m * m) / (s + m);
#if MATS_USE_SSE4
    c.m = _mm_blendv_ps(c.m, F32_4x(-pio2_lo).m, lessThanZero.m);
#else
    c = select(c, lessThanZero, F32_4x(-pio2_lo));
#endif
    f32_4x w = r * s + c;
    f32_4x mult = m + w;
#if MATS_USE_SSE4
    mult.m = _mm_blendv_ps(mult.m, (lessThanZero & F32_4x(-pio2_lo)).m, isOne.m);
#else
    mult = select(mult, isOne, lessThanZero & F32_4x(-pio2_lo));
#endif
    f32_4x z = F32_4x(2.0f) * mult;

    // NOTE(michiel): Result for |x| < 0.5
    f32_4x smallResult = x * r;
    smallResult = smallResult & greaterThanMin;
    smallResult = F32_4x(pio2_lo) - smallResult;
    smallResult = x - smallResult;
    smallResult = F32_4x(pio2_hi) - smallResult;

#if MATS_USE_SSE4
    f32_4x normalResult; normalResult.m = _mm_blendv_ps(z.m, (F32_4x(pi) - z).m, lessThanZero.m);
    f32_4x result; result.m = _mm_blendv_ps(normalResult.m, smallResult.m, lessThanHalf.m);
    result.m = _mm_blendv_ps(result.m, r.m, greaterThanOne.m);
#else
    f32_4x normalResult = select(z, lessThanZero, F32_4x(pi) - z);
    f32_4x result = select(normalResult, lessThanHalf, smallResult);
    result = select(result, greaterThanOne, r);
#endif

    return result;
#undef pio2_hi
#undef pio2_lo
#undef pi
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

#if MATS_ASINCOS_USE_SMALL_POLY
    f32_4x p = F32_4x(-8.6563630030e-03f) * polyIn;
    f32_4x q = F32_4x(-7.0662963390e-01f) * polyIn;
    p = (p - F32_4x(4.2743422091e-02f)) * polyIn;
    q = (q + F32_4x(1.0f));
    p = (p + F32_4x(1.6666586697e-01f)) * polyIn;
#else
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
#endif

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
    result = select(result, greaterThanOne, r);
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

#if MATS_ATAN_USE_SMALL_POLY
    f32_4x s1 = F32_4x(6.1687607318e-02f) * x4;
    f32_4x s2 = F32_4x(-1.0648017377e-01f) * x4;
    s1 = (s1 + F32_4x(1.4253635705e-01f)) * x4;
    s2 = (s2 - F32_4x(1.9999158382e-01f)) * x4;
    s1 = (s1 + F32_4x(3.3333328366e-01f)) * x2;
#else
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
#endif

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
// TODO(michiel): When atan32_4x is faster we can look into widening this function
internal f32_4x
atan2_32_4x(f32_4x y, f32_4x x)
{

}
#endif
