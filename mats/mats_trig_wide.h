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

    f64_2x xmLo = select2x(xModLo, smallMaskLo, xdLo);
    f64_2x xmHi = select2x(xModHi, smallMaskHi, xdHi);
    f32_4x n4x = and_not(n4xMod, smallMask);

    f32_4x sel0 = s32_4x_equal(s32_4x_and(n4x, S32_4x(1)), S32_4x(1));
    f32_4x sel1 = s32_4x_equal(s32_4x_and(n4x, S32_4x(2)), S32_4x(2));

    f64_2x x2Lo = xmLo * xmLo;
    f64_2x x2Hi = xmHi * xmHi;

    f64_2x sinLo = sinf_poly_q0_4x(xmLo, x2Lo);
    f64_2x sinHi = sinf_poly_q0_4x(xmHi, x2Hi);
    f64_2x cosLo = sinf_poly_q1_4x(x2Lo);
    f64_2x cosHi = sinf_poly_q1_4x(x2Hi);

    f32_4x sins = F32_4x(sinLo, sinHi);
    f32_4x coss = F32_4x(cosLo, cosHi);

    f32_4x result = select4x(coss, sel0, sins);

    f32_4x signMask = S32_4x(MATS_F32_SIGN_MASK);
    f32_4x selSign = s32_4x_xor(sel0, sel1);
    f32_4x signMod = s32_4x_and(selSign, signMask);
    result = result ^ signMod;
    result = select4x(result, smallestMask, ones);

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

    f64_2x xmLo = select2x(xModLo, smallMaskLo, xdLo);
    f64_2x xmHi = select2x(xModHi, smallMaskHi, xdHi);
    f32_4x n4x = and_not(n4xMod, smallMask);

    f32_4x sel0 = s32_4x_equal(s32_4x_and(n4x, S32_4x(1)), S32_4x(1));
    f32_4x sel1 = s32_4x_equal(s32_4x_and(n4x, S32_4x(2)), S32_4x(2));

    f64_2x x2Lo = xmLo * xmLo;
    f64_2x x2Hi = xmHi * xmHi;

    f64_2x sinLo = sinf_poly_q0_4x(xmLo, x2Lo);
    f64_2x sinHi = sinf_poly_q0_4x(xmHi, x2Hi);
    f64_2x cosLo = sinf_poly_q1_4x(x2Lo);
    f64_2x cosHi = sinf_poly_q1_4x(x2Hi);

    f32_4x sins = F32_4x(sinLo, sinHi);
    f32_4x coss = F32_4x(cosLo, cosHi);

    f32_4x result = select4x(sins, sel0, coss);

    f32_4x signMask = S32_4x(MATS_F32_SIGN_MASK);
    f32_4x signMod = s32_4x_and(sel1, signMask);
    result = result ^ signMod;
    result = select4x(result, smallestMask, y);

    return result;
}

internal SinCos32_4x
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

    f64_2x xmLo = select2x(xModLo, smallMaskLo, xdLo);
    f64_2x xmHi = select2x(xModHi, smallMaskHi, xdHi);
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
    cosResult = select4x(coss, sel0, sins);
    sinResult = select4x(sins, sel0, coss);

    f32_4x signMask = S32_4x(MATS_F32_SIGN_MASK);
    f32_4x selCosSign = s32_4x_xor(sel0, sel1);
    f32_4x selSinSign = sel1;
    f32_4x minCos = cosResult ^ signMask;
    f32_4x minSin = sinResult ^ signMask;
    cosResult = select4x(cosResult, selCosSign, minCos);
    cosResult = select4x(cosResult, smallestMask, ones);
    sinResult = select4x(sinResult, selSinSign, minSin);
    sinResult = select4x(sinResult, smallestMask, y);

    SinCos32_4x result;
    result.cos = cosResult;
    result.sin = sinResult;
    return result;
}

internal f32_4x
tan32_kernel_4x(f32_4x x, f32_4x mod)
{
    f32_4x xOrig = x;
    f32_4x hx = x & S32_4x(MATS_F32_ABS_MASK);

    f32_4x specialMask = s32_4x_less(hx, S32_4x(0x31800000));
    specialMask = s32_4x_and(s32_4x_equal(s32_4x_from_f32_trunc(x), zero_f32_4x()), specialMask);
    f32_4x special = zero_f32_4x();
    {
        f32_4x hxModP1Mask = s32_4x_equal(s32_4x_or(hx, s32_4x_add(mod, S32_4x(1))), zero_f32_4x());
        hxModP1Mask = s32_4x_and(hxModP1Mask, specialMask);
        special = (F32_4x(1.0f) / absolute(x)) & hxModP1Mask;
        f32_4x workMask = and_not(specialMask, hxModP1Mask);
        f32_4x modIsOneMask = s32_4x_equal(mod, S32_4x(1));
        f32_4x modIsOne = select4x(F32_4x(-1.0f) / x, modIsOneMask, x);
        special = special | (modIsOne & workMask);
    }

    f32_4x hxGreat = s32_4x_less(S32_4x(0x3F2CA140), hx);
    f32_4x xLessZero = s32_4x_less(x, zero_f32_4x());
    xLessZero = s32_4x_and(hxGreat, s32_4x_and(xLessZero, S32_4x(MATS_F32_SIGN_MASK)));
    x = x ^ xLessZero;
    f32_4x zPre = F32_4x(gPiOver4F32_hi) - x;
    f32_4x wPre = F32_4x(gPiOver4F32_lo) - (F32_4x(gPiOver4F32_hi) - zPre - x);
    f32_4x xPre = zPre + wPre;

    x = select4x(x, hxGreat, xPre);
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

    f32_4x oneOrMinus = F32_4x(1.0f) ^ (xOrig & F32_4x(MATS_F32_SIGN_MASK));
    v = f32_4x_from_s32(mod);
    f32_4x greatRes = oneOrMinus * (v - F32_4x(2.0f) * (x - (result * result / (result + v) - r)));

    f32_4x nonModdedMask = s32_4x_equal(mod, S32_4x(1));
    f32_4x modded = result & S32_4x(0xFFFFF000);
    v = r - (modded - x);
    f32_4x a = F32_4x(-1.0f) / result;
    f32_4x t = a & S32_4x(0xFFFFF000);
    f32_4x s = F32_4x(1.0f) + t * modded;
    modded = t + a * (s + t * v);

    result = select4x(modded, nonModdedMask, result);
    result = select4x(result, hxGreat, greatRes);
    result = select4x(result, specialMask, special);

    return result;
}

internal f32_4x
tan32_4x(f32_4x x)
{
    // NOTE(michiel): absolute value of y should be smaller than 100.0f
    f32_4x absX = x & S32_4x(MATS_F32_ABS_MASK);
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

    f32_4x xIn = select4x(xMod, inputSmall, x);
    f32_4x nIn = select4x(nMod, inputSmall, S32_4x(1));

    f32_4x result = tan32_kernel_4x(xIn, nIn);

    return result;
}

// NOTE(michiel): 4 calls to acos32/asin32 is faster at the moment

#if 0
// NOTE(michiel): Templates used to convert to 4x
internal f32
acos32_temp(f32 x)
{
#define pi        3.1415925026e+00f
#define pio2_hi   1.5707962513e+00f
#define pio2_lo   7.5497894159e-08f

	s32 ix = (s32)MATS_F32U(x).u;
    u32 hx = ix & MATS_F32_ABS_MASK;

    b32 hxEqOne = hx == 0x3F800000;
    b32 hxLtHaf = hx <  0x3F000000;
    b32 xLtZero = ix <  0;

    f32 polyIn = (1.0f - absolute32(x)) * 0.5f;
    if (hxLtHaf) {
        polyIn = x * x;
    }
    poly(polyIn);

    f32 result;
    if (hxLtHaf)
    {
        result = pio2_hi - (x - (pio2_lo - x * r));
    }
    else
    {
        f32 s = sqrt32(polyIn);
        f32 w = s * r;
        f32 c = -pio2_lo;
        if (hxEqOne)
        {
            w = 0.0f;
            s = 0.0f;
        }
        else if (!xLtZero)
        {
            f32 df = MATS_F32U(MATS_F32U(s).u & 0xFFFFF000).f;
            c = (polyIn - df * df) / (s + df);
            s = df;
        }

        w = w + c;
        result = 2.0f * (s + w);
        f32 lessTZero = pi - result;
        result = hxEqOne ? 0.0f : result;
        result = xLtZero ? lessTZero : result;
    }

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

    s32 ix = (s32)MATS_F32U(x).u;
    u32 hx = ix & MATS_F32_ABS_MASK;

    b32 hxEqOne = hx == 0x3F800000;
    b32 hxGtOne = hx >  0x3F800000;
    b32 hxLtMax = hx <  0x3F79999A;
    b32 hxLtHaf = hx <  0x3F000000;

    f32 preW = 1.0f - MATS_F32U(hx).f;
    f32 polyIn = preW * 0.5f;
    if (hxLtHaf) {
        polyIn = x * x;
    }
    poly(polyIn);

    f32 result;
    if (hxGtOne)
    {
        result = (x - x) / (x - x);
    }
    else if (hxLtHaf)
    {
        result = x + x * r;
    }
    else
    {
        f32 s = sqrt32(polyIn);
        f32 sr = s * r;
        f32 a = pio2_hi;
        f32 b = 2.0f * (s + sr);
        f32 c = pio2_lo;
        if (hxEqOne)
        {
            b = 0;
        }
        else if (hxLtMax)
        {
            f32 w = s;
            u32 iw = MATS_F32U(w).u;
            w = MATS_F32U(iw & 0xFFFFF000).f;
            f32 correct = (polyIn - w * w) / (s + w);
            b = 2.0f * sr - (pio2_lo - 2.0f * correct);
            c = pio4_hi - 2.0f * w;
            a = pio4_hi;
        }
        result = a - (b - c);

        result = (ix < 0) ? -result : result;
    }

    return result;
#undef pio4_hi
#undef pio2_hi
#undef pio2_lo
}
#endif

internal f32_4x
acos32_4x(f32_4x x)
{
    // NOTE(michiel): Function is unspecified out of the [-1, 1] domain
    f32_4x pi4x     = F32_4x(3.1415925026e+00f);
    f32_4x pio2Hi4x = F32_4x(1.5707962513e+00f);
    f32_4x pio2Lo4x = F32_4x(7.5497894159e-08f);
    f32_4x zero4x   = F32_4x(0.0f);
    f32_4x one4x    = F32_4x(1.0f);
    f32_4x two4x    = F32_4x(2.0f);
    f32_4x half4x   = F32_4x(0.5f);

    f32_4x preMask = S32_4x(0x7FFFFFFF);
    f32_4x hx = x & preMask;

    preMask.mi = _mm_slli_epi32(preMask.mi, 12);

    f32_4x hxEqOne  = s32_4x_equal(hx, one4x);
    f32_4x hxLtHalf = s32_4x_less(hx, half4x);
    f32_4x xLtZero  = s32_4x_less(x, zero4x);

    f32_4x polyIn   = (one4x - hx) * half4x;
    polyIn = select4x(polyIn, hxLtHalf, x * x);

    f32_4x s = sqrt32_4x(polyIn);
    //f32_4x df = s & S32_4x(0xFFFFF000);
    f32_4x df = s & preMask;
    f32_4x cMod = (polyIn - df * df) / (s + df);

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

    f32_4x r = p / q;

    f32_4x smallResult = pio2Hi4x - (x - (pio2Lo4x - x * r));

    f32_4x w = s * r;
    f32_4x c = -pio2Lo4x;

    f32_4x xLtZeroOrOne = s32_4x_or(xLtZero, hxEqOne);
    s = select4x(df, xLtZeroOrOne, s);
    c = select4x(cMod, xLtZeroOrOne, c);

    // NOTE(michiel): Zero out w and s, if |x| is 1.0
    w = and_not(w, hxEqOne);
    s = and_not(s, hxEqOne);

    w = w + c;
    f32_4x result = two4x * (s + w);
    f32_4x negResult = pi4x - result;
    result = and_not(result, hxEqOne);
    result = select4x(result, xLtZero, negResult);
    result = select4x(result, hxLtHalf, smallResult);

    return result;
}

internal f32_4x
asin32_4x(f32_4x x)
{
    // NOTE(michiel): Function is unspecified out of the [-1, 1] domain
    f32_4x one4x  = F32_4x(1.0f);
    f32_4x two4x  = F32_4x(2.0f);
    f32_4x half4x = F32_4x(0.5f);
    f32_4x piOver2Hi = F32_4x(1.57079637050628662109375f);
    f32_4x piOver2Lo = F32_4x(-4.37113900018624283e-8f);
    f32_4x piOver4Hi = F32_4x(0.785398185253143310546875f);

    f32_4x preMask = S32_4x(0x7FFFFFFF);

    f32_4x hx = x & preMask;
    f32_4x xSign = and_not(x, preMask);

    preMask.mi = _mm_slli_epi32(preMask.mi, 12);

    f32_4x hxEqOne  = s32_4x_equal(hx, one4x);
    f32_4x hxLtMax  = s32_4x_less(hx, S32_4x(0x3F79999A));
    f32_4x hxLtHalf = s32_4x_less(hx, half4x);

    f32_4x polyIn = select4x(half4x * (one4x - hx), hxLtHalf, x * x);

#if MATS_ASINCOS_USE_SMALL_POLY
    f32_4x p = F32_4x(-8.6563630030e-03f) * polyIn;
    f32_4x q = F32_4x(-7.0662963390e-01f) * polyIn;
    p = (p - F32_4x(4.2743422091e-02f)) * polyIn;
    q = (q + one4x);
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
    q = (q + one4x);
    p = (p + F32_4x(1.6666667163e-01f)) * polyIn;
#endif

    f32_4x r = p / q;

    f32_4x ltHalfRes = x + x * r;

    f32_4x extraRes = ltHalfRes;
    f32_4x extraMask = hxLtHalf;

    f32_4x s = sqrt32_4x(polyIn);
    f32_4x sr = s * r;

    //f32_4x w = s & S32_4x(0xFFFFF000);
    f32_4x w = s & preMask;
    f32_4x correct = (polyIn - w * w) / (s + w);

    //preMask.mi = _mm_srli_epi32(_mm_slli_epi32(preMask.mi, 19), 8);

    f32_4x bLarge = and_not(two4x * (s + sr), hxEqOne);
    f32_4x bSmall = two4x * sr - (piOver2Lo - two4x * correct);
    f32_4x cSmall = piOver4Hi - two4x * w;

    //f32_4x a = s32_4x_or(piOver4Hi, s32_4x_and_not(S32_4x(0x00800000), hxLtMax)); // select4x(piOver2Hi, hxLtMax, piOver4Hi);
    //f32_4x a = piOver4Hi | and_not(preMask, hxLtMax);
    f32_4x a = select4x(piOver2Hi, hxLtMax, piOver4Hi);
    f32_4x b = select4x(bLarge, hxLtMax, bSmall);
    f32_4x c = select4x(piOver2Lo, hxLtMax, cSmall);

    f32_4x result = a - (b - c);
    result = result ^ xSign;
    result = select4x(result, extraMask, extraRes);
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

    f32_4x xSign          = x & F32_4x(MATS_F32_SIGN_MASK);
    f32_4x lessThanMax    = s32_4x_less(hx, S32_4x(0x50800000));

    // NOTE(michiel): x >= 0x50800000
    f32_4x aBig = F32_4x(1.5707962513e+00f);
    f32_4x bBig = F32_4x(7.5497894159e-08f);
    aBig = aBig | xSign; // select(aBig, lessThanZero, -aBig);
    bBig = bBig | xSign; // select(bBig, lessThanZero, -bBig);

#if 0
    f32_4x isNan          = s32_4x_greater(hx, S32_4x(MATS_F32_EXP_MASK));
    aBig = select(aBig, isNan, x);
    bBig = select(bBig, isNan, x);
#endif

    f32_4x bigResult = aBig + bBig;

    // NOTE(michiel): Small result
    f32_4x smallestMask = s32_4x_less(hx, S32_4x(0x31000000));
    x  = and_not(x, smallestMask);
    hx = and_not(hx, smallestMask);

    f32_4x smallMask    = s32_4x_less(hx, S32_4x(0x3EE00000));
    f32_4x middlestMask = s32_4x_less(hx, S32_4x(0x3F300000));
    f32_4x middleMask   = s32_4x_less(hx, S32_4x(0x3F980000));
    f32_4x notBigMask   = s32_4x_less(hx, S32_4x(0x401C0000));

#if 0
    s32 id = 0;
    f32 a = 0.0f;
    f32 b = 1.0f;
    if (middleMask) {
        a += 1.0f;
    }
    if (middlestMask) {
        a += 1.0f;
    }
    if (notBigMask && !middleMask) {
        a += 1.0f;
        b += 0.5f;
    }

    f32 nom = a * x - b;
    f32 den = a + b * x;
    x = nom / den;
    id = (middleMask ? 0x00 : 0x08) | ((middlestMask | notBigMask) ? 0x00 : 0x04);
#endif

    f32_4x aSmall = zero_f32_4x();
    f32_4x bSmall = F32_4x(1.0f);

    f32_4x half4x = F32_4x(0.5f);

    aSmall = aSmall + (bSmall & notBigMask);

    notBigMask = s32_4x_and_not(notBigMask, middleMask);
    aSmall = aSmall + (bSmall & middlestMask);

    bSmall = bSmall + (half4x & notBigMask);

    f32_4x mod = S32_4x(0x04040404);
    f32_4x ids = S32_4x(0x03020100);
    ids = s32_4x_or(ids, s32_4x_and_not(mod, s32_4x_or(middlestMask, notBigMask)));
    mod.mi = _mm_slli_epi32(mod.mi, 1);
    ids = s32_4x_or(ids, s32_4x_and_not(mod, middleMask));

    f32_4x nom = aSmall * hx - bSmall;
    f32_4x den = aSmall + bSmall * hx;
    f32_4x xMod = nom / den;

    x = select4x(xMod, smallMask, x);

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
    result = result ^ xSign;

    result = select4x(result, smallMask, x - smallResult);
    result = select4x(bigResult, lessThanMax, result);

    return result;
}

#if 0
// NOTE(michiel): Atan2 template used for simd
internal f32
atan2_32_temp(f32 y, f32 x)
{
    s32 ix = MATS_S32_FROM_F32(x);
    u32 hx = ix & MATS_F32_ABS_MASK;
    s32 iy = MATS_S32_FROM_F32(y);
    u32 hy = iy & MATS_F32_ABS_MASK;

    b32 xIsOne  = ix == 0x3F800000;
    b32 xLtZero = ix < 0;
    b32 yLtZero = iy < 0;
    b32 xIsZero = hx == 0;
    b32 yIsZero = hy == 0;
    b32 xIsInf  = hx == MATS_F32_EXP_MASK;
    b32 yIsInf  = hy == MATS_F32_EXP_MASK;

    f32 atanInput  = xIsOne ? y : absolute32(y / x);
    f32 atanOutput = atan32(atanInput);

    f32 result;
    if (xIsOne)
    {
        result = atanOutput;
    }
    else
    {
        if (yIsZero || xIsZero || xIsInf || yIsInf)
        {
            f32 mult = 2.0f;
            if (xIsInf || yIsZero)
            {
                mult = xLtZero ? 4.0f : 0.0f;
            }
            if (xIsInf && yIsInf)
            {
                mult = xLtZero ? 3.0f : 1.0f;
            }
            result = mult * gPiOver4F32;
        }
        else
        {
            result = atanOutput;

            if (xLtZero) {
                result = gPiF32 - (result - gPiF32_lo);
            }
        }
        result = yLtZero ? -result : result;
    }

    return result;
}
#endif

internal f32_4x
atan2_32_4x(f32_4x y, f32_4x x)
{
    f32_4x hx = x & F32_4x(0x7FFFFFFFU);
    f32_4x hy = y & F32_4x(0x7FFFFFFFU);
    f32_4x ySign = and_not(y, F32_4x(0x7FFFFFFFU));

    f32_4x xIsOne  = s32_4x_equal(x, S32_4x(0x3F800000));
    f32_4x xLtZero = s32_4x_less(x, zero_s32_4x());
    f32_4x xIsZero = s32_4x_equal(hx, zero_f32_4x());
    f32_4x yIsZero = s32_4x_equal(hy, zero_f32_4x());
    f32_4x xIsInf  = s32_4x_equal(hx, S32_4x(MATS_F32_EXP_MASK));
    f32_4x yIsInf  = s32_4x_equal(hy, S32_4x(MATS_F32_EXP_MASK));

    f32_4x atanIn  = select4x(absolute(y / x), xIsOne, y);
    f32_4x atanOut = atan32_4x(atanIn);

    f32_4x doSpecial = yIsZero | xIsZero | xIsInf | yIsInf;
    f32_4x do40      = xIsInf | yIsZero;
    f32_4x do31      = xIsInf & yIsInf;
    f32_4x piMult    = F32_4x(2.0f);
    f32_4x do40mul   = s32_4x_and(F32_4x(4.0f), xLtZero);
    f32_4x do31mul   = s32_4x_and(F32_4x(2.0f), xLtZero) + F32_4x(1.0f);
    piMult = select4x(piMult, do40, do40mul);
    piMult = select4x(piMult, do31, do31mul);
    f32_4x specialRes = piMult * F32_4x(gPiOver4F32);

    f32_4x result = atanOut;
    f32_4x negResult = F32_4x(gPiF32) - (result - F32_4x(gPiF32_lo));
    result = select4x(result, xLtZero, negResult);
    result = select4x(result, doSpecial, specialRes);
    result = result ^ ySign;
    result = select4x(result, xIsOne, atanOut);

    return result;
}
