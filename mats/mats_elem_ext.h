//
// NOTE(michiel): Pow
//

internal f64
pow32_log2(u32 xu)
{
    u32 temp = xu - 0x3F330000;
    u32 index = (temp >> (MATS_F32_EXP_SHIFT - POW32_LOG2_TABLE_BITS)) % (1 << POW32_LOG2_TABLE_BITS);
    u32 top = temp & 0xFF800000;
    u32 iz = xu - top;
    s32 k = (s32)temp >> MATS_F32_EXP_SHIFT;
    f64 invC = gPowF32_Log2Table[index].invC;
    f64 logC = gPowF32_Log2Table[index].logC;
    f64 z = (f64)MATS_F32U(iz).f;
    f64 r = z * invC - 1.0;
    f64 y0 = logC + (f64)k;

    f64 r2 = r * r;
    f64 result = gPowF32_Log2Poly[0] * r + gPowF32_Log2Poly[1];
    f64 p = gPowF32_Log2Poly[2] * r + gPowF32_Log2Poly[3];
    f64 r4 = r2 * r2;
    f64 q = gPowF32_Log2Poly[4] * r + y0;
    q = p * r2 + q;
    result = result * r4 + q;
    return result;
}

internal f64
pow32_exp2(f64 xd, u32 signBias)
{
    f64 kd = xd + gExp2F32_ShiftScaled;
    u64 ki = MATS_F64U(kd).u;
    kd -= gExp2F32_ShiftScaled;

    f64 r = xd - kd;
    u64 t = gExp2F32_Table[ki % (1 << EXP2_32_TABLE_BITS)];
    u64 ski = ki + signBias;
    t += ski << (52 - EXP2_32_TABLE_BITS);
    f64 s = MATS_F64U(t).f;
    f64 z = gExp2F32_Poly[0] * r + gExp2F32_Poly[1];
    f64 r2 = r * r;
    f64 result = gExp2F32_Poly[2] * r + 1.0;
    result = z * r2 + result;
    result = result * s;
    return result;
}

internal s32
pow32_checkint(u32 xu)
{
    s32 e = (xu & MATS_F32_EXP_MASK) >> MATS_F32_EXP_SHIFT;
    if (e < MATS_F32_EXP_BIAS) {
        return 0;
    }
    if (e > MATS_F32_EXP_BIAS + MATS_F32_EXP_SHIFT) {
        return 2;
    }
    if (xu & ((1 << (MATS_F32_EXP_BIAS + MATS_F32_EXP_SHIFT - e)) - 1)) {
        return 0;
    }
    if (xu & (1 << (MATS_F32_EXP_BIAS + MATS_F32_EXP_SHIFT - e))) {
        return 1;
    }
    return 2;
}

internal s32
pow32_zeroinfnan(u32 ix)
{
    return 2 * ix - 1 >= 2u * MATS_F32_EXP_MASK - 1;
}

internal f32
pow32(f32 x, f32 y)
{
    u32 signBias = 0;

    u32 xu = MATS_F32U(x).u;
    u32 yu = MATS_F32U(y).u;

    if (((xu - 0x00800000) >= (MATS_F32_EXP_MASK - 0x00800000)) || pow32_zeroinfnan(yu))
    {
        /* Either (x < 0x1p-126 or inf or nan) or (y is 0 or inf or nan).  */
        if (pow32_zeroinfnan(yu))
        {
            if (2 * yu == 0) {
                return issignaling32(x) ? x + y : 1.0f;
            }
            if (xu == 0x3F800000) {
                return issignaling32(y) ? x + y : 1.0f;
            }
            if (((2 * xu) > (2u * MATS_F32_EXP_MASK)) ||
                ((2 * yu) > (2u * MATS_F32_EXP_MASK))) {
                return x + y;
            }
            if ((2 * xu) == 0x3F800000) {
                return 1.0f;
            }
            if (((2 * xu) < (2 * 0x3F800000)) == !(yu & MATS_F32_SIGN_MASK)) {
                return 0.0f; // NOTE(michiel): |x| < 1 && y == inf or |x| > 1 && y == -inf
            }
            return y * y;
        }
        if (pow32_zeroinfnan(xu))
        {
            f32 x2 = x * x;
            if ((xu & MATS_F32_SIGN_MASK) &&
                (pow32_checkint(yu) == 1))
            {
                x2 = -x2;
                signBias = 1;
            }
            return yu & MATS_F32_SIGN_MASK? 1.0f / x2 : x2;
        }
        /* x and y are non-zero finite.  */
        if (xu & MATS_F32_SIGN_MASK)
        {
            s32 yInt = pow32_checkint(yu);
            if (yInt == 0) {
                return mats_invalid32(x);
            }
            if (yInt == 1) {
                signBias = (1 << (EXP2_32_TABLE_BITS + 11));
            }
            xu &= MATS_F32_ABS_MASK;
        }
        if (xu < 0x00800000)
        {
            xu = MATS_F32U(x * 0x1p23f).u;
            xu &= MATS_F32_ABS_MASK;
            xu -= MATS_F32_EXP_SHIFT << MATS_F32_EXP_SHIFT;
        }
    }

    f64 logX = pow32_log2(xu);
    f64 yLogX = (f64)y * logX; /* Note: cannot overflow, y is single prec.  */
    if ((MATS_F64U(yLogX).u >> 47 & 0xFFFF) >= (MATS_F64U(126.0 * POW32_SCALE).u >> 47))
    {
        if (yLogX > 0x1.fffffffd1d571p+6 * POW32_SCALE) {
            return mats_overflow32(signBias);
        }
        if (yLogX <= -150.0 * POW32_SCALE) {
            return mats_underflow32(signBias);
        }
    }
    return (f32)pow32_exp2(yLogX, signBias);
}

internal f32
pow10_32(f32 x)
{
    return pow32(10.0f, x);
}

internal f32
exp10_32(f32 x)
{
    return pow32(10.0f, x);
}

/* Compute y+TAIL = log(x) where the rounded result is y and TAIL has about
   additional 15 bits precision.  IX is the bit representation of x, but
   normalized in the subnormal range using the sign bit for the exponent.  */
internal f64
pow64_log(u64 ix, f64 *tail)
{
#define POW64_LOG_N (1 << POW64_LOG_TABLE_BITS)
    /* x = 2^k z; where z is in range [OFF,2*OFF) and exact. (OFF = 0x3fe6955500000000)
       The range is split into N subintervals.
       The ith subinterval contains z and c is near its center.  */
    u64 tmp = ix - 0x3FE6955500000000ULL;
    s32 i = (tmp >> (MATS_F64_EXP_SHIFT - POW64_LOG_TABLE_BITS)) % POW64_LOG_N;
    s32 k = (s64)tmp >> MATS_F64_EXP_SHIFT; /* arithmetic shift */
    u64 iz = ix - (tmp & (0xFFFULL << 52));
    f64 z = MATS_F64_FROM_U64(iz);
    f64 kd = (f64)k;

    /* log(x) = k*Ln2 + log(c) + log1p(z/c-1).  */
    f64 invc = gPowLogData64.tab[i].invc;
    f64 logc = gPowLogData64.tab[i].logc;
    f64 logctail = gPowLogData64.tab[i].logctail;

    /* Note: 1/c is j/N or j/N/2 where j is an integer in [N,2N) and
       |z/c - 1| < 1/N, so r = z/c - 1 is exactly representible.  */
#if MATS_HAVE_FAST_FMA
    f64 r = mats_fma(z, invc, -1.0);
#else
    /* Split z such that rhi, rlo and rhi*rhi are exact and |rlo| <= |r|.  */
    f64 zhi = MATS_F64_FROM_U64((iz + (1ULL << 31)) & (-1ULL << 32));
    f64 zlo = z - zhi;
    f64 rhi = zhi * invc - 1.0;
    f64 rlo = zlo * invc;
    f64 r = rhi + rlo;
#endif

    /* k*Ln2 + log(c) + r.  */
    f64 t1 = kd * gPowLogData64.ln2hi + logc;
    f64 t2 = t1 + r;
    f64 lo1 = kd * gPowLogData64.ln2lo + logctail;
    f64 lo2 = t1 - t2 + r;

    /* Evaluation is optimized assuming superscalar pipelined execution.  */
    f64 ar = gPowLogData64.poly[0] * r; /* A[0] = -0.5.  */
    f64 ar2 = r * ar;
    f64 ar3 = r * ar2;
    /* k*Ln2 + log(c) + r + A[0]*r*r.  */
#if MATS_HAVE_FAST_FMA
    f64 hi = t2 + ar2;
    f64 lo3 = mats_fma(ar, r, -ar2);
    f64 lo4 = t2 - hi + ar2;
#else
    f64 arhi = gPowLogData64.poly[0] * rhi;
    f64 arhi2 = rhi * arhi;
    f64 hi = t2 + arhi2;
    f64 lo3 = rlo * (ar + arhi);
    f64 lo4 = t2 - hi + arhi2;
#endif
    /* p = log1p(r) - r - A[0]*r*r.  */

    f64 p = (ar3 * (gPowLogData64.poly[1] + r * gPowLogData64.poly[2] + ar2 * (gPowLogData64.poly[3] + r * gPowLogData64.poly[4] + ar2 * (gPowLogData64.poly[5] + r * gPowLogData64.poly[6]))));

    f64 lo = lo1 + lo2 + lo3 + lo4 + p;
    f64 y = hi + lo;
    *tail = hi - y + lo;
    return y;
#undef POW64_LOG_N
}

/* Computes sign*exp(x+xtail) where |xtail| < 2^-8/N and |xtail| <= |x|.
   The sign_bias argument is SIGN_BIAS or 0 and sets the sign to -1 or 1.  */
internal f64
pow64_exp(f64 x, f64 xtail, u32 signBias)
{
#define POW64_EXP_N   (1 << EXP64_TABLE_BITS)

    u32 abstop = mats_top12(x) & 0x7FF;
    if ((abstop - mats_top12(0x1p-54)) >= (mats_top12(512.0) - mats_top12(0x1p-54)))
    {
        if ((abstop - mats_top12(0x1p-54)) >= 0x80000000)
        {
            /* Avoid spurious underflow for tiny x.  */
            /* Note: 0 is common input.  */
            f64 one = 1.0 + x;
            return signBias ? -one : one;
        }
        if (abstop >= mats_top12(1024.0))
        {
            /* Note: inf and nan are already handled.  */
            if (MATS_U64_FROM_F64(x) >> 63)
                return mats_underflow64(signBias);
            else
                return mats_overflow64(signBias);
        }
        /* Large x is special cased below.  */
        abstop = 0;
    }

    /* exp(x) = 2^(k/N) * exp(r), with exp(r) in [2^(-1/2N),2^(1/2N)].  */
    /* x = ln2/N*k + r, with int k and r in [-ln2/2N, ln2/2N].  */
    f64 z = gExp64Data.invln2N * x;
#if 0 && TOINT_INTRINSICS
    f64 kd = roundtoint(z);
    u64 ki = converttoint(z);
#elif EXP64_USE_TOINT_NARROW
    /* z - kd is in [-0.5-2^-16, 0.5] in all rounding modes.  */
    f64 kd = z + gExp64Data.shift;
    u64 ki = MATS_U64_FROM_F64(kd) >> 16;
    kd = (f64)(s32)ki;
#else
    /* z - kd is in [-1, 1] in non-nearest rounding modes.  */
    f64 kd = z + gExp64Data.shift;
    u64 ki = MATS_U64_FROM_F64(kd);
    kd -= gExp64Data.shift;
#endif

    f64 r = x + kd * gExp64Data.negln2hiN+ kd * gExp64Data.negln2loN;
    /* The code assumes 2^-200 < |xtail| < 2^-8/N.  */
    r += xtail;
    /* 2^(k/N) ~= scale * (1 + tail).  */
    u64 idx = 2 * (ki % POW64_EXP_N);
    u64 top = (ki + signBias) << (MATS_F64_EXP_SHIFT - EXP64_TABLE_BITS);
    f64 tail = MATS_F64_FROM_U64(gExp64Data.tab[idx]);
    /* This is only a valid scale when -1023*N < k < 1024*N.  */
    u64 sbits = gExp64Data.tab[idx + 1] + top;
    /* exp(x) = 2^(k/N) * exp(r) ~= scale + scale * (tail + exp(r) - 1).  */
    /* Evaluation is optimized assuming superscalar pipelined execution.  */
    f64 r2 = r * r;
    /* Without fma the worst case error is 0.25/N ulp larger.  */
    /* Worst case error is less than 0.5+1.11/N+(abs poly error * 2^53) ulp.  */

#define C2 gExp64Data.poly[5 - EXP64_POLY_ORDER]
#define C3 gExp64Data.poly[6 - EXP64_POLY_ORDER]
#define C4 gExp64Data.poly[7 - EXP64_POLY_ORDER]
#define C5 gExp64Data.poly[8 - EXP64_POLY_ORDER]
#define C6 gExp64Data.poly[9 - EXP64_POLY_ORDER]

#if EXP64_POLY_ORDER == 4
    f64 tmp = tail + r + r2 * C2 + r * r2 * (C3 + r * C4);
#elif EXP64_POLY_ORDER == 5
    f64 tmp = tail + r + r2 * (C2 + r * C3) + r2 * r2 * (C4 + r * C5);
#elif EXP64_POLY_ORDER == 6
    f64 tmp = tail + r + r2 * (0.5 + r * C3) + r2 * r2 * (C4 + r * C5 + r2 * C6);
#endif

#undef C2
#undef C3
#undef C4
#undef C5
#undef C6

    if (abstop == 0) {
        /* Handle cases that may overflow or underflow when computing the result that
           is scale*(1+TMP) without intermediate rounding.  The bit representation of
           scale is in SBITS, however it has a computed exponent that may have
           overflown into the sign bit so that needs to be adjusted before using it as
           a double.  (int32_t)KI is the k used in the argument reduction and exponent
           adjustment of scale, positive k here means the result may overflow and
           negative k means the result may underflow.  */
        if ((ki & 0x80000000) == 0)
        {
            /* k > 0, the exponent of scale might have overflowed by <= 460.  */
            sbits -= 1009ULL << 52;
            f64 scale = MATS_F64_FROM_U64(sbits);
            f64 y = 0x1p1009 * (scale + scale * tmp);
            return y;
        }
        /* k < 0, need special care in the subnormal range.  */
        sbits += 1022ULL << 52;
        /* Note: sbits is signed scale.  */
        f64 scale = MATS_F64_FROM_U64(sbits);
        f64 y = scale + scale * tmp;
        if (absolute64(y) < 1.0)
        {
            /* Round y to the right precision before scaling it into the subnormal
           range to avoid double rounding that can cause 0.5+E/2 ulp error where
           E is the worst-case ulp error outside the subnormal range.  So this
           is only useful if the goal is better than 1 ulp worst-case error.  */
            f64 one = (y < 0.0) ? -1.0 : 1.0;
            f64 lo = scale - y + scale * tmp;
            f64 hi = one + y;
            lo = one - hi + y + lo;
            y = (hi + lo) - one;
            /* Fix the sign of 0.  */
            if (y == 0.0) {
                y = MATS_F64_FROM_U64(sbits & 0x8000000000000000);
            }
        }
        y = 0x1p-1022 * y;
        return y;
    }
    f64 scale = MATS_F64_FROM_U64(sbits);
    /* Note: tmp == 0 or |tmp| > 2^-200 and scale > 2^-739, so there
       is no spurious underflow here even without fma.  */
    return scale + scale * tmp;
#undef POW64_EXP_N
}

/* Returns 0 if not int, 1 if odd int, 2 if even int.  The argument is
   the bit representation of a non-zero finite floating-point value.  */
internal s32
pow64_checkint(u64 iy)
{
    s32 e = (iy >> 52) & 0x7ff;
    if (e < 0x3ff) {
        return 0;
    } else if (e > 0x3ff + 52) {
        return 2;
    } else if (iy & ((1ULL << (0x3ff + 52 - e)) - 1)) {
        return 0;
    } else if (iy & (1ULL << (0x3ff + 52 - e))) {
        return 1;
    }
    return 2;
}

/* Returns 1 if input is the bit representation of 0, infinity or nan.  */
internal s32
pow64_zeroinfnan(u64 i)
{
    return (2 * i - 1) >= (2 * MATS_U64_FROM_F64(F64_INF) - 1);
}

internal f64
pow64(f64 x, f64 y)
{
    u32 signBias = 0;

    u64 ix = MATS_U64_FROM_F64(x);
    u64 iy = MATS_U64_FROM_F64(y);
    u32 topx = mats_top12(x);
    u32 topy = mats_top12(y);

    if (((topx - 0x001) >= (0x7ff - 0x001) ||
         ((topy & 0x7ff) - 0x3be) >= (0x43e - 0x3be)))
    {
        /* Note: if |y| > 1075 * ln2 * 2^53 ~= 0x1.749p62 then pow(x,y) = inf/0
       and if |y| < 2^-54 / 1075 ~= 0x1.e7b6p-65 then pow(x,y) = +-1.  */
        /* Special cases: (x < 0x1p-126 or inf or nan) or
       (|y| < 0x1p-65 or |y| >= 0x1p63 or nan).  */
        if (pow64_zeroinfnan(iy))
        {
            if (2 * iy == 0) {
                return issignaling64(x) ? x + y : 1.0;
            }
            if (ix == MATS_U64_FROM_F64(1.0)) {
                return issignaling64(y) ? x + y : 1.0;
            }
            if (((2 * ix) > (2 * MATS_U64_FROM_F64(F64_INF))) ||
                ((2 * iy) > (2 * MATS_U64_FROM_F64(F64_INF)))) {
                return x + y;
            }
            if ((2 * ix) == (2 * MATS_U64_FROM_F64(1.0))) {
                return 1.0;
            }
            if (((2 * ix) < (2 * MATS_U64_FROM_F64(1.0))) == !(iy >> 63)) {
                return 0.0; /* |x|<1 && y==inf or |x|>1 && y==-inf.  */
            }
            return y * y;
        }
        if (pow64_zeroinfnan(ix))
        {
            f64 x2 = x * x;
            if ((ix >> 63) && (pow64_checkint(iy) == 1))
            {
                x2 = -x2;
                signBias = 1;
            }
            // TODO(michiel): Apparently:
            // Without the barrier some versions of clang hoist the 1/x2 and
            // thus division by zero exception can be signaled spuriously.
            // opt_barrier_double(1.0 / x2)
            return (iy >> 63) ? 1.0 / x2 : x2;
        }
        /* Here x and y are non-zero finite.  */
        if (ix >> 63)
        {
            /* Finite x < 0.  */
            s32 yint = pow64_checkint(iy);
            if (yint == 0) {
                return mats_invalid64(x);
            }
            if (yint == 1) {
                signBias = (0x800 << EXP64_TABLE_BITS);
            }
            ix &= 0x7FFFFFFFFFFFFFFFULL;
            topx &= 0x7FF;
        }

        if (((topy & 0x7FF) - 0x3BE) >= (0x43E - 0x3BE))
        {
            /* Note: signBias == 0 here because y is not odd.  */
            if (ix == MATS_U64_FROM_F64(1.0)) {
                return 1.0;
            }
            if ((topy & 0x7FF) < 0x3BE)
            {
                /* |y| < 2^-65, x^y ~= 1 + y*log(x).  */
                return (ix > MATS_U64_FROM_F64(1.0)) ? 1.0 + y : 1.0 - y;
            }
            return (ix > MATS_U64_FROM_F64(1.0)) == (topy < 0x800) ? mats_overflow64(0) : mats_underflow64(0);
        }
        if (topx == 0)
        {
            /* Normalize subnormal x so exponent becomes negative.  */
            ix = MATS_U64_FROM_F64(x * 0x1p52);
            ix &= 0x7FFFFFFFFFFFFFFFULL;
            ix -= 52ULL << 52;
        }
    }

    f64 lo;
    f64 hi = pow64_log(ix, &lo);
#if MATS_HAVE_FAST_FMA
    f64 ehi = y * hi;
    f64 elo = y * lo + mats_fma(y, hi, -ehi);
#else
    f64 yhi = MATS_F64_FROM_U64(iy & (-1ULL << 27));
    f64 ylo = y - yhi;
    f64 lhi = MATS_F64_FROM_U64(MATS_U64_FROM_F64(hi) & (-1ULL << 27));
    f64 llo = hi - lhi + lo;
    f64 ehi = yhi * lhi;
    f64 elo = ylo * lhi + y * llo; /* |elo| < |ehi| * 2^-25.  */
#endif
    return pow64_exp(ehi, elo, signBias);
}
