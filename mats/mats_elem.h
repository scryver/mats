//
// NOTE(michiel): Sqrt/Hypot/Exp/Exp2/Log/Log2/Log10
//

//
// NOTE(michiel): 32-bit
//

internal f32
sqrt32(f32 x)
{
    f32 result;
#if MATS_USE_SSE2
    WideMath m;
    m.m = _mm_sqrt_ss(_mm_set_ss(x));
    result = m.e[0];
#else
	s32 ix = (s32)MATS_F32U(x).u;
    u32 hx = ix & MATS_F32_ABS_MASK;

    /* take care of Inf and NaN */
	if(!MATS_F32_UWORD_IS_FINITE(hx)) {
	    return x * x + x;		// sqrt(NaN)=NaN, sqrt(+inf)=+inf, sqrt(-inf)=sNaN
    }
    /* take care of zero and -ves */
	if (MATS_F32_UWORD_IS_ZERO(hx)) {
        return x; // sqrt(+-0) = +-0
    }
	if (ix < 0) {
        return (x - x) / (x - x);		/* sqrt(-ve) = sNaN */
    }

    /* normalize x */
	s32 m = (ix >> MATS_F32_EXP_SHIFT);
	if (MATS_F32_UWORD_IS_SUBNORMAL(hx))
    {   /* subnormal x */
        s32 i;
	    for (i = 0; (ix & 0x00800000) == 0; ++i) {
            ix <<= 1;
        }
	    m -= i - 1;
	}
	m -= MATS_F32_EXP_BIAS;	/* unbias exponent */
	ix = (ix & MATS_F32_MANT_MASK) | 0x00800000;
	if (m & 1) {	/* odd m, double x to make it even */
	    ix += ix;
    }
	m >>= 1;	/* m = [m/2] */

    /* generate sqrt(x) bit by bit */
	ix += ix;
	s32 q = 0;        /* q = sqrt(x) */
    s32 s = 0;
	u32 r = 0x01000000;  /* r = moving bit from right to left */

	while (r != 0)
    {
	    s32 t = s + r;
	    if (t <= ix)
        {
            s    = t + r;
            ix  -= t;
            q   += r;
	    }
	    ix += ix;
	    r >>= 1;
	}

    /* use floating add to find out rounding direction */
	if (ix != 0)
    {
	    f32 z = 1.0f - gTinyF32; /* trigger inexact flag */
	    if (z >= 1.0f)
        {
	        z = 1.0f + gTinyF32;
            if (z > 1.0f) {
                q += 2;
            } else {
                q += (q & 1);
            }
	    }
	}
	ix = (q >> 1) + 0x3f000000L;
	ix += (m << MATS_F32_EXP_SHIFT);

    result = MATS_F32U((u32)ix).f;
#endif
    return result;
}

internal f32
hypot32(f32 x, f32 y)
{
    u32 ux = MATS_F32U(x).u & MATS_F32_ABS_MASK;
    u32 uy = MATS_F32U(y).u & MATS_F32_ABS_MASK;

	if (ux < uy) {
        u32 temp = ux;
        ux = uy;
        uy = temp;
	}

    x = MATS_F32U(ux).f;
    y = MATS_F32U(uy).f;
    if (uy == MATS_F32_EXP_MASK) {
        return y;
    }
    if ((ux >= MATS_F32_EXP_MASK) ||
        (uy == 0) ||
        ((ux - uy) >= (25 << MATS_F32_EXP_SHIFT))) {
        return x + y;
    }

    f32 z = 1.0f;
    if (ux >= ((MATS_F32_EXP_BIAS + 60) << MATS_F32_EXP_SHIFT))
    {
        z = 0x1p90f;
        x *= 0x1p-90f;
        y *= 0x1p-90f;
    }
    else if (uy < ((MATS_F32_EXP_BIAS - 60) << MATS_F32_EXP_SHIFT))
    {
        z = 0x1p-90f;
        x *= 0x1p90f;
        y *= 0x1p90f;
    }
    f32 result = z * sqrt32((f32)((f64)x * x + (f64)y * y));
    //f32 result = z * sqrt32(x * x + y * y); // NOTE(michiel): Less accurate, but faster
    return result;
}

internal f32
exp32(f32 x)
{
    u32 xu = MATS_F32U(x).u;

    u32 abstop = abstop12_(x);
    if (abstop >= abstop12_(88.0f))
    {
        if (xu == MATS_F32U(-F32_INF).u) {
            return 0.0f;
        }
        if (abstop >= abstop12_(F32_INF)) {
            return x + x;
        }
        if (x > 0x1.62E42Ep6f) {
            return mats_overflow32(0);
        }
        if (x < -0x1.9FE368p6f) {
            return mats_underflow32(0);
        }
    }

    f64 xd = x;
    f64 z = gExp2F32_InvLn2Scaled * xd;

    f64 kd = z + gExp2F32_Shift; /* Rounding to double precision is required.  */
    u64 ki = MATS_F64U(kd).u;
    kd -= gExp2F32_Shift;

    f64 r = z - kd;

    /* exp(x) = 2^(k/N) * 2^(r/N) ~= s * (C0*r^3 + C1*r^2 + C2*r + 1) */
    u64 t = gExp2F32_Table[ki % (1 << EXP2_32_TABLE_BITS)];
    t += ki << (52 - EXP2_32_TABLE_BITS);

    f64 s = MATS_F64U(t).f;
    f64 res = gExp2F32_PolyScaled[0] * r + gExp2F32_PolyScaled[1];
    f64 r2 = r * r;
    f64 y = gExp2F32_PolyScaled[2] * r + 1.0;
    res = res * r2 + y;
    res = res * s;

    return (f32)res;
}

internal f32
exp2_32(f32 x)
{
    // NOTE(michiel): result = 2^x
    u32 xu = MATS_F32U(x).u;

    u32 abstop = abstop12_(x);
    if (abstop >= abstop12_(128.0f))
    {
        if (xu == MATS_F32U(-F32_INF).u) {
            return 0.0f;
        }
        if (abstop >= abstop12_(F32_INF)) {
            return x + x;
        }
        if (x > 0.0f) {
            return mats_overflow32(0);
        }
        if (x <= -150.0f) {
            return mats_underflow32(0);
        }
    }

    f64 xd = x;
    f64 kd = xd + gExp2F32_ShiftScaled; /* Rounding to double precision is required.  */
    u64 ki = MATS_F64U(kd).u;
    kd -= gExp2F32_ShiftScaled;

    f64 r = xd - kd;
    /* exp(x) = 2^(k/N) * 2^(r/N) ~= s * (C0*r^3 + C1*r^2 + C2*r + 1) */
    u64 t = gExp2F32_Table[ki % (1 << EXP2_32_TABLE_BITS)];
    t += ki << (52 - EXP2_32_TABLE_BITS);

    f64 s = MATS_F64U(t).f;
    f64 z = gExp2F32_Poly[0] * r + gExp2F32_Poly[1];
    f64 r2 = r * r;
    f64 y = gExp2F32_Poly[2] * r + 1.0;
    y = z * r2 + y;
    y = y * s;
    return (f32)y;
}

internal f32
pow2_32(f32 x)
{
    return exp2_32(x);
}

internal f32
log32(f32 x)
{
    u32 xu = MATS_F32U(x).u;

    if ((xu - 0x00800000) >= (MATS_F32_EXP_MASK - 0x00800000))
    {
        /* x < 0x1p-126 or inf or nan.  */
        if (xu * 2 == 0) {
            return mats_divzero32(1);
        }
        if (xu == MATS_F32_EXP_MASK) {
            return x;
        }
        if ((xu & MATS_F32_SIGN_MASK) || ((xu * 2) >= 0xFF000000)) {
            return mats_invalid32(x);
        }
        // NOTE(michiel): Normalize a subnormal
        xu = MATS_F32U(x * 0x1p23f).u;
        xu -= MATS_F32_EXP_SHIFT << MATS_F32_EXP_SHIFT;
    }

    /* x = 2^k z; where z is in range [OFF,2*OFF] and exact.
            The range is split into N subintervals.
            The ith subinterval contains z and c is near its center.  */
    u32 temp = xu - 0x3F330000;
    u32 index = (temp >> (MATS_F32_EXP_SHIFT - LOG32_TABLE_BITS)) % (1 << LOG32_TABLE_BITS);
    s32 k = (s32)temp >> MATS_F32_EXP_SHIFT;
    u32 iz = xu - (temp & (0x1FF << MATS_F32_EXP_SHIFT));
    f64 invC = gLogF32_Table[index].invC;
    f64 logC = gLogF32_Table[index].logC;
    f64 z = (f64)MATS_F32U(iz).f;

    /* log(x) = log1p(z/c-1) + log(c) + k*Ln2 */
    f64 r = z * invC - 1.0;
    f64 y0 = logC + (f64)k * gLogF32_Ln2;

    f64 r2 = r * r;
    f64 result = gLogF32_Poly[1] * r + gLogF32_Poly[2];
    result = gLogF32_Poly[0] * r2 + result;
    result = result * r2 + (y0 + r);
    return (f32)result;
}

internal f32
log2_32(f32 x)
{
    u32 xu = MATS_F32U(x).u;

    if ((xu - 0x00800000) >= (MATS_F32_EXP_MASK - 0x00800000))
    {
        /* x < 0x1p-126 or inf or nan.  */
        if (xu * 2 == 0) {
            return mats_divzero32(1);
        }
        if (xu == MATS_F32_EXP_MASK) {
            return x;
        }
        if ((xu & MATS_F32_SIGN_MASK) || ((xu * 2) >= 0xFF000000)) {
            return mats_invalid32(x);
        }
        // NOTE(michiel): Normalize a subnormal
        xu = MATS_F32U(x * 0x1p23f).u;
        xu -= MATS_F32_EXP_SHIFT << MATS_F32_EXP_SHIFT;
    }

    /* x = 2^k z; where z is in range [OFF,2*OFF] and exact.
            The range is split into N subintervals.
            The ith subinterval contains z and c is near its center.  */
    u32 temp = xu - 0x3F330000;
    u32 index = (temp >> (MATS_F32_EXP_SHIFT - LOG2_32_TABLE_BITS)) % (1 << LOG2_32_TABLE_BITS);
    u32 top = temp & 0xFF800000;
    u32 iz = xu - top;
    s32 k = (s32)temp >> MATS_F32_EXP_SHIFT;
    f64 invC = gLog2F32_Table[index].invC;
    f64 logC = gLog2F32_Table[index].logC;
    f64 z = (f64)MATS_F32U(iz).f;

    /* log(x) = log1p(z/c-1) + log(c) + k*Ln2 */
    f64 r = z * invC - 1.0;
    f64 y0 = logC + (f64)k;

    f64 r2 = r * r;
    f64 result = gLog2F32_Poly[1] * r + gLog2F32_Poly[2];
    result = gLog2F32_Poly[0] * r2 + result;
    f64 p = gLog2F32_Poly[3] * r + y0;
    result = result * r2 + p;
    return (f32)result;
}

internal f32
log10_32(f32 x)
{
    s32 hx = MATS_S32_FROM_F32(x);

    s32 k = 0;
    if (MATS_F32_UWORD_IS_ZERO(hx & MATS_F32_ABS_MASK)) {
        return -F32_INF; // -g2pow25F32 / 0.0f; // NOTE(michiel): log10(+/-0) = -inf
    } else if (hx < 0) {
        return F32_NAN;     // NOTE(michiel): log10(-X) = NaN
    } else if (!MATS_F32_UWORD_IS_FINITE(hx)) {
        return x + x;
    } else if (MATS_F32_UWORD_IS_SUBNORMAL(hx)) {
        k -= 25;
        x *= g2pow25F32; // NOTE(michiel): x * 2^25
        hx = MATS_S32_FROM_F32(x);
    }

    k += (hx >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;
    s32 i = (u32)k >> 31;
    hx = (hx & MATS_F32_MANT_MASK) | ((MATS_F32_EXP_BIAS - i) << MATS_F32_EXP_SHIFT);
    f32 y = (f32)(k + i);
    x = MATS_F32_FROM_S32(hx);
    f32 z = y * gLog10F32_2_lo + gInvLn10F32 * log32(x);
    return z + y * gLog10F32_2_hi;
}

internal f32
expm1_32(f32 x)
{
    // NOTE(michiel): result = exp(x) - 1
    u32 hx = MATS_F32U(x).u;
    f32 ln2_hi = gLn2HighF32s[0];
    f32 ln2_lo = gLn2LowF32s[0];

    s32 xSign = hx & MATS_F32_SIGN_MASK;

    hx &= MATS_F32_ABS_MASK;

    // NOTE(michiel): Filter out huge and non-finite values
    if (hx >= 0x4195B844)
    {
        // NOTE(michiel): |x| >= 27 * ln(2)
        if (MATS_F32_UWORD_IS_NAN(hx)) {
            return x + x;
        }
        if (MATS_F32_UWORD_IS_INFINITE(hx)) {
            return xSign ? -1.0f : x; // NOTE(michiel): exp(+inf) = inf - 1, exp(-inf) = 0 - 1
        }
        if (!xSign && (hx > MATS_F32_UWORD_LOG_MAX)) {
            return mats_overflow32(0);
        }
        if (xSign) {
            // NOTE(michiel): This did raise an inexact signal
            return -1.0f;
        }
    }

    // NOTE(michiel): Argument reduction
    f32 c = 0.0f;
    s32 k = 0;
    if (hx > 0x3EB17218)
    {
        f32 hi, lo;
        // NOTE(michiel): |x| > 0.5*ln(2)
        if (hx < 0x3F851592)
        {
            // NOTE(michiel): |x| < 1.5*ln(2)
            if (xSign) {
                hi = x + ln2_hi;
                lo = -ln2_lo;
                k  = -1;
            } else {
                hi = x - ln2_hi;
                lo = ln2_lo;
                k  = 1;
            }
        }
        else
        {
            k = (s32)(gInvLn2F32 * x + (xSign ? -0.5f : 0.5f));
            f32 t = (f32)k;
            hi = x - t * ln2_hi;
            lo = t * ln2_lo;
        }
        x = hi - lo;
        c = (hi - x) - lo;
    }
    else if (hx < 0x33000000)
    {
        // NOTE(michiel): |x| < 2^-25
        // NOTE(michiel): This did raise an inexact signal
        return x;
    }

    // NOTE(michiel): x is in primary range
    f32 hfx = 0.5f * x;
    f32 hxs = x * hfx;
#if 1
    f32 r1 = 1.5807170421e-3f * hxs;
    r1 = (r1 - 3.3333212137e-2f) * hxs;
#else
    f32 r1 = -2.0109921195e-07f * hxs;
    r1 = (r1 + 4.0082177293e-06f) * hxs;
    r1 = (r1 - 7.9365076090e-05f) * hxs;
    r1 = (r1 + 1.5873016091e-03f) * hxs;
    r1 = (r1 - 3.3333335072e-02f) * hxs;
#endif
    r1 = (r1 + 1.0f);
    f32 t = 3.0f - r1 * hfx;
    f32 e = hxs * ((r1 - t) / (6.0f - x * t));

    f32 result;
    if (k == 0)
    {
        result = x - (x * e - hxs); // NOTE(michiel): c == 0
    }
    else
    {
        e = (x * (e - c) - c);
        e -= hxs;
        if (k == -1)
        {
            result = 0.5f * (x - e) - 0.5f;
        }
        else if (k == 1)
        {
            if (x < -0.25f) {
                result = -2.0f * (e - (x + 0.5f));
            } else {
                result = 1.0f + 2.0f * (x - e);
            }
        }
        else
        {
            if ((k <= -2) || (k > 56))
            {
                f32 y = 1.0f - (e - x);
                u32 i = MATS_F32U(y).u;
                y = MATS_F32U(i + (k << MATS_F32_EXP_SHIFT)).f;
                result = y - 1.0f;
            }
            else if (k < MATS_F32_EXP_SHIFT)
            {
                f32 y = MATS_F32U((u32)0x3F800000 - (0x01000000 >> k)).f;
                y = y - (e - x);
                u32 i = MATS_F32U(y).u;
                result = MATS_F32U(i + (k << MATS_F32_EXP_SHIFT)).f;
            }
            else
            {
                f32 y = MATS_F32U((u32)(MATS_F32_EXP_BIAS - k) << MATS_F32_EXP_SHIFT).f;
                y = x - (e + y);
                y += 1.0f;
                u32 i = MATS_F32U(y).u;
                result = MATS_F32U(i + (k << MATS_F32_EXP_SHIFT)).f;
            }
        }
    }
    return result;
}

internal f32
log1p32(f32 x)
{
    // NOTE(michiel): result = log(1 + x)
    s32 hx = (s32)MATS_F32U(x).u;
    f32 ln2_hi = gLn2HighF32s[0];
    f32 ln2_lo = gLn2LowF32s[0];

    s32 ax = hx & MATS_F32_ABS_MASK;

    s32 k = 1;
    if (!MATS_F32_UWORD_IS_FINITE(hx)) {
        return x + x;
    }

    s32 hu = 0;
    f32 f = 0;
    if (hx < 0x3ED413D7)
    {
        // NOTE(michiel): x < 0.41422
        if (ax >= 0x3F800000)
        {
            // NOTE(michiel): x <= -1.0
            if (x == -1.0f) {
                return mats_divzero32(1); // NOTE(michiel): log1p(-1) = -inf
            } else {
                return mats_invalid32(x); // NOTE(michiel): log1p(x < -1) = NaN
            }
        }
        if (ax < 0x31000000)
        {
            // NOTE(michiel): |x| < 2^-29
            if (ax < 0x24800000) {
                // NOTE(michiel): |x| < 2^-54
                return x;
            } else {
                return x - x * x * 0.5f;
            }
        }
        if ((hx > 0) || (hx < (s32)0xBE95F61F))
        {
            // NOTE(michiel): -0.2929 < x < 0.41422
            k = 0;
            f = x;
            hu = 1;
        }
    }

    f32 c = 0;
    if (k != 0)
    {
        if (hx < 0x5A000000)
        {
            f32 u = 1.0f + x;
            hu = (s32)MATS_F32U(u).u;
            k  = (hu >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;
            // NOTE(michiel): Correction term
            c = (k > 0) ? 1.0f - (u - x) : x - (u - 1.0f);
            c /= u;
        }
        else
        {
            f32 u = x;
            hu = (s32)MATS_F32U(u).u;
            k  = (hu >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;
        }
        hu &= MATS_F32_MANT_MASK;
        f32 u;
        if (hu < 0x003504F7)
        {
            u = MATS_F32U((u32)hu | 0x3F800000).f; // NOTE(michiel): Normalize
        }
        else
        {
            ++k;
            u = MATS_F32U((u32)hu | 0x3F000000).f;
            hu = (0x00800000 - hu) >> 2;
        }
        f = u - 1.0f;
    }

    f32 hfSq = 0.5f * f * f;
    if (hu == 0)
    {
        // NOTE(michiel): |f| < 2^-20
        if (f == 0.0f)
        {
            if (k == 0) {
                return 0.0f;
            } else {
                c += (f32)k * ln2_lo;
                return (f32)k * ln2_hi + c;
            }
        }
        f32 r = hfSq * (1.0f - 0.66666666666666666f * f);
        if (k == 0) {
            return f - r;
        } else {
            return (f32)k * ln2_hi - ((r - ((f32)k * ln2_lo + c)) - f);
        }
    }
    f32 s = f / (2.0f + f);
    f32 z = s * s;
    f32 r = 1.4798198640e-01f * z;
    r = (r + 1.5313838422e-01f) * z;
    r = (r + 1.8183572590e-01f) * z;
    r = (r + 2.2222198546e-01f) * z;
    r = (r + 2.8571429849e-01f) * z;
    r = (r + 4.0000000596e-01f) * z;
    r = (r + 6.6666668653e-01f) * z;

    if (k == 0) {
        return f - (hfSq - s * (hfSq + r));
    } else {
        return (f32)k * ln2_hi - ((hfSq - (s * (hfSq + r) + ((f32)k * ln2_lo + c))) - f);
    }
}

internal f32
log1p_fast32(f32 x)
{
    // NOTE(michiel): Musl libc version
    /* |(log(1+s)-log(1-s))/s - Lg(s)| < 2**-34.24 (~[-4.95e-11, 4.97e-11]). */
    f32 ln2_hi = gLn2HighF32s[0];
    f32 ln2_lo = gLn2LowF32s[0];
    u32 ix = MATS_F32U(x).u;

    f32 c;
    f32 f;
    s32 k = 1;
    if ((ix < 0x3ED413D0) || (ix >> 31))
    {
        // NOTE(michiel): 1 + x < sqrt(2)
        if (ix >= 0xBF800000)
        {
            // NOTE(michiel): x <= -1
            if (x == -1.0f) {
                return x / 0.0f;
            } else {
                return (x - x) / 0.0f;
            }
        }
        if ((ix << 1) < (0x33800000 << 1))
        {
            // NOTE(michiel): |x| < 2^-24
            return x;
        }
        if (ix <= 0xBE95F619)
        {
            // NOTE(michiel): sqrt(2)/2 <= 1 + x < sqrt(2)
            k = 0;
            c = 0;
            f = x;
        }
    }
    else if (ix >= MATS_F32_EXP_MASK)
    {
        return x;
    }

    if (k)
    {
        mats_f32u u = MATS_F32U(1.0f + x);
        u32 iu = u.u;
        iu += 0x3F800000 - 0x3F3504F3;
        k = (s32)(iu >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;
        // NOTE(michiel): Correction term
        if (k < 25)
        {
            c = k >= 2 ? 1.0f - (u.f - x) : x - (u.f - 1.0f);
            c /= u.f;
        }
        else
        {
            c = 0.0f;
        }

        // NOTE(michiel): Reduce u into [sqrt(2)/2, sqrt(2)]
        iu = (iu & MATS_F32_MANT_MASK) + 0x3F3504F3;
        u.u = iu;
        f = u.f - 1.0f;
    }

    f32 s = f / (2.0f + f);
    f32 z = s * s;
    f32 w = z * z;
    f32 t1 = w * (0xccce13.0p-25f + w * 0xf89e26.0p-26f);
    f32 t2 = z * (0xaaaaaa.0p-24f + w * 0x91e9ee.0p-25f);
    f32 r  = t1 + t2;
    f32 hfSq = 0.5f * f * f;
    f32 dk = (f32)k;
    return s * (hfSq + r) + (dk * ln2_lo + c) - hfSq + f + dk * ln2_hi;
}

//
// NOTE(michiel): 64-bit
//

internal f64
sqrt64(f64 x)
{
    f64 result = 0.0;
#if MATS_USE_SSE2
    WideMath m;
    m.md = _mm_sqrt_sd(_mm_set_sd(x), _mm_set_sd(x));
    result = m.ed[0];
#else
    INVALID_CODE_PATH;
#endif
    return result;
}

internal f64
hypot64(f64 x, f64 y)
{
    s64 ax = MATS_S64_FROM_F64(x) & MATS_F64_ABS_MASK;
    s64 ay = MATS_S64_FROM_F64(y) & MATS_F64_ABS_MASK;

    if (ay > ax)
    {
        s64 temp = ax;
        ax = ay;
        ay = temp;
    }

    f64 a = MATS_F64_FROM_S64(ax);
    f64 b = MATS_F64_FROM_S64(ay);

    if ((ax - ay) > 0x3C00000000000000LL)
    {
        // NOTE(michiel): x / y > 2^60
        return a + b;
    }

    s32 k = 0;
    if (ax > 0x5F300000FFFFFFFFLL)
    {
        // NOTE(michiel): a > 2^500
        if (ax >= 0x7FF0000000000000LL) {
            // NOTE(michiel): Infinity or NaN
            f64 w = a + b;
            if ((ax & 0x000FFFFFFFFFFFFFLL) == 0) {
                w = a;
            }
            if ((ay ^ 0x7FF0000000000000LL) == 0) {
                w = b;
            }
            return w;
        }
        // NOTE(michiel): Scale a and b by 2^-600
        ax -= 0x2580000000000000LL;
        ay -= 0x2580000000000000LL;
        k  += 600;

        a = MATS_F64_FROM_S64(ax);
        b = MATS_F64_FROM_S64(ay);
    }

    if (ay < 0x20B0000000000000LL)
    {
        // NOTE(michiel): b < 2^-500
        if (ay < 0x0010000000000000LL)
        {
            // NOTE(michiel): Subnormal or zero
            if (ay == 0) {
                return a;
            }
            f64 mult = MATS_F64_FROM_U64(0x7FD0000000000000ULL); // NOTE(michiel): 2^1022
            b *= mult;
            a *= mult;
            k -= 1022;
        }
        else
        {
            // NOTE(michiel): Scale a and b by 2^600
            ax += 0x2580000000000000LL; // NOTE(michiel): a *= 2^600
            ay += 0x2580000000000000LL;
            k  -= 600;

            a = MATS_F64_FROM_S64(ax);
            b = MATS_F64_FROM_S64(ay);
        }
    }

    // NOTE(michiel): Medium sized a and b
    f64 w = a - b;
    if (w > b)
    {
        f64 t1 = MATS_F64_FROM_U64(ax & 0xFFFFFFFF00000000LL);
        f64 t2 = a - t1;
        w = sqrt64(t1 * t1 - (b * (-b) - t2 * (a + t1)));
    }
    else
    {
        a = a + a;
        f64 y1 = MATS_F64_FROM_U64(ay & 0xFFFFFFFF00000000LL);
        f64 y2 = b - y1;
        f64 t1 = MATS_F64_FROM_U64((ax & 0xFFFFFFFF00000000LL) + 0x0010000000000000LL);
        f64 t2 = a - t1;
        w = sqrt64(t1 * y1 - (w * (-w) - (t1 * y2 + t2 * b)));
    }

    if (k)
    {
        f64 t1 = MATS_F64_FROM_S64(0x3FF0000000000000LL + ((s64)k << 52));
        return t1 * w;
    }
    else
    {
        return w;
    }
}


internal f64
exp64(f64 x)
{
#define EXP64_N  (1 << EXP64_TABLE_BITS)
    u32 abstop = mats_top12(x) & 0x7FF;
    if ((abstop - mats_top12(0x1p-54)) >= (mats_top12(512.0) - mats_top12(0x1p-54)))
    {
        if ((abstop - mats_top12(0x1p-54)) >= 0x80000000)
        {
            return 1.0 + x;
        }
        if (abstop >= mats_top12(1024.0))
        {
            if (MATS_U64_FROM_F64(x) == MATS_U64_FROM_F64(-F64_INF)) {
                return 0.0;
            }
            if (abstop >= mats_top12(F64_INF)) {
                return 1.0 + x;
            }
            if (MATS_U64_FROM_F64(x) >> 63) {
                return mats_underflow64(0);
            } else {
                return mats_overflow64(0);
            }
        }
        abstop = 0;
    }
    // NOTE(michiel): exp(x) = 2^(k/N) * exp(r), with exp(r) in [2^(-1/2N),2^(1/2N)].
    // x = ln2/N*k + r, with int k and r in [-ln2/2N, ln2/2N].
    f64 z = gExp64Data.invln2N * x;
#if 0 && TOINT_INTRINSICS
    f64 kd = roundtoint(z);
    u64 ki = converttoint(z);
#else
#if EXP64_USE_TOINT_NARROW
    // NOTE(michiel): z - kd is in [-0.5-2^-16, 0.5] in all rounding modes.
    f64 kd = z + gExp64Data.shift;
    u64 ki = MATS_U64_FROM_F64(kd) >> 16;
    kd = (f64)(s32)ki;
#else
    // NOTE(michiel): z - kd is in [-1, 1] in non-nearest rounding modes.
    f64 kd = z + gExp64Data.shift;
    u64 ki = MATS_U64_FROM_F64(kd);
    kd -= gExp64Data.shift;
#endif
#endif
    f64 r = x + kd * gExp64Data.negln2hiN + kd * gExp64Data.negln2loN;
    // NOTE(michiel): 2^(k/N) ~= scale * (1 + tail).
    u64 idx = 2 * (ki % EXP64_N);
    u64 top = ki << (52 - EXP64_TABLE_BITS);
    f64 tail = MATS_F64_FROM_U64(gExp64Data.tab[idx]);
    // NOTE(michiel): This is only a valid scale when -1023*N < k < 1024*N.
    u64 sbits = gExp64Data.tab[idx + 1] + top;
    // NOTE(michiel): exp(x) = 2^(k/N) * exp(r) ~= scale + scale * (tail + exp(r) - 1).
    // Evaluation is optimized assuming superscalar pipelined execution.
    f64 r2 = r * r;
    // NOTE(michiel): Without fma the worst case error is 0.25/N ulp larger.
    // Worst case error is less than 0.5+1.11/N+(abs poly error * 2^53) ulp.
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
    if (abstop == 0)
    {
        f64 scale;
        f64 y;

        if ((ki & 0x80000000) == 0)
        {
            // NOTE(michiel): k > 0, the exponent of scale might have overflowed by <= 460
            sbits -= 1009ull << 52;
            scale = MATS_F64_FROM_U64(sbits);
            y = 0x1p1009 * (scale + scale * tmp);
            return y;
        }

        sbits += 1022ull << 52;
        scale = MATS_F64_FROM_U64(sbits);
        y = scale + scale * tmp;
        if (y < 1.0)
        {
            // NOTE(michiel): Round y to the right precision before scaling it into the subnormal
            // range to avoid double rounding that can cause 0.5+E/2 ulp error where
            // E is the worst-case ulp error outside the subnormal range.  So this
            // is only useful if the goal is better than 1 ulp worst-case error.
            f64 lo = scale - y + scale * tmp;
            f64 hi = 1.0 + y;
            lo = 1.0 - hi + y + lo;
            y = (hi + lo) - 1.0;
        }
        y = 0x1p-1022 * y;
        return y;
    }
    f64 scale = MATS_F64_FROM_U64(sbits);
    // NOTE(michiel): tmp == 0 or |tmp| > 2^-65 and scale > 2^-739, so there
    // is no spurious underflow here even without fma.
    return scale + scale * tmp;
#undef EXP64_N
}

internal f64
exp2_64(f64 x)
{
#define EXP2_64_N  (1 << EXP64_TABLE_BITS)
    u32 abstop = mats_top12(x) & 0x7FF;
    if ((abstop - mats_top12(0x1p-54)) >= (mats_top12(512.0) - mats_top12(0x1p-54)))
    {
        if ((abstop - mats_top12(0x1p-54)) >= 0x80000000) {
            /* Avoid spurious underflow for tiny x.  */
            /* Note: 0 is common input.  */
            return 1.0 + x;
        }
        if (abstop >= mats_top12(1024.0))
        {
            if (MATS_U64_FROM_F64(x) == MATS_U64_FROM_F64(-F64_INF)) {
                return 0.0;
            }
            if (abstop >= mats_top12(F64_INF)) {
                return 1.0 + x;
            }
            if (!(MATS_U64_FROM_F64(x) >> 63)) {
                return mats_overflow64(0);
            } else if (MATS_U64_FROM_F64(x) >= MATS_U64_FROM_F64(-1075.0)) {
                return mats_underflow64(0);
            }
        }
        if ((2 * MATS_U64_FROM_F64(x)) > (2 * MATS_U64_FROM_F64(928.0))) {
            /* Large x is special cased below.  */
            abstop = 0;
        }
    }

    /* exp2(x) = 2^(k/N) * 2^r, with 2^r in [2^(-1/2N),2^(1/2N)].  */
    /* x = k/N + r, with int k and r in [-1/2N, 1/2N].  */
    f64 kd = x + gExp64Data.exp2_shift;
    u64 ki = MATS_U64_FROM_F64(kd); /* k.  */
    kd -= gExp64Data.exp2_shift; /* k/N for int k.  */
    f64 r = x - kd;
    /* 2^(k/N) ~= scale * (1 + tail).  */
    u64 idx = 2 * (ki % EXP2_64_N);
    u64 top = ki << (52 - EXP64_TABLE_BITS);
    f64 tail = MATS_F64_FROM_U64(gExp64Data.tab[idx]);
    /* This is only a valid scale when -1023*N < k < 1024*N.  */
    u64 sbits = gExp64Data.tab[idx + 1] + top;
    /* exp2(x) = 2^(k/N) * 2^r ~= scale + scale * (tail + 2^r - 1).  */
    /* Evaluation is optimized assuming superscalar pipelined execution.  */
    f64 r2 = r * r;
    /* Without fma the worst case error is 0.5/N ulp larger.  */
    /* Worst case error is less than 0.5+0.86/N+(abs poly error * 2^53) ulp.  */

#define C1 gExp64Data.exp2_poly[0]
#define C2 gExp64Data.exp2_poly[1]
#define C3 gExp64Data.exp2_poly[2]
#define C4 gExp64Data.exp2_poly[3]
#define C5 gExp64Data.exp2_poly[4]
#define C6 gExp64Data.exp2_poly[5]

#if EXP2_64_POLY_ORDER == 4
    f64 tmp = tail + r * C1 + r2 * C2 + r * r2 * (C3 + r * C4);
#elif EXP2_64_POLY_ORDER == 5
    f64 tmp = tail + r * C1 + r2 * (C2 + r * C3) + r2 * r2 * (C4 + r * C5);
#elif EXP2_64_POLY_ORDER == 6
    f64 tmp = tail + r * C1 + r2 * (0.5 + r * C3) + r2 * r2 * (C4 + r * C5 + r2 * C6);
#endif

#undef C1
#undef C2
#undef C3
#undef C4
#undef C5
#undef C6

    if (abstop == 0)
    {
        /* Handle cases that may overflow or underflow when computing the result that
           is scale*(1+TMP) without intermediate rounding.  The bit representation of
           scale is in SBITS, however it has a computed exponent that may have
           overflown into the sign bit so that needs to be adjusted before using it as
           a double.  (int32_t)KI is the k used in the argument reduction and exponent
           adjustment of scale, positive k here means the result may overflow and
           negative k means the result may underflow.  */
        f64 scale, y;

        if ((ki & 0x80000000) == 0)
        {
            /* k > 0, the exponent of scale might have overflowed by 1.  */
            sbits -= 1ull << 52;
            scale = MATS_F64_FROM_U64(sbits);
            y = 2 * (scale + scale * tmp);
            return y;
        }
        /* k < 0, need special care in the subnormal range.  */
        sbits += 1022ull << 52;
        scale = MATS_F64_FROM_U64(sbits);
        y = scale + scale * tmp;
        if (y < 1.0)
        {
            /* Round y to the right precision before scaling it into the subnormal
           range to avoid double rounding that can cause 0.5+E/2 ulp error where
           E is the worst-case ulp error outside the subnormal range.  So this
           is only useful if the goal is better than 1 ulp worst-case error.  */
            f64 hi, lo;
            lo = scale - y + scale * tmp;
            hi = 1.0 + y;
            lo = 1.0 - hi + y + lo;
            y = (hi + lo) - 1.0;
        }
        y = 0x1p-1022 * y;
        return y;
    }
    f64 scale = MATS_F64_FROM_U64(sbits);
    /* Note: tmp == 0 or |tmp| > 2^-65 and scale > 2^-928, so there
       is no spurious underflow here even without fma.  */
    return scale + scale * tmp;
#undef EXP2_64_N
}

internal f64
log64(f64 x)
{
#define LOG64_N (1 << LOG64_TABLE_BITS)

    u64 ix = MATS_U64_FROM_F64(x);
    u32 top = mats_top16(x);

#if (LOG64_POLY1_ORDER == 10) || (LOG64_POLY1_ORDER == 11)
#define LOG64_LO   MATS_U64_FROM_F64(1.0 - 0x1p-5)
#define LOG64_HI   MATS_U64_FROM_F64(1.0 + 0x1.1p-5)
#elif LOG64_POLY1_ORDER == 12
#define LOG64_LO   MATS_U64_FROM_F64(1.0 - 0x1p-4)
#define LOG64_HI   MATS_U64_FROM_F64(1.0 + 0x1.09p-4)
#endif
    if ((ix - LOG64_LO) < (LOG64_HI - LOG64_LO))
    {
        // NOTE(michiel): Handle close to 1.0 inputs separately.
        // Fix sign of zero with downward rounding when x==1.
        if (ix == MATS_U64_FROM_F64(1.0)) {
            return 0;
        }
        f64 r = x - 1.0;
        f64 r2 = r * r;
        f64 r3 = r * r2;
#if LOG64_POLY1_ORDER == 10
        /* Worst-case error is around 0.516 ULP.  */
        f64 y = r3 * (gLogData64.poly1[1] + r * gLogData64.poly1[2] + r2 * gLogData64.poly1[3]
                      + r3 * (gLogData64.poly1[4] + r * gLogData64.poly1[5] + r2 * gLogData64.poly1[6] + r3 * (gLogData64.poly1[7] + r * gLogData64.poly1[8])));
        f64 w = gLogData64.poly1[0] * r2; /* gLogData64.poly1[0] == -0.5.  */
        f64 hi = r + w;
        y += r - hi + w;
        y += hi;
#elif LOG64_POLY1_ORDER == 11
        /* Worst-case error is around 0.516 ULP.  */
        f64 y = r3 * (gLogData64.poly1[1] + r * gLogData64.poly1[2]
                      + r2 * (gLogData64.poly1[3] + r * gLogData64.poly1[4] + r2 * gLogData64.poly1[5]
                              + r3 * (gLogData64.poly1[6] + r * gLogData64.poly1[7] + r2 * gLogData64.poly1[8] + r3 * gLogData64.poly1[9])));
        f64 w = gLogData64.poly1[0] * r2; /* gLogData64.poly1[0] == -0.5.  */
        f64 hi = r + w;
        y += r - hi + w;
        y += hi;
#elif LOG64_POLY1_ORDER == 12
        f64 y = r3 * (gLogData64.poly1[1] + r * gLogData64.poly1[2] + r2 * gLogData64.poly1[3]
                      + r3 * (gLogData64.poly1[4] + r * gLogData64.poly1[5] + r2 * gLogData64.poly1[6]
                              + r3 * (gLogData64.poly1[7] + r * gLogData64.poly1[8] + r2 * gLogData64.poly1[9] + r3 * gLogData64.poly1[10])));
#if LOG64_N <= 64
        /* Worst-case error is around 0.532 ULP.  */
        f64 w = gLogData64.poly1[0] * r2; /* gLogData64.poly1[0] == -0.5.  */
        f64 hi = r + w;
        y += r - hi + w;
        y += hi;
#else
        /* Worst-case error is around 0.507 ULP.  */
        f64 w = r * 0x1p27;
        f64 rhi = r + w - w;
        f64 rlo = r - rhi;
        w = rhi * rhi * gLogData64.poly1[0]; /* gLogData64.poly1[0] == -0.5.  */
        f64 hi = r + w;
        f64 lo = r - hi + w;
        lo += gLogData64.poly1[0] * rlo * (rhi + r);
        y += lo;
        y += hi;
#endif
#endif
        return y;
    }
    if ((top - 0x0010) >= (0x7ff0 - 0x0010))
    {
        /* x < 0x1p-1022 or inf or nan.  */
        if (ix * 2 == 0) {
            return mats_divzero64(1);
        }
        if (ix == MATS_U64_FROM_F64(F64_INF)) { /* log(inf) == inf.  */
            return x;
        }
        if ((top & 0x8000) || ((top & 0x7ff0) == 0x7ff0)) {
            return mats_invalid64(x);
        }
        /* x is subnormal, normalize it.  */
        ix = MATS_U64_FROM_F64(x * 0x1p52);
        ix -= 52ULL << 52;
    }

    /* x = 2^k z; where z is in range [OFF,2*OFF) and exact. (OFF = 0x3fe6000000000000)
       The range is split into N subintervals.
       The ith subinterval contains z and c is near its center.  */
    u64 tmp = ix - 0x3fe6000000000000ULL;
    s32 i = (tmp >> (52 - LOG64_TABLE_BITS)) % LOG64_N;
    s32 k = (s64)tmp >> 52; /* arithmetic shift */
    u64 iz = ix - (tmp & (0xfffULL << 52));
    f64 invc = gLogData64.tab[i].invc;
    f64 logc = gLogData64.tab[i].logc;
    f64 z = MATS_F64_FROM_U64(iz);

    /* log(x) = log1p(z/c-1) + log(c) + k*Ln2.  */
    /* r ~= z/c - 1, |r| < 1/(2*N).  */
#if MATS_HAVE_FAST_FMA
    /* rounding error: 0x1p-55/N.  */
    f64 r = mats_fma(z, invc, -1.0);
#else
    /* rounding error: 0x1p-55/N + 0x1p-66.  */
    f64 r = (z - gLogData64.tab2[i].chi - gLogData64.tab2[i].clo) * invc;
#endif
    f64 kd = (f64) k;

    /* hi + lo = r + log(c) + k*Ln2.  */
    f64 w = kd * gLogData64.ln2hi+ logc;
    f64 hi = w + r;
    f64 lo = w - hi + r + kd * gLogData64.ln2lo;

    /* log(x) = lo + (log1p(r) - r) + hi.  */
    f64 r2 = r * r; /* rounding error: 0x1p-54/N^2.  */
    /* Worst case error if |y| > 0x1p-5:
       0.5 + 4.13/N + abs-poly-error*2^57 ULP (+ 0.002 ULP without fma)
       Worst case error if |y| > 0x1p-4:
       0.5 + 2.06/N + abs-poly-error*2^56 ULP (+ 0.001 ULP without fma).  */
#if LOG64_POLY_ORDER == 6
    f64 y = lo + r2 * gLogData64.poly[0] + r * r2 * (gLogData64.poly[1] + r * gLogData64.poly[2] + r2 * (gLogData64.poly[3] + r * gLogData64.poly[4])) + hi;
#elif LOG64_POLY_ORDER == 7
    f64 y = lo
        + r2 * (gLogData64.poly[0] + r * gLogData64.poly[1] + r2 * (gLogData64.poly[2] + r * gLogData64.poly[3])
                + r2 * r2 * (gLogData64.poly[4] + r * gLogData64.poly[5]))
        + hi;
#else
#error LOG64_POLY_ORDER invalid
#endif
    return y;
#undef LOG64_LO
#undef LOG64_HI
#undef LOG64_N
}

internal f64
log2_64(f64 x)
{
#define LOG2_64_N (1 << LOG2_64_TABLE_BITS)
    u64 ix = MATS_U64_FROM_F64(x);
    u32 top = mats_top16(x);

#define LOG2_64_LO MATS_U64_FROM_F64(1.0 - 0x1.5b51p-5)
#define LOG2_64_HI MATS_U64_FROM_F64(1.0 + 0x1.6ab2p-5)
    if ((ix - LOG2_64_LO) < (LOG2_64_HI - LOG2_64_LO))
    {
        /* Handle close to 1.0 inputs separately.  */
        /* Fix sign of zero with downward rounding when x==1.  */
        if (ix == MATS_U64_FROM_F64(1.0)) {
            return 0;
        }
        f64 r = x - 1.0;
#if MATS_HAVE_FAST_FMA
        f64 hi = r * gLog2Data64.invln2hi;
        f64 lo = r * gLog2Data64.invln2lo + fma (r, gLog2Data64.invln2hi, -hi);
#else
        f64 rhi = MATS_F64_FROM_U64(MATS_U64_FROM_F64(r) & -1ULL << 32);
        f64 rlo = r - rhi;
        f64 hi = rhi * gLog2Data64.invln2hi;
        f64 lo = rlo * gLog2Data64.invln2hi+ r * gLog2Data64.invln2lo;
#endif
        f64 r2 = r * r; /* rounding error: 0x1p-62.  */
        f64 r4 = r2 * r2;
        /* Worst-case error is less than 0.54 ULP (0.55 ULP without fma).  */
        f64 p = r2 * (gLog2Data64.poly1[0] + r * gLog2Data64.poly1[1]);
        f64 y = hi + p;
        lo += hi - y + p;
        lo += r4 * (gLog2Data64.poly1[2] + r * gLog2Data64.poly1[3] + r2 * (gLog2Data64.poly1[4] + r * gLog2Data64.poly1[5])
                    + r4 * (gLog2Data64.poly1[6] + r * gLog2Data64.poly1[7] + r2 * (gLog2Data64.poly1[8] + r * gLog2Data64.poly1[9])));
        y += lo;
        return y;
    }
    if ((top - 0x0010) >= (0x7ff0 - 0x0010))
    {
        /* x < 0x1p-1022 or inf or nan.  */
        if (ix * 2 == 0) {
            return mats_divzero64(1);
        }
        if (ix == MATS_U64_FROM_F64(F64_INF)) { /* log(inf) == inf.  */
            return x;
        }
        if ((top & 0x8000) || ((top & 0x7ff0) == 0x7ff0)) {
            return mats_invalid64(x);
        }
        /* x is subnormal, normalize it.  */
        ix = MATS_U64_FROM_F64(x * 0x1p52);
        ix -= 52ULL << 52;
    }

    /* x = 2^k z; where z is in range [OFF,2*OFF) and exact. (OFF = 0x3fe6000000000000)
       The range is split into N subintervals.
       The ith subinterval contains z and c is near its center.  */
    u64 tmp = ix - 0x3fe6000000000000;
    s32 i = (tmp >> (52 - LOG2_64_TABLE_BITS)) % LOG2_64_N;
    s32 k = (s64)tmp >> 52; /* arithmetic shift */
    u64 iz = ix - (tmp & (0xfffULL << 52));
    f64 invc = gLog2Data64.tab[i].invc;
    f64 logc = gLog2Data64.tab[i].logc;
    f64 z = MATS_F64_FROM_U64(iz);
    f64 kd = (f64)k;

    /* log2(x) = log2(z/c) + log2(c) + k.  */
    /* r ~= z/c - 1, |r| < 1/(2*N).  */
#if MATS_HAVE_FAST_FMA
    /* rounding error: 0x1p-55/N.  */
    f64 r = mats_fma(z, invc, -1.0);
    f64 t1 = r * gLog2Data64.invln2hi;
    f64 t2 = r * gLog2Data64.invln2lo + mats_fma(r, gLog2Data64.invln2hi, -t1);
#else
    /* rounding error: 0x1p-55/N + 0x1p-65.  */
    f64 r = (z - gLog2Data64.tab2[i].chi - gLog2Data64.tab2[i].clo) * invc;
    f64 rhi = MATS_F64_FROM_U64(MATS_U64_FROM_F64(r) & -1ULL << 32);
    f64 rlo = r - rhi;
    f64 t1 = rhi * gLog2Data64.invln2hi;
    f64 t2 = rlo * gLog2Data64.invln2hi + r * gLog2Data64.invln2lo;
#endif

    /* hi + lo = r/ln2 + log2(c) + k.  */
    f64 t3 = kd + logc;
    f64 hi = t3 + t1;
    f64 lo = t3 - hi + t1 + t2;

    /* log2(r+1) = r/ln2 + r^2*poly(r).  */
    /* Evaluation is optimized assuming superscalar pipelined execution.  */
    f64 r2 = r * r; /* rounding error: 0x1p-54/N^2.  */
    f64 r4 = r2 * r2;
    /* Worst-case error if |y| > 0x1p-4: 0.547 ULP (0.550 ULP without fma).
       ~ 0.5 + 2/N/ln2 + abs-poly-error*0x1p56 ULP (+ 0.003 ULP without fma).  */
    f64 p = gLog2Data64.poly[0] + r * gLog2Data64.poly[1] + r2 * (gLog2Data64.poly[2] + r * gLog2Data64.poly[3]) + r4 * (gLog2Data64.poly[4] + r * gLog2Data64.poly[5]);
    f64 y = lo + r2 * p + hi;

    return y;
#undef LOG2_64_LO
#undef LOG2_64_HI
#undef LOG2_64_N
}

internal f64
log10_64(f64 x)
{
    s64 sx = MATS_S64_FROM_F64(x);

    s32 k = 0;
    if (sx < 0x0010000000000000LL)
    {
        // NOTE(michiel): x < 2**-1022
        if ((sx & MATS_F64_ABS_MASK) == 0) {
            return -F32_INF; // -g2pow54F64 / 0.0;       /* log(+-0)=-inf */
        } else if (sx < 0) {
            return (x - x) / 0.0;           /* log(-#) = NaN */
        }
        k -= 54;
        x *= g2pow54F64; /* subnormal number, scale up x */
        sx = MATS_S64_FROM_F64(x);
    }
    if (sx >= 0x7FF0000000000000LL) {
        return x + x;
    }

    k += (sx >> MATS_F64_EXP_SHIFT) - MATS_F64_EXP_BIAS;
    s32 i  = (u32)k >> 31;
    sx = (sx & MATS_F64_MANT_MASK) | ((s64)(MATS_F64_EXP_BIAS - i) << MATS_F64_EXP_SHIFT);
    f64 y  = (f64)(k + i);
    x = MATS_F64_FROM_S64(sx);
    f64 z = y * gLog10F64_2_lo + gInvLn10F64 * log64(x);
    return  z + y * gLog10F64_2_hi;
}

internal f64
expm1_64(f64 x)
{
    // NOTE(michiel): result = exp(x) - 1
    f64 ln2_hi = gLn2HighF64s[0];
    f64 ln2_lo = gLn2LowF64s[0];

    s64 sx = MATS_S64_FROM_F64(x);
    u32 hx = sx >> 32;
    s32 xSign = hx & 0x80000000;   /* sign bit of x */
    hx &= 0x7FFFFFFF;            /* high word of |x| */
    u64 ux = sx & MATS_F64_ABS_MASK;

    /* filter out huge and non-finite argument */
	if (ux >= 0x4043687A00000000ULL)
    {
        // NOTE(michiel): |x| >= 56 * ln2
	    if (ux >= 0x40862E4200000000ULL)
        {
            // NOTE(michiel): |x| >= 709.78... */
            if (ux >= 0x7FF0000000000000ULL)
            {
                if ((sx & MATS_F64_MANT_MASK) != 0) {
                    return x + x; // NOTE(michiel): NaN
                } else {
                    return xSign ? -1.0 : x; // NOTE(michiel): exp(+inf) = inf - 1, exp(-inf) = 0 - 1
                }
	        }
	        if (x > 7.09782712893383973096e+02) {
                return mats_overflow64(0);
            }
	    }
	    if (xSign) {
            // NOTE(michiel): x < -56*ln2, return -1.0 with inexact
            // NOTE(michiel): This did raise an inexact signal
            return -1.0;
	    }
	}

    f64 c = 0.0;
    s32 k = 0;
    /* argument reduction */
	if (ux > 0x3FD62E42FFFFFFFFULL)
    {
        f64 hi, lo;
		// NOTE(michiel): |x| > 0.5 * ln2
	    if (ux < 0x3FF0A2B200000000ULL)
        {
            // NOTE(michiel): |x| < 1.5 * ln2
            if (xSign) {
                hi = x + ln2_hi;
                lo = -ln2_lo;
                k = -1;
            } else {
                hi = x - ln2_hi;
                lo =  ln2_lo;
                k =  1;
            }
	    }
        else
        {
            k  = (s32)(gInvLn2F64 * x + (xSign ? -0.5 : 0.5));
            f64 t  = (f64)k;
            hi = x - t * ln2_hi;	/* t * ln2_hi is exact here */
            lo = t * ln2_lo;
	    }
	    x  = hi - lo;
	    c  = (hi - x) - lo;
	}
	else if (ux < 0x3C90000000000000ULL)
    {
        // NOTE(michiel): |x|<2**-54, return x
	    return x;
	}

#define Q1  -3.33333333333331316428e-02  /* 0xBFA11111111110F4 */
#define Q2   1.58730158725481460165e-03  /* 0x3F5A01A019FE5585 */
#define Q3  -7.93650757867487942473e-05  /* 0xBF14CE199EAADBB7 */
#define Q4   4.00821782732936239552e-06  /* 0x3ED0CFCA86E65239 */
#define Q5  -2.01099218183624371326e-07  /* 0xBE8AFDB76E09C32D */

    /* x is now in primary range */
	f64 hfx = 0.5 * x;
	f64 hxs = x * hfx;
    f64 r1 = Q5 * hxs;
    r1 = (r1 + Q4) * hxs;
    r1 = (r1 + Q3) * hxs;
    r1 = (r1 + Q2) * hxs;
    r1 = (r1 + Q1) * hxs;
    r1 = r1 + 1.0;
	f64 t  = 3.0 - r1 * hfx;
	f64 e  = hxs * ((r1 - t) / (6.0 - x * t));

#undef Q1
#undef Q2
#undef Q3
#undef Q4
#undef Q5

    f64 result;
    if (k == 0)
    {
        result = x - (x * e - hxs); /* c is 0 */
    }
	else
    {
	    e  = (x * (e - c) - c);
	    e -= hxs;
	    if (k == -1)
        {
            result = 0.5 * (x - e) - 0.5;
        }
        else if (k == 1)
        {
            if (x < -0.25) {
                result = -2.0 * (e - (x + 0.5));
            } else {
                result = 1.0 + 2.0 * (x - e);
            }
        }
        else
        {
            if ((k <= -2) || (k > 56))
            {   /* suffice to return exp(x)-1 */
                f64 y = 1.0 - (e - x);
                u64 i = MATS_U64_FROM_F64(y);
                y = MATS_F64_FROM_U64(i + ((s64)k << MATS_F64_EXP_SHIFT)); /* add k to y's exponent */
                result = y - 1.0;
            }
            else if (k < 20)
            {
                f64 y = MATS_F64_FROM_U64(0x3FF0000000000000ULL - (0x0020000000000000ULL >> k)); // 1 - 2^-k
                y = y - (e - x);
                u64 i = MATS_U64_FROM_F64(y);
                result = MATS_F64_FROM_U64(i + ((s64)k << MATS_F64_EXP_SHIFT));
            }
            else
            {
                f64 y = MATS_F64_FROM_U64((u64)(MATS_F64_EXP_BIAS - k) << MATS_F64_EXP_SHIFT); // 2^-k
                y = x - (e + y);
                y += 1.0;
                u64 i = MATS_U64_FROM_F64(y);
                result = MATS_F64_FROM_U64(i + ((s64)k << MATS_F64_EXP_SHIFT));
            }
        }
    }
    return result;
}

internal f64
log1p64(f64 x)
{
    f64 ln2_hi = gLn2HighF64s[0];
    f64 ln2_lo = gLn2LowF64s[0];

	s64 sx = MATS_S64_FROM_F64(x);
    s32 hx = sx >> 32;
    s32 ax = hx & 0x7FFFFFFF;

    f64 f = 0.0;
    s32 hu = 0;
    s32 k = 1;
    if (hx < 0x3FDA827A)
    {
        // NOTE(michiel): x < 0.41422
        if (ax >= 0x3FF00000)
        {
            // NOTE(michiel): x <= -1.0
            if (x == -1.0) {
                return mats_divzero64(1); // log1p(-1) = -inf
            } else {
                return mats_invalid64(x); // log1p(x < -1) = NaN
            }
        }
        if (ax < 0x3E200000)
        {
            // NOTE(michiel): |x| < 2^-29
            if (((g2pow54F64 + x) > 0.0) &&
                (ax < 0x3C900000)) {
                // NOTE(michiel): |x| < 2^-54
                return x;
            } else {
                return x - x * x * 0.5;
            }
        }
        if ((hx > 0) || (hx <= ((s32)0xBFD2BEC3)))
        {
            k = 0;
            f = x;
            hu = 1;
        }
    }
    if (hx >= 0x7FF00000) {
        return x + x;
    }

    f64 c = 0;
	if (k != 0)
    {
        s64 shu;
	    if (hx < 0x43400000)
        {
            f64 u  = 1.0 + x;
            shu = MATS_S64_FROM_F64(u);
	        k  = (shu >> MATS_F64_EXP_SHIFT) - MATS_F64_EXP_BIAS;
	        c  = (k > 0) ? 1.0 - (u - x) : x - (u - 1.0); /* correction term */
            c /= u;
	    }
        else
        {
            f64 u = x;
            shu = MATS_S64_FROM_F64(u);
	        k  = (shu >> MATS_F64_EXP_SHIFT) - MATS_F64_EXP_BIAS;
	    }
	    hu = (shu >> 32) & 0x000fffff;
        f64 u;
	    if (hu < 0x0006A09E) {
            u = MATS_F64_FROM_S64(((s64)(hu | 0x3FF00000) << 32) | (shu & 0xFFFFFFFF)); /* normalize u */
	    } else {
	        ++k;
            u = MATS_F64_FROM_S64(((s64)(hu | 0x3FE00000) << 32) | (shu & 0xFFFFFFFF)); /* normalize u/2 */
	        hu = (0x00100000 - hu) >> 2;
	    }
	    f = u - 1.0;
	}

	f64 hfSq = 0.5 * f * f;
	if (hu == 0)
    {
        // NOTE(michiel): |f| < 2^-20
        if (f == 0.0)
        {
            if (k == 0) {
                return 0.0;
            } else {
                c += (f64)k * ln2_lo;
                return (f64)k * ln2_hi + c;
            }
        }
	    f64 r = hfSq * (1.0 - 0.66666666666666666 * f);
	    if (k == 0) {
            return f - r;
        } else {
            return k * ln2_hi - ((r - (k * ln2_lo + c)) - f);
        }
	}
    f64 s = f / (2.0 + f);
	f64 z = s * s;

#define Lp1 6.666666666666735130e-01  /* 0x3FE5555555555593 */
#define Lp2 3.999999999940941908e-01  /* 0x3FD999999997FA04 */
#define Lp3 2.857142874366239149e-01  /* 0x3FD2492494229359 */
#define Lp4 2.222219843214978396e-01  /* 0x3FCC71C51D8E78AF */
#define Lp5 1.818357216161805012e-01  /* 0x3FC7466496CB03DE */
#define Lp6 1.531383769920937332e-01  /* 0x3FC39A09D078C69F */
#define Lp7 1.479819860511658591e-01  /* 0x3FC2F112DF3E5244 */

    f64 r = Lp7 * z;
    r = (r + Lp6) * z;
    r = (r + Lp5) * z;
    r = (r + Lp4) * z;
    r = (r + Lp3) * z;
    r = (r + Lp2) * z;
    r = (r + Lp1) * z;

#undef Lp1
#undef Lp2
#undef Lp3
#undef Lp4
#undef Lp5
#undef Lp6
#undef Lp7

	if (k == 0) {
        return f - (hfSq - s * (hfSq + r));
    } else {
        return (f64)k * ln2_hi - ((hfSq - (s * (hfSq + r) + ((f64)k * ln2_lo + c))) - f);
    }
}
