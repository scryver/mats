//
// NOTE(michiel): Sqrt/Hypot/Exp/Exp2/Log/Log2/Log10
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
	s32 ix = (s32)u32f32(x).u;
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

    result = u32f32((u32)ix).f;
#endif
    return result;
}

internal f32
hypot32(f32 x, f32 y)
{
    u32 ux = u32f32(x).u & MATS_F32_ABS_MASK;
    u32 uy = u32f32(y).u & MATS_F32_ABS_MASK;

	if (ux < uy) {
        u32 temp = ux;
        ux = uy;
        uy = temp;
	}

    x = u32f32(ux).f;
    y = u32f32(uy).f;
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
    f32 result = z * sqrt32((f64)x * x + (f64)y * y);
    //f32 result = z * sqrt32(x * x + y * y); // NOTE(michiel): Less accurate, but faster
    return result;
}

internal f32
exp32(f32 x)
{
    u32 xu = u32f32(x).u;

    u32 abstop = abstop12_(x);
    if (abstop >= abstop12_(88.0f))
    {
        if (xu == u32f32(-F32_INF).u) {
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
    u64 ki = u64f64(kd).u;
    kd -= gExp2F32_Shift;

    f64 r = z - kd;

    /* exp(x) = 2^(k/N) * 2^(r/N) ~= s * (C0*r^3 + C1*r^2 + C2*r + 1) */
    u64 t = gExp2F32_Table[ki % (1 << EXP2F_TABLE_BITS)];
    t += ki << (52 - EXP2F_TABLE_BITS);

    f64 s = u64f64(t).f;
    f64 res = gExp2F32_PolyScaled[0] * r + gExp2F32_PolyScaled[1];
    f64 r2 = r * r;
    f64 y = gExp2F32_PolyScaled[2] * r + 1.0;
    res = res * r2 + y;
    res = res * s;

    return res;
}

internal f32
exp2_32(f32 x)
{
    // NOTE(michiel): result = 2^x
    u32 xu = u32f32(x).u;

    u32 abstop = abstop12_(x);
    if (abstop >= abstop12_(128.0f))
    {
        if (xu == u32f32(-F32_INF).u) {
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
    u64 ki = u64f64(kd).u;
    kd -= gExp2F32_ShiftScaled;

    f64 r = xd - kd;
    /* exp(x) = 2^(k/N) * 2^(r/N) ~= s * (C0*r^3 + C1*r^2 + C2*r + 1) */
    u64 t = gExp2F32_Table[ki % (1 << EXP2F_TABLE_BITS)];
    t += ki << (52 - EXP2F_TABLE_BITS);

    f64 s = u64f64(t).f;
    f64 z = gExp2F32_Poly[0] * r + gExp2F32_Poly[1];
    f64 r2 = r * r;
    f64 y = gExp2F32_Poly[2] * r + 1.0;
    y = z * r2 + y;
    y = y * s;
    return y;
}

internal f32
pow2_32(f32 x)
{
    return exp2_32(x);
}

internal f32
log32(f32 x)
{
    u32 xu = u32f32(x).u;

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
        xu = u32f32(x * 0x1p23f).u;
        xu -= MATS_F32_EXP_SHIFT << MATS_F32_EXP_SHIFT;
    }

    /* x = 2^k z; where z is in range [OFF,2*OFF] and exact.
            The range is split into N subintervals.
            The ith subinterval contains z and c is near its center.  */
    u32 temp = xu - 0x3F330000;
    u32 index = (temp >> (MATS_F32_EXP_SHIFT - LOGF_TABLE_BITS)) % (1 << LOGF_TABLE_BITS);
    s32 k = (s32)temp >> MATS_F32_EXP_SHIFT;
    u32 iz = xu - (temp & (0x1FF << MATS_F32_EXP_SHIFT));
    f64 invC = gLogF32_Table[index].invC;
    f64 logC = gLogF32_Table[index].logC;
    f64 z = (f64)u32f32(iz).f;

    /* log(x) = log1p(z/c-1) + log(c) + k*Ln2 */
    f64 r = z * invC - 1.0;
    f64 y0 = logC + (f64)k * gLogF32_Ln2;

    f64 r2 = r * r;
    f64 result = gLogF32_Poly[1] * r + gLogF32_Poly[2];
    result = gLogF32_Poly[0] * r2 + result;
    result = result * r2 + (y0 + r);
    return result;
}

internal f32
log2_32(f32 x)
{
    u32 xu = u32f32(x).u;

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
        xu = u32f32(x * 0x1p23f).u;
        xu -= MATS_F32_EXP_SHIFT << MATS_F32_EXP_SHIFT;
    }

    /* x = 2^k z; where z is in range [OFF,2*OFF] and exact.
            The range is split into N subintervals.
            The ith subinterval contains z and c is near its center.  */
    u32 temp = xu - 0x3F330000;
    u32 index = (temp >> (MATS_F32_EXP_SHIFT - LOG2F_TABLE_BITS)) % (1 << LOG2F_TABLE_BITS);
    u32 top = temp & 0xFF800000;
    u32 iz = xu - top;
    s32 k = (s32)temp >> MATS_F32_EXP_SHIFT;
    f64 invC = gLog2F32_Table[index].invC;
    f64 logC = gLog2F32_Table[index].logC;
    f64 z = (f64)u32f32(iz).f;

    /* log(x) = log1p(z/c-1) + log(c) + k*Ln2 */
    f64 r = z * invC - 1.0;
    f64 y0 = logC + (f64)k;

    f64 r2 = r * r;
    f64 result = gLog2F32_Poly[1] * r + gLog2F32_Poly[2];
    result = gLog2F32_Poly[0] * r2 + result;
    f64 p = gLog2F32_Poly[3] * r + y0;
    result = result * r2 + p;
    return result;
}

internal f32
log10_32(f32 x)
{
    s32 hx = (s32)u32f32(x).u;

    s32 k = 0;
    if (MATS_F32_UWORD_IS_ZERO(hx & MATS_F32_ABS_MASK)) {
        return -g2pow25F32 / 0.0f; // NOTE(michiel): log10(+/-0) = -inf
    } else if (hx < 0) {
        return (x - x) / 0.0f;     // NOTE(michiel): log10(-X) = NaN
    } else if (!MATS_F32_UWORD_IS_FINITE(hx)) {
        return x + x;
    } else if (MATS_F32_UWORD_IS_SUBNORMAL(hx)) {
        k -= 25;
        x *= g2pow25F32; // NOTE(michiel): x * 2^25
        hx = (s32)u32f32(x).u;
    }

    k += (hx >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;
    s32 i = (u32)k >> 31;
    hx = (hx & MATS_F32_MANT_MASK) | ((MATS_F32_EXP_BIAS - i) << MATS_F32_EXP_SHIFT);
    f32 y = (f32)(k + i);
    x = u32f32((u32)hx).f;
    f32 z = y * gLog10_2_lo + gInvLn10F32 * log32(x);
    return z + y * gLog10_2_hi;
}

internal f32
expm1_32(f32 x)
{
    // NOTE(michiel): result = exp(x) - 1
    u32 hx = u32f32(x).u;
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
    r1 = (r1 - 3.3333212137e-2) * hxs;
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
                u32 i = u32f32(y).u;
                y = u32f32(i + (k << MATS_F32_EXP_SHIFT)).f;
                result = y - 1.0f;
            }
            else if (k < MATS_F32_EXP_SHIFT)
            {
                f32 y = u32f32((u32)0x3F800000 - (0x01000000 >> k)).f;
                y = y - (e - x);
                u32 i = u32f32(y).u;
                result = u32f32(i + (k << MATS_F32_EXP_SHIFT)).f;
            }
            else
            {
                f32 y = u32f32((u32)(MATS_F32_EXP_BIAS - k) << MATS_F32_EXP_SHIFT).f;
                y = x - (e + y);
                y += 1.0f;
                u32 i = u32f32(y).u;
                result = u32f32(i + (k << MATS_F32_EXP_SHIFT)).f;
            }
        }
    }
    return result;
}

internal f32
log1p32(f32 x)
{
    // NOTE(michiel): result = log(1 + x)
    s32 hx = (s32)u32f32(x).u;
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
            hu = (s32)u32f32(u).u;
            k  = (hu >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;
            // NOTE(michiel): Correction term
            c = (k > 0) ? 1.0f - (u - x) : x - (u - 1.0f);
            c /= u;
        }
        else
        {
            f32 u = x;
            hu = (s32)u32f32(u).u;
            k  = (hu >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;
        }
        hu &= MATS_F32_MANT_MASK;
        f32 u;
        if (hu < 0x003504F7)
        {
            u = u32f32((u32)hu | 0x3F800000).f; // NOTE(michiel): Normalize
        }
        else
        {
            ++k;
            u = u32f32((u32)hu | 0x3F000000).f;
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
    u32 ix = u32f32(x).u;

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
        U32F32 u = u32f32(1.0f + x);
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

