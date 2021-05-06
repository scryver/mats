//
// NOTE(michiel): Sqrt/Exp/Exp2/Log/Log2
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
    u32 hx = ix & 0x7FFFFFFF;

    /* take care of Inf and NaN */
	if(!FLT_UWORD_IS_FINITE(hx)) {
	    return x * x + x;		// sqrt(NaN)=NaN, sqrt(+inf)=+inf, sqrt(-inf)=sNaN
    }
    /* take care of zero and -ves */
	if (FLT_UWORD_IS_ZERO(hx)) {
        return x; // sqrt(+-0) = +-0
    }
	if (ix < 0) {
        return (x - x) / (x - x);		/* sqrt(-ve) = sNaN */
    }

    /* normalize x */
	s32 m = (ix >> 23);
	if (FLT_UWORD_IS_SUBNORMAL(hx))
    {   /* subnormal x */
        s32 i;
	    for (i = 0; (ix & 0x00800000) == 0; ++i) {
            ix <<= 1;
        }
	    m -= i - 1;
	}
	m -= 127;	/* unbias exponent */
	ix = (ix & 0x007FFFFF) | 0x00800000;
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
	ix += (m << 23);

    result = u32f32((u32)ix).f;
#endif
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
    z = gExp2F32_PolyScaled[0] * r + gExp2F32_PolyScaled[1];
    f64 r2 = r * r;
    f64 y = gExp2F32_PolyScaled[2] * r + 1.0;
    y = z * r2 + y;
    y = y * s;
    return y;
}

internal f32
exp2_32(f32 x)
{
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
#if 0 // TOINT_INTRINSICS
    // TODO(michiel): In 4x we can use round/cvt
    f64 kd = roundtoint(z);
    u64 ki = converttoint(z);
#else
    f64 kd = xd + gExp2F32_ShiftScaled; /* Rounding to double precision is required.  */
    u64 ki = u64f64(kd).u;
    kd -= gExp2F32_ShiftScaled;
#endif

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
log32(f32 x)
{
    u32 xu = u32f32(x).u;

    if ((xu - 0x00800000) >= (0x7F800000 - 0x00800000))
    {
        /* x < 0x1p-126 or inf or nan.  */
        if (xu * 2 == 0) {
            return mats_divzero32(1);
        }
        if (xu == 0x7F800000) {
            return x;
        }
        if ((xu & 0x80000000) || ((xu * 2) >= 0xFF000000)) {
            return mats_invalid32(x);
        }
        // NOTE(michiel): Normalize a subnormal
        xu = u32f32(x * 0x1p23f).u;
        xu -= 23 << 23;
    }

    /* x = 2^k z; where z is in range [OFF,2*OFF] and exact.
            The range is split into N subintervals.
            The ith subinterval contains z and c is near its center.  */
    u32 temp = xu - 0x3F330000;
    u32 index = (temp >> (23 - LOGF_TABLE_BITS)) % (1 << LOGF_TABLE_BITS);
    s32 k = (s32)temp >> 23;
    u32 iz = xu - (temp & 0x1FF << 23);
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

    if ((xu - 0x00800000) >= (0x7F800000 - 0x00800000))
    {
        /* x < 0x1p-126 or inf or nan.  */
        if (xu * 2 == 0) {
            return mats_divzero32(1);
        }
        if (xu == 0x7F800000) {
            return x;
        }
        if ((xu & 0x80000000) || ((xu * 2) >= 0xFF000000)) {
            return mats_invalid32(x);
        }
        // NOTE(michiel): Normalize a subnormal
        xu = u32f32(x * 0x1p23f).u;
        xu -= 23 << 23;
    }

    /* x = 2^k z; where z is in range [OFF,2*OFF] and exact.
            The range is split into N subintervals.
            The ith subinterval contains z and c is near its center.  */
    u32 temp = xu - 0x3F330000;
    u32 index = (temp >> (23 - LOG2F_TABLE_BITS)) % (1 << LOG2F_TABLE_BITS);
    u32 top = temp & 0xFF800000;
    u32 iz = xu - top;
    s32 k = (s32)temp >> 23;
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
