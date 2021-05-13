//
// NOTE(michiel): Pow
//

internal f64
pow32_log2(u32 xu)
{
    u32 temp = xu - 0x3F330000;
    u32 index = (temp >> (23 - POWF_LOG2_TABLE_BITS)) % (1 << POWF_LOG2_TABLE_BITS);
    u32 top = temp & 0xFF800000;
    u32 iz = xu - top;
    s32 k = (s32)temp >> 23;
    f64 invC = gPowF32_Log2Table[index].invC;
    f64 logC = gPowF32_Log2Table[index].logC;
    f64 z = (f64)u32f32(iz).f;
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
    u64 ki = u64f64(kd).u;
    kd -= gExp2F32_ShiftScaled;

    f64 r = xd - kd;
    u64 t = gExp2F32_Table[ki % (1 << EXP2F_TABLE_BITS)];
    u64 ski = ki + signBias;
    t += ski << (52 - EXP2F_TABLE_BITS);
    f64 s = u64f64(t).f;
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
    s32 e = (xu >> 23) & 0xFF;
    if (e < 0x7F) {
        return 0;
    }
    if (e > 0x7F + 23) {
        return 2;
    }
    if (xu & ((1 << (0x7F + 23 - e)) - 1)) {
        return 0;
    }
    if (xu & (1 << (0x7F + 23 - e))) {
        return 1;
    }
    return 2;
}

internal s32
pow32_zeroinfnan(u32 ix)
{
    return 2 * ix - 1 >= 2u * 0x7F800000 - 1;
}

internal f32
pow32(f32 x, f32 y)
{
    u32 signBias = 0;

    u32 xu = u32f32(x).u;
    u32 yu = u32f32(y).u;

    if (((xu - 0x00800000) >= (0x7f800000 - 0x00800000)) || pow32_zeroinfnan(yu))
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
            if (((2 * xu) > (2u * 0x7F800000)) ||
                ((2 * yu) > (2u * 0x7F800000))) {
                return x + y;
            }
            if ((2 * xu) == 0x3F800000) {
                return 1.0f;
            }
            if (((2 * xu) < (2 * 0x3F800000)) == !(yu & 0x80000000)) {
                return 0.0f; // NOTE(michiel): |x| < 1 && y == inf or |x| > 1 && y == -inf
            }
            return y * y;
        }
        if (pow32_zeroinfnan(xu))
        {
            f32 x2 = x * x;
            if ((xu & 0x80000000) &&
                (pow32_checkint(yu) == 1))
            {
                x2 = -x2;
                signBias = 1;
            }
            return yu & 0x80000000 ? 1.0f / x2 : x2;
        }
        /* x and y are non-zero finite.  */
        if (xu & 0x80000000)
        {
            s32 yInt = pow32_checkint(yu);
            if (yInt == 0) {
                return mats_invalid32(x);
            }
            if (yInt == 1) {
                signBias = (1 << (EXP2F_TABLE_BITS + 11));
            }
            xu &= 0x7FFFFFFF;
        }
        if (xu < 0x00800000)
        {
            xu = u32f32(x * 0x1p23f).u;
            xu &= 0x7FFFFFFF;
            xu -= 23 << 23;
        }
    }

    f64 logX = pow32_log2(xu);
    f64 yLogX = (f64)y * logX; /* Note: cannot overflow, y is single prec.  */
    if ((u64f64(yLogX).u >> 47 & 0xFFFF) >= (u64f64(126.0 * POWF_SCALE).u >> 47))
    {
        if (yLogX > 0x1.fffffffd1d571p+6 * POWF_SCALE) {
            return mats_overflow32(signBias);
        }
        if (yLogX <= -150.0 * POWF_SCALE) {
            return mats_underflow32(signBias);
        }
    }
    return pow32_exp2(yLogX, signBias);
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
