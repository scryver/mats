
internal f32_4x
sqrt32_4x(f32_4x x)
{
    f32_4x result;
    result.m = _mm_sqrt_ps(x.m);
    return result;
}

internal f32_4x
exp32_4x(f32_4x x)
{
    f32_4x topMask   = S32_4x(0x7FF00000);
    f32_4x f88top    = s32_4x_and(F32_4x(88.0f), topMask);

    f32_4x infOut    = x + x;
    f32_4x overflow  = square(F32_4x(0x1p97f));
    f32_4x underflow = square(F32_4x(0x1p-95f));

    f32_4x abstop        = s32_4x_and(x, topMask);
    f32_4x errorMask     = s32_4x_less(f88top, abstop);

    f32_4x zeroMask      = s32_4x_and(s32_4x_equal(x, F32_4x(-F32_INF)), errorMask);
    f32_4x infOutMask    = s32_4x_and(s32_4x_less(s32_4x_and(F32_4x(F32_INF), topMask), abstop), errorMask);
    f32_4x overflowMask  = s32_4x_and(x > F32_4x(0x1.62E42Ep6f), errorMask);
    f32_4x underflowMask = s32_4x_and(x < F32_4x(-0x1.9FE368p6f), errorMask);

    f64_2x xdLo = F64_2x(x);
    f64_2x xdHi = F64_2x(x, true);

    f64_2x zLo = F64_2x(gExp2F32_InvLn2Scaled) * xdLo;
    f64_2x zHi = F64_2x(gExp2F32_InvLn2Scaled) * xdHi;

    f64_2x kdLo = round(zLo);
    f64_2x kdHi = round(zHi);
    // TODO(michiel): Should this truncate or round?
    u64 ki0 = (u64)s64_from_f64_2x(zLo);
    u64 ki1 = (u64)s64_from_f64_2x(zLo, true);
    u64 ki2 = (u64)s64_from_f64_2x(zHi);
    u64 ki3 = (u64)s64_from_f64_2x(zHi, true);

    u64 t0 = gExp2F32_Table[ki0 % (1 << EXP2F_TABLE_BITS)];
    u64 t1 = gExp2F32_Table[ki1 % (1 << EXP2F_TABLE_BITS)];
    t0 += ki0 << (52 - EXP2F_TABLE_BITS);
    t1 += ki1 << (52 - EXP2F_TABLE_BITS);

    u64 t2 = gExp2F32_Table[ki2 % (1 << EXP2F_TABLE_BITS)];
    u64 t3 = gExp2F32_Table[ki3 % (1 << EXP2F_TABLE_BITS)];
    t2 += ki2 << (52 - EXP2F_TABLE_BITS);
    t3 += ki3 << (52 - EXP2F_TABLE_BITS);

    f64_2x rLo = zLo - kdLo;
    f64_2x rHi = zHi - kdHi;

    f64_2x sLo = S64_2x(t0, t1);
    f64_2x sHi = S64_2x(t2, t3);

#define EXP2F_N  ((f64)(1 << EXP2F_TABLE_BITS))
    f64_2x poly0 = F64_2x(0x1.c6af84b912394p-5 / EXP2F_N / EXP2F_N / EXP2F_N);
    f64_2x poly1 = F64_2x(0x1.ebfce50fac4f3p-3 / EXP2F_N / EXP2F_N);
    f64_2x poly2 = F64_2x(0x1.62e42ff0c52d6p-1 / EXP2F_N);
#undef EXP2F_N
    zLo = poly0 * rLo + poly1;
    zHi = poly0 * rHi + poly1;

    f64_2x r2Lo = square(rLo);
    f64_2x r2Hi = square(rHi);

    f64_2x yLo = poly2 * rLo + F64_2x(1.0);
    f64_2x yHi = poly2 * rHi + F64_2x(1.0);
    yLo += zLo * r2Lo;
    yHi += zHi * r2Hi;
    yLo *= sLo;
    yHi *= sHi;

    f32_4x y2Lo = F32_4x(yLo);
    f32_4x y2Hi = F32_4x(yHi);
    f32_4x result;
    result.m = _mm_shuffle_ps(y2Lo.m, y2Hi.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));

#if MATS_USE_SSE4
    result.m = _mm_blendv_ps(result.m, underflow.m, underflowMask.m);
    result.m = _mm_blendv_ps(result.m, overflow.m, overflowMask.m);
    result.m = _mm_blendv_ps(result.m, infOut.m, infOutMask.m);
    result.m = _mm_blendv_ps(result.m, _mm_setzero_ps(), zeroMask.m);
#else
    result = select(result, underflowMask, underflow);
    result = select(result, overflowMask, overflow);
    result = select(result, infOutMask, infOut);
    result = select(result, zeroMask, zero_f32_4x());
#endif

    return result;
}
