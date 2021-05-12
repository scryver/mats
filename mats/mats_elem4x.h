
internal f32_4x
sqrt32_4x(f32_4x x)
{
    f32_4x result;
    result.m = _mm_sqrt_ps(x.m);
    return result;
}

//
// NOTE(michiel): Exp
//

internal f32_4x
exp32_fast_4x(f32_4x x)
{
    // NOTE(michiel): |x| < 88.0f
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

    f32_4x result = exp32_fast_4x(x);
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

//
// NOTE(michiel): Exp2
//

internal f32_4x
exp2_32_fast_4x(f32_4x x)
{
    // NOTE(michiel): |x| < 128.0f

    f64_2x xdLo = F64_2x(x);
    f64_2x xdHi = F64_2x(x, true);

    f64_2x zLo = xdLo;
    f64_2x zHi = xdHi;

#if 0
    f64_2x kdLo = round(zLo);
    f64_2x kdHi = round(zHi);
    // TODO(michiel): Should this truncate or round?
    u64 ki0 = (u64)s64_from_f64_2x(zLo);
    u64 ki1 = (u64)s64_from_f64_2x(zLo, true);
    u64 ki2 = (u64)s64_from_f64_2x(zHi);
    u64 ki3 = (u64)s64_from_f64_2x(zHi, true);
#else
    f64_2x kdLo = xdLo + F64_2x(gExp2F32_ShiftScaled);
    f64_2x kdHi = xdHi + F64_2x(gExp2F32_ShiftScaled);
    u64 ki0 = kdLo.u[0];
    u64 ki1 = kdLo.u[1];
    u64 ki2 = kdHi.u[0];
    u64 ki3 = kdHi.u[1];
    kdLo = kdLo - F64_2x(gExp2F32_ShiftScaled);
    kdHi = kdHi - F64_2x(gExp2F32_ShiftScaled);;
#endif

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

    f64_2x poly0 = F64_2x(0x1.c6af84b912394p-5);
    f64_2x poly1 = F64_2x(0x1.ebfce50fac4f3p-3);
    f64_2x poly2 = F64_2x(0x1.62e42ff0c52d6p-1);
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

    return result;
}

internal f32_4x
exp2_32_4x(f32_4x x)
{
    f32_4x topMask   = S32_4x(0x7FF00000);
    f32_4x f128top   = s32_4x_and(F32_4x(128.0f), topMask);

    f32_4x infOut    = x + x;
    f32_4x overflow  = square(F32_4x(0x1p97f));
    f32_4x underflow = square(F32_4x(0x1p-95f));

    f32_4x abstop        = s32_4x_and(x, topMask);
    f32_4x errorMask     = s32_4x_less(f128top, abstop);

    f32_4x zeroMask      = s32_4x_and(s32_4x_equal(x, F32_4x(-F32_INF)), errorMask);
    f32_4x infOutMask    = s32_4x_and(s32_4x_less(s32_4x_and(F32_4x(F32_INF), topMask), abstop), errorMask);
    f32_4x overflowMask  = s32_4x_and(x > zero_f32_4x(), errorMask);
    f32_4x underflowMask = s32_4x_and(x < F32_4x(-150.0f), errorMask);

    f32_4x result = exp2_32_fast_4x(x);

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

//
// NOTE(michiel): Log
//

internal f32_4x
log32_fast_4x(f32_4x x)
{
    // NOTE(michiel): x >= 0x1p-126, and not inf or nan
    f64_2x one2x = F64_2x(1.0);
    f64_2x ln2_2x = F64_2x(gLogF32_Ln2);
    f64_2x poly0 = F64_2x(gLogF32_Poly[0]);
    f64_2x poly1 = F64_2x(gLogF32_Poly[1]);
    f64_2x poly2 = F64_2x(gLogF32_Poly[2]);

    f32_4x temp = s32_4x_sub(x, S32_4x(0x3F330000));
    f32_4x index = s32_4x_and(s32_4x_srl(temp, 23 - LOGF_TABLE_BITS),
                              S32_4x((1 << LOGF_TABLE_BITS) - 1));
    f32_4x k = s32_4x_sra(temp, 23);
    f32_4x iz = s32_4x_sub(x, s32_4x_and(temp, S32_4x(0x1FF << 23)));

    // TODO(michiel): Hmm...
    f64_2x invCLo = F64_2x(gLogF32_Table[index.u[0]].invC, gLogF32_Table[index.u[1]].invC);
    f64_2x invCHi = F64_2x(gLogF32_Table[index.u[2]].invC, gLogF32_Table[index.u[3]].invC);
    f64_2x logCLo = F64_2x(gLogF32_Table[index.u[0]].logC, gLogF32_Table[index.u[1]].logC);
    f64_2x logCHi = F64_2x(gLogF32_Table[index.u[2]].logC, gLogF32_Table[index.u[3]].logC);

    f64_2x zLo = F64_2x(iz);
    f64_2x zHi = F64_2x(iz, true);

    f64_2x rLo = zLo * invCLo - one2x;
    f64_2x rHi = zHi * invCHi - one2x;
    f64_2x k64Lo; k64Lo.md = _mm_cvtepi32_pd(k.mi);
    f64_2x k64Hi; k64Hi.md = _mm_cvtepi32_pd(_mm_castps_si128(_mm_shuffle_ps(k.m, k.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3))));
    f64_2x y0Lo = logCLo + k64Lo * ln2_2x;
    f64_2x y0Hi = logCHi + k64Hi * ln2_2x;

    f64_2x r2Lo = rLo * rLo;
    f64_2x r2Hi = rHi * rHi;
    f64_2x resultLo = poly1 * rLo + poly2;
    f64_2x resultHi = poly1 * rHi + poly2;
    resultLo = poly0 * r2Lo + resultLo;
    resultHi = poly0 * r2Hi + resultHi;
    resultLo = resultLo * r2Lo + (y0Lo + rLo);
    resultHi = resultHi * r2Hi + (y0Hi + rHi);

    f32_4x resLo = F32_4x(resultLo);
    f32_4x resHi = F32_4x(resultHi);
    f32_4x result;
    result.m = _mm_shuffle_ps(resLo.m, resHi.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));

    return result;
}

internal f32_4x
log32_4x(f32_4x x)
{
    f32_4x lowExpBit = S32_4x(0x00800000);
    f32_4x fullExp   = S32_4x(0x7F800000);
    f32_4x signMask  = S32_4x(0x80000000);
    f32_4x maxExp    = S32_4x(0xFF000000);
    f32_4x fullBits  = S32_4x(0xFFFFFFFF);

    f32_4x maxX      = S32_4x(0x7F800000 - 0x00800000);

    f32_4x xuTimes2  = s32_4x_add(x, x);
    f32_4x xuMinLowExp = s32_4x_sub(x, lowExpBit);

    f32_4x divZeroResult = zero_f32_4x();
    f32_4x infResult     = x;
    f32_4x invalidResult = F32_4x(F32_INF); // TODO(michiel): Or should we return 0 here

    f32_4x errorMask   = s32_4x_xor(s32_4x_less(xuMinLowExp, maxX), fullBits);
    f32_4x divZeroMask = s32_4x_and(s32_4x_equal(xuTimes2, zero_f32_4x()), errorMask);
    f32_4x infMask     = s32_4x_and(s32_4x_equal(x, fullExp), errorMask);
    f32_4x invalidMask = s32_4x_and(s32_4x_or(s32_4x_and(x, signMask),
                                              s32_4x_xor(s32_4x_less(xuTimes2, maxExp), fullBits)),
                                    errorMask);

    f32_4x subNormalFix = x * F32_4x(0x1p23f);
    subNormalFix = s32_4x_sub(subNormalFix, S32_4x(23 << 23));

    x = select(x, errorMask, subNormalFix);

    f32_4x result = log32_fast_4x(x);
    result = select(result, invalidMask, invalidResult);
    result = select(result, infMask, infResult);
    result = select(result, divZeroMask, divZeroResult);

    return result;
}

//
// NOTE(michiel): Log2
//

internal f32_4x
log2_32_fast_4x(f32_4x x)
{
    // NOTE(michiel): x >= 0x1p-126, and not inf or nan
    f64_2x one2x = F64_2x(1.0);
    f64_2x poly0 = F64_2x(gLog2F32_Poly[0]);
    f64_2x poly1 = F64_2x(gLog2F32_Poly[1]);
    f64_2x poly2 = F64_2x(gLog2F32_Poly[2]);
    f64_2x poly3 = F64_2x(gLog2F32_Poly[3]);

    f32_4x temp = s32_4x_sub(x, S32_4x(0x3F330000));
    f32_4x index = s32_4x_and(s32_4x_srl(temp, 23 - LOG2F_TABLE_BITS),
                              S32_4x((1 << LOG2F_TABLE_BITS) - 1));
    f32_4x k = s32_4x_sra(temp, 23);
    f32_4x iz = s32_4x_sub(x, s32_4x_and(temp, S32_4x(0x1FF << 23)));

    // TODO(michiel): Hmm...
    f64_2x invCLo = F64_2x(gLog2F32_Table[index.u[0]].invC, gLog2F32_Table[index.u[1]].invC);
    f64_2x invCHi = F64_2x(gLog2F32_Table[index.u[2]].invC, gLog2F32_Table[index.u[3]].invC);
    f64_2x logCLo = F64_2x(gLog2F32_Table[index.u[0]].logC, gLog2F32_Table[index.u[1]].logC);
    f64_2x logCHi = F64_2x(gLog2F32_Table[index.u[2]].logC, gLog2F32_Table[index.u[3]].logC);

    f64_2x zLo = F64_2x(iz);
    f64_2x zHi = F64_2x(iz, true);

    f64_2x rLo = zLo * invCLo - one2x;
    f64_2x rHi = zHi * invCHi - one2x;
    f64_2x k64Lo; k64Lo.md = _mm_cvtepi32_pd(k.mi);
    f64_2x k64Hi; k64Hi.md = _mm_cvtepi32_pd(_mm_castps_si128(_mm_shuffle_ps(k.m, k.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3))));
    f64_2x y0Lo = logCLo + k64Lo;
    f64_2x y0Hi = logCHi + k64Hi;

    f64_2x r2Lo = rLo * rLo;
    f64_2x r2Hi = rHi * rHi;
    f64_2x resultLo = poly1 * rLo + poly2;
    f64_2x resultHi = poly1 * rHi + poly2;
    resultLo = poly0 * r2Lo + resultLo;
    resultHi = poly0 * r2Hi + resultHi;
    f64_2x pLo = poly3 * rLo + y0Lo;
    f64_2x pHi = poly3 * rHi + y0Hi;
    resultLo = resultLo * r2Lo + pLo;
    resultHi = resultHi * r2Hi + pHi;

    f32_4x resLo = F32_4x(resultLo);
    f32_4x resHi = F32_4x(resultHi);
    f32_4x result;
    result.m = _mm_shuffle_ps(resLo.m, resHi.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));

    return result;
}

internal f32_4x
log2_32_4x(f32_4x x)
{
    f32_4x lowExpBit = S32_4x(0x00800000);
    f32_4x fullExp   = S32_4x(0x7F800000);
    f32_4x signMask  = S32_4x(0x80000000);
    f32_4x maxExp    = S32_4x(0xFF000000);
    f32_4x fullBits  = S32_4x(0xFFFFFFFF);

    f32_4x maxX      = S32_4x(0x7F800000 - 0x00800000);

    f32_4x xuTimes2  = s32_4x_add(x, x);
    f32_4x xuMinLowExp = s32_4x_sub(x, lowExpBit);

    f32_4x divZeroResult = zero_f32_4x();
    f32_4x infResult     = x;
    f32_4x invalidResult = F32_4x(F32_INF); // TODO(michiel): Or should we return 0 here

    f32_4x errorMask   = s32_4x_xor(s32_4x_less(xuMinLowExp, maxX), fullBits);
    f32_4x divZeroMask = s32_4x_and(s32_4x_equal(xuTimes2, zero_f32_4x()), errorMask);
    f32_4x infMask     = s32_4x_and(s32_4x_equal(x, fullExp), errorMask);
    f32_4x invalidMask = s32_4x_and(s32_4x_or(s32_4x_and(x, signMask),
                                              s32_4x_xor(s32_4x_less(xuTimes2, maxExp), fullBits)),
                                    errorMask);

    f32_4x subNormalFix = x * F32_4x(0x1p23f);
    subNormalFix = s32_4x_sub(subNormalFix, S32_4x(23 << 23));

    x = select(x, errorMask, subNormalFix);

    f32_4x result = log2_32_fast_4x(x);
    result = select(result, invalidMask, invalidResult);
    result = select(result, infMask, infResult);
    result = select(result, divZeroMask, divZeroResult);

    return result;
}

#if 0
//
// NOTE(mich2iel): Pow
//

internal void
pow32_log2_4x(f32_4x x, f64_2x *dstLo, f64_2x *dstHi)
{
    f64_2x poly0 = F64_2x(gPowF32_Log2PolyScaled[0]);
    f64_2x poly1 = F64_2x(gPowF32_Log2PolyScaled[1]);
    f64_2x poly2 = F64_2x(gPowF32_Log2PolyScaled[2]);
    f64_2x poly3 = F64_2x(gPowF32_Log2PolyScaled[3]);
    f64_2x poly4 = F64_2x(gPowF32_Log2PolyScaled[4]);
    f32_4x temp = s32_4x_sub(x, S32_4x(0x3F30000));
    f32_4x index = s32_4x_and(s32_4x_srl(temp, 23 - POWF_LOG2_TABLE_BITS),
                              S32_4x((1 << POWF_LOG2_TABLE_BITS) - 1));
    f32_4x iz = s32_4x_sub(x, s32_4x_and(temp, S32_4x(0x1FF << 23)));
    f32_4x k = s32_4x_sra(temp, (23 - EXP2F_TABLE_BITS));
    f64_2x k64Lo; k64Lo.md = _mm_cvtepi32_pd(k.mi);
    f64_2x k64Hi; k64Hi.md = _mm_cvtepi32_pd(_mm_castps_si128(_mm_shuffle_ps(k.m, k.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3))));

    // TODO(michiel): Hmm...
    f64_2x invCLo = F64_2x(gPowF32_Log2TableScaled[index.u[0]].invC, gPowF32_Log2TableScaled[index.u[1]].invC);
    f64_2x invCHi = F64_2x(gPowF32_Log2TableScaled[index.u[2]].invC, gPowF32_Log2TableScaled[index.u[3]].invC);
    f64_2x logCLo = F64_2x(gPowF32_Log2TableScaled[index.u[0]].logC, gPowF32_Log2TableScaled[index.u[1]].logC);
    f64_2x logCHi = F64_2x(gPowF32_Log2TableScaled[index.u[2]].logC, gPowF32_Log2TableScaled[index.u[3]].logC);

    f64_2x zLo = F64_2x(iz);
    f64_2x zHi = F64_2x(iz, true);
    f64_2x rLo = zLo * invCLo - F64_2x(1.0);
    f64_2x rHi = zHi * invCHi - F64_2x(1.0);
    f64_2x y0Lo = logCLo + k64Lo;
    f64_2x y0Hi = logCHi + k64Hi;

    f64_2x r2Lo = rLo * rLo;
    f64_2x r2Hi = rHi * rHi;
    f64_2x resultLo = poly0 * rLo + poly1;
    f64_2x resultHi = poly0 * rHi + poly1;
    f64_2x pLo = poly2 * rLo + poly3;
    f64_2x pHi = poly2 * rHi + poly3;
    f64_2x r4Lo = r2Lo * r2Lo;
    f64_2x r4Hi = r2Hi * r2Hi;
    f64_2x qLo = poly4 * rLo + y0Lo;
    f64_2x qHi = poly4 * rHi + y0Hi;
    qLo = pLo * r2Lo + qLo;
    qHi = pHi * r2Hi + qHi;
    resultLo = resultLo * r4Lo + qLo;
    resultHi = resultHi * r4Hi + qHi;
    *dstLo = resultLo;
    *dstHi = resultHi;
}

internal void
pow32_exp2_4x(f64_2x xLo, f64_2x xHi, f32_4x signBias, f64_2x *dstLo, f64_2x *dstHi)
{
    f64_2x poly0 = F64_2x(gExp2F32_PolyScaled[0]);
    f64_2x poly1 = F64_2x(gExp2F32_PolyScaled[1]);
    f64_2x poly2 = F64_2x(gExp2F32_PolyScaled[2]);

    f64_2x kdLo = round(xLo);
    f64_2x kdHi = round(xHi);
    // TODO(michiel): Should this truncate or round?
    u64 ki0 = (u64)s64_from_f64_2x(xLo);
    u64 ki1 = (u64)s64_from_f64_2x(xLo, true);
    u64 ki2 = (u64)s64_from_f64_2x(xHi);
    u64 ki3 = (u64)s64_from_f64_2x(xHi, true);

    f64_2x rLo = xLo - kdLo;
    f64_2x rHi = xHi - kdHi;

    u64 t0 = gExp2F32_Table[ki0 % (1 << EXP2F_TABLE_BITS)];
    u64 t1 = gExp2F32_Table[ki1 % (1 << EXP2F_TABLE_BITS)];
    u64 t2 = gExp2F32_Table[ki2 % (1 << EXP2F_TABLE_BITS)];
    u64 t3 = gExp2F32_Table[ki3 % (1 << EXP2F_TABLE_BITS)];
    t0 += (ki0 + signBias.u[0]) << (52 - EXP2F_TABLE_BITS);
    t1 += (ki1 + signBias.u[1]) << (52 - EXP2F_TABLE_BITS);
    t2 += (ki2 + signBias.u[2]) << (52 - EXP2F_TABLE_BITS);
    t3 += (ki3 + signBias.u[3]) << (52 - EXP2F_TABLE_BITS);

    f64_2x sLo = S64_2x(t0, t1);
    f64_2x sHi = S64_2x(t2, t3);

    f64_2x zLo = poly0 * rLo + poly1;
    f64_2x zHi = poly0 * rHi + poly1;
    f64_2x r2Lo = rLo * rLo;
    f64_2x r2Hi = rHi * rHi;
    f64_2x resultLo = poly2 * rLo + F64_2x(1.0);
    f64_2x resultHi = poly2 * rHi + F64_2x(1.0);
    resultLo = zLo * r2Lo + resultLo;
    resultHi = zHi * r2Hi + resultHi;
    resultLo = resultLo * sLo;
    resultHi = resultHi * sHi;
    *dstLo = resultLo;
    *dstHi = resultHi;
}

internal f32_4x
pow32_checkint_4x(f32_4x x)
{
    f32_4x result;
    result.u[0] = pow32_checkint(x.u[0]);
    result.u[1] = pow32_checkint(x.u[1]);
    result.u[2] = pow32_checkint(x.u[2]);
    result.u[3] = pow32_checkint(x.u[3]);
    return result;
}

internal f32_4x
pow32_zeroinfnan_4x(f32_4x x)
{
    // NOTE(michiel): We bias the comparison to bring 0 to the most negative value
    // we only have a signed comparison in simd land. Also we xor the result with
    // U32_MAX to flip the bits.
    f32_4x result = s32_4x_xor(s32_4x_less(s32_4x_xor(s32_4x_sub(s32_4x_add(x, x), S32_4x(1)), S32_4x(0x80000000)),
                                           S32_4x((2u * 0x7F800000 - 1) ^ 0x80000000)),
                               S32_4x(0xFFFFFFFF));
    return result;
}

internal f32
pow32_temp(f32 x, f32 y)
{
    u32 signBias = 0;

    u32 xu = u32f32(x).u;
    u32 yu = u32f32(y).u;

    b32 yIsZeroInfNan = pow32_zeroinfnan(yu);
    b32 xIsZeroInfNan = pow32_zeroinfnan(xu);
    s32 yInt = pow32_checkint(yu);

    if (xIsZeroInfNan)
    {
        f32 x2 = x * x;
        if ((xu & 0x80000000) && (yInt == 1)) {
            x2 = -x2;
        }
        return yu & 0x80000000 ? 1.0f / x2 : x2;
    }
    else if (yIsZeroInfNan)
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
    else if ((yInt == 0) && (xu & 0x80000000))
    {
        return mats_invalid32(x);
    }
    else if (((xu - 0x00800000) >= (0x7f800000 - 0x00800000)))
    {
        /* x and y are non-zero finite.  */
        if (xu & 0x80000000)
        {
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
#endif
