//
// NOTE(michiel): Sqrt
//

internal f32_4x
sqrt32_4x(f32_4x x)
{
    f32_4x result;
    result.m = _mm_sqrt_ps(x.m);
    return result;
}

//
// NOTE(michiel): Hypot
//

internal f32_4x
hypot32_fast_4x(f32_4x x, f32_4x y)
{
    return sqrt32_4x(x * x + y * y);
}

internal f32_4x
hypot32_4x(f32_4x x, f32_4x y)
{
    f32_4x xu = s32_4x_and(x, S32_4x(MATS_F32_ABS_MASK));
    f32_4x yu = s32_4x_and(y, S32_4x(MATS_F32_ABS_MASK));

    f32_4x xLessY = s32_4x_less(xu, yu);

    x = select4x(xu, xLessY, yu);
    y = select4x(yu, xLessY, xu);

    xu = x;
    yu = y;

    f32_4x z = F32_4x(1.0f);

    f32_4x yZero = s32_4x_equal(yu, zero_s32_4x());
    f32_4x zeroRes = x;

    f32_4x zMax = F32_4x(0x1p90f);
    f32_4x zMin = F32_4x(0x1p-90f);

    f32_4x xBig = x * zMax;
    f32_4x xSmall = x * zMin;
    f32_4x yBig = y * zMax;
    f32_4x ySmall = y * zMin;

    f32_4x xOkay = s32_4x_less(xu, S32_4x((MATS_F32_EXP_BIAS + 60) << MATS_F32_EXP_SHIFT));
    f32_4x yLess = s32_4x_less(yu, S32_4x((MATS_F32_EXP_BIAS - 60) << MATS_F32_EXP_SHIFT));

    z = select4x(z, yLess, zMin);
    z = select4x(zMax, xOkay, z);

    x = select4x(x, yLess, xBig);
    x = select4x(xSmall, xOkay, x);

    y = select4x(y, yLess, yBig);
    y = select4x(ySmall, xOkay, y);

    f64_2x xLo = F64_2x(x);
    f64_2x xHi = F64_2x(x, true);

    f64_2x yLo = F64_2x(y);
    f64_2x yHi = F64_2x(y, true);

    f64_2x x2Lo = xLo * xLo;
    f64_2x x2Hi = xHi * xHi;
    f64_2x y2Lo = yLo * yLo;
    f64_2x y2Hi = yHi * yHi;

    f64_2x preLo = x2Lo + y2Lo;
    f64_2x preHi = x2Hi + y2Hi;

    f32_4x sqrtIn = F32_4x(preLo, preHi);
    f32_4x result = z * sqrt32_4x(sqrtIn);
    result = select4x(result, yZero, zeroRes);

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
    result = select4x(result, underflowMask, underflow);
    result = select4x(result, overflowMask, overflow);
    result = select4x(result, infOutMask, infOut);
    result = and_not(result, zeroMask);

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

    f64_2x kdLo = xdLo + F64_2x(gExp2F32_ShiftScaled);
    f64_2x kdHi = xdHi + F64_2x(gExp2F32_ShiftScaled);
    u64 ki0 = kdLo.u[0];
    u64 ki1 = kdLo.u[1];
    u64 ki2 = kdHi.u[0];
    u64 ki3 = kdHi.u[1];
    kdLo = kdLo - F64_2x(gExp2F32_ShiftScaled);
    kdHi = kdHi - F64_2x(gExp2F32_ShiftScaled);;

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

    result = select4x(result, underflowMask, underflow);
    result = select4x(result, overflowMask, overflow);
    result = select4x(result, infOutMask, infOut);
    result = select4x(result, zeroMask, zero_f32_4x());

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
    f32_4x index = s32_4x_and(s32_4x_srl(temp, MATS_F32_EXP_SHIFT - LOGF_TABLE_BITS),
                              S32_4x((1 << LOGF_TABLE_BITS) - 1));
    f32_4x k = s32_4x_sra(temp, MATS_F32_EXP_SHIFT);
    f32_4x iz = s32_4x_sub(x, s32_4x_and(temp, S32_4x(0x1FF << MATS_F32_EXP_SHIFT)));

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
    f32_4x fullExp   = S32_4x(MATS_F32_EXP_MASK);
    f32_4x signMask  = S32_4x(MATS_F32_SIGN_MASK);
    f32_4x maxExp    = S32_4x(0xFF000000);
    f32_4x fullBits  = S32_4x(0xFFFFFFFF);

    f32_4x maxX      = S32_4x(MATS_F32_EXP_MASK - 0x00800000);

    f32_4x xuTimes2  = s32_4x_add(x, x);
    f32_4x xuMinLowExp = s32_4x_sub(x, lowExpBit);

    f32_4x divZeroResult = zero_f32_4x();
    f32_4x infResult     = x;
    f32_4x invalidResult = F32_4x(F32_INF); // TODO(michiel): Or should we return 0 here

    f32_4x noErrorMask = s32_4x_less(xuMinLowExp, maxX);
    f32_4x divZeroMask = s32_4x_and_not(s32_4x_equal(xuTimes2, zero_f32_4x()), noErrorMask);
    f32_4x infMask     = s32_4x_and_not(s32_4x_equal(x, fullExp), noErrorMask);
    f32_4x invalidMask = s32_4x_and_not(s32_4x_or(s32_4x_and(x, signMask),
                                                  s32_4x_xor(s32_4x_less(xuTimes2, maxExp), fullBits)),
                                        noErrorMask);

    f32_4x subNormalFix = x * F32_4x(0x1p23f);
    subNormalFix = s32_4x_sub(subNormalFix, S32_4x(MATS_F32_EXP_SHIFT << MATS_F32_EXP_SHIFT));

    x = select4x(subNormalFix, noErrorMask, x);

    f32_4x result = log32_fast_4x(x);
    result = select4x(result, invalidMask, invalidResult);
    result = select4x(result, infMask, infResult);
    result = select4x(result, divZeroMask, divZeroResult);

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
    f32_4x index = s32_4x_and(s32_4x_srl(temp, MATS_F32_EXP_SHIFT - LOG2F_TABLE_BITS),
                              S32_4x((1 << LOG2F_TABLE_BITS) - 1));
    f32_4x k = s32_4x_sra(temp, MATS_F32_EXP_SHIFT);
    f32_4x iz = s32_4x_sub(x, s32_4x_and(temp, S32_4x(0x1FF << MATS_F32_EXP_SHIFT)));

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
    f32_4x fullExp   = S32_4x(MATS_F32_EXP_MASK);
    f32_4x signMask  = S32_4x(MATS_F32_SIGN_MASK);
    f32_4x maxExp    = S32_4x(0xFF000000);
    f32_4x fullBits  = S32_4x(0xFFFFFFFF);

    f32_4x maxX      = S32_4x(MATS_F32_EXP_MASK - 0x00800000);

    f32_4x xuTimes2  = s32_4x_add(x, x);
    f32_4x xuMinLowExp = s32_4x_sub(x, lowExpBit);

    f32_4x divZeroResult = zero_f32_4x();
    f32_4x infResult     = x;
    f32_4x invalidResult = F32_4x(F32_INF); // TODO(michiel): Or should we return 0 here

    f32_4x noErrorMask = s32_4x_less(xuMinLowExp, maxX);
    f32_4x divZeroMask = s32_4x_and_not(s32_4x_equal(xuTimes2, zero_f32_4x()), noErrorMask);
    f32_4x infMask     = s32_4x_and_not(s32_4x_equal(x, fullExp), noErrorMask);
    f32_4x invalidMask = s32_4x_and_not(s32_4x_or(s32_4x_and(x, signMask),
                                                  s32_4x_xor(s32_4x_less(xuTimes2, maxExp), fullBits)),
                                        noErrorMask);

    f32_4x subNormalFix = x * F32_4x(0x1p23f);
    subNormalFix = s32_4x_sub(subNormalFix, S32_4x(MATS_F32_EXP_SHIFT << MATS_F32_EXP_SHIFT));

    x = select4x(subNormalFix, noErrorMask, x);

    f32_4x result = log2_32_fast_4x(x);
    result = select4x(result, invalidMask, invalidResult);
    result = select4x(result, infMask, infResult);
    result = select4x(result, divZeroMask, divZeroResult);

    return result;
}

internal f32_4x
log10_32_fast_4x(f32_4x x, f32_4x k = zero_f32_4x())
{
    f32_4x invLn10    = F32_4x(gInvLn10F32);
    f32_4x log10_2hi  = F32_4x(gLog10_2_hi);
    f32_4x log10_2lo  = F32_4x(gLog10_2_lo);

    k = s32_4x_add(k, s32_4x_sub(s32_4x_sra(x, MATS_F32_EXP_SHIFT), S32_4x(MATS_F32_EXP_BIAS)));
    f32_4x i = s32_4x_srl(k, 31);
    x = s32_4x_and(x, S32_4x(MATS_F32_MANT_MASK)) | s32_4x_sll(s32_4x_sub(S32_4x(MATS_F32_EXP_BIAS), i), MATS_F32_EXP_SHIFT);
    f32_4x y = f32_4x_from_s32(s32_4x_add(k, i));
    f32_4x z = y * log10_2lo + invLn10 * log32_fast_4x(x);

    f32_4x result = z + y * log10_2hi;
    return result;
}

internal f32_4x
log10_32_4x(f32_4x x)
{
    f32_4x zero4x     = zero_f32_4x();
    f32_4x fullBits   = S32_4x(0xFFFFFFFF);
    f32_4x fullExp    = S32_4x(MATS_F32_EXP_MASK);
    f32_4x lowExpBit  = S32_4x(0x00800000);
    f32_4x pow2_25    = F32_4x(g2pow25F32);

    f32_4x zeroResult = F32_4x(-F32_INF);
    f32_4x negResult  = F32_4x(F32_NAN) ^ (x & S32_4x(MATS_F32_SIGN_MASK));
    f32_4x infResult  = x + x;

    f32_4x zeroMask   = s32_4x_equal(s32_4x_and(x, S32_4x(MATS_F32_ABS_MASK)), zero4x);
    f32_4x errorMask  = zeroMask;
    f32_4x negMask    = and_not(s32_4x_less(x, zero4x), errorMask);
    errorMask         = s32_4x_or(errorMask, negMask);
    f32_4x infMask    = and_not(s32_4x_xor(s32_4x_less(x, fullExp), fullBits), errorMask);
    errorMask         = s32_4x_or(errorMask, infMask);

    f32_4x errorRes   = (zeroResult & zeroMask) | (negResult & negMask) | (infResult & infMask);

    f32_4x k          = zero4x;
    f32_4x subNrmMask = s32_4x_less(x, lowExpBit);

    f32_4x kMod       = s32_4x_sub(k, S32_4x(25));
    f32_4x xMod       = x * pow2_25;

    k = select4x(k, subNrmMask, kMod);
    x = select4x(x, subNrmMask, xMod);

    f32_4x result = log10_32_fast_4x(x, k);
    result = select4x(result, errorMask, errorRes);
    return result;
}

//
// NOTE(mich2iel): Pow
//

internal void
pow32_log2_4x(f32_4x x, f64_2x *dstLo, f64_2x *dstHi)
{
    f64_2x poly0 = F64_2x(gPowF32_Log2Poly[0]);
    f64_2x poly1 = F64_2x(gPowF32_Log2Poly[1]);
    f64_2x poly2 = F64_2x(gPowF32_Log2Poly[2]);
    f64_2x poly3 = F64_2x(gPowF32_Log2Poly[3]);
    f64_2x poly4 = F64_2x(gPowF32_Log2Poly[4]);
    f32_4x temp = s32_4x_sub(x, S32_4x(0x3F330000));
    f32_4x index = s32_4x_and(s32_4x_srl(temp, MATS_F32_EXP_SHIFT - POWF_LOG2_TABLE_BITS),
                              S32_4x((1 << POWF_LOG2_TABLE_BITS) - 1));
    f32_4x iz = s32_4x_sub(x, s32_4x_and(temp, S32_4x(0xFF800000)));
    f32_4x k = s32_4x_sra(temp, MATS_F32_EXP_SHIFT);
    f64_2x k64Lo; k64Lo.md = _mm_cvtepi32_pd(k.mi);
    f64_2x k64Hi; k64Hi.md = _mm_cvtepi32_pd(_mm_castps_si128(_mm_shuffle_ps(k.m, k.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3))));

    // TODO(michiel): Hmm...
    f64_2x invCLo = F64_2x(gPowF32_Log2Table[index.u[0]].invC, gPowF32_Log2Table[index.u[1]].invC);
    f64_2x invCHi = F64_2x(gPowF32_Log2Table[index.u[2]].invC, gPowF32_Log2Table[index.u[3]].invC);
    f64_2x logCLo = F64_2x(gPowF32_Log2Table[index.u[0]].logC, gPowF32_Log2Table[index.u[1]].logC);
    f64_2x logCHi = F64_2x(gPowF32_Log2Table[index.u[2]].logC, gPowF32_Log2Table[index.u[3]].logC);

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

internal f32_4x
pow32_exp2_4x(f64_2x xLo, f64_2x xHi, f32_4x signBias)
{
    f64_2x poly0 = F64_2x(gExp2F32_Poly[0]);
    f64_2x poly1 = F64_2x(gExp2F32_Poly[1]);
    f64_2x poly2 = F64_2x(gExp2F32_Poly[2]);

    f64_2x kdLo = xLo + F64_2x(gExp2F32_ShiftScaled);
    f64_2x kdHi = xHi + F64_2x(gExp2F32_ShiftScaled);
    u64 ki0 = kdLo.u[0];
    u64 ki1 = kdLo.u[1];
    u64 ki2 = kdHi.u[0];
    u64 ki3 = kdHi.u[1];
    kdLo = kdLo - F64_2x(gExp2F32_ShiftScaled);
    kdHi = kdHi - F64_2x(gExp2F32_ShiftScaled);

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

    f32_4x result = F32_4x(resultLo, resultHi);
    return result;
}

internal f32_4x
pow32_4x(f32_4x x, f32_4x y)
{
    // NOTE(michiel): We don't really care about the results for bad inputs!
    f32_4x signBias = zero_f32_4x();
    f32_4x expManMask = S32_4x(MATS_F32_ABS_MASK);

    f32_4x xu = x & expManMask;
    f32_4x yu = y & expManMask;

    f32_4x xIsZero = s32_4x_equal(xu, zero_s32_4x());
    f32_4x yIsZero = s32_4x_equal(yu, zero_s32_4x());

    f32_4x one4x = F32_4x(1.0f);
    f32_4x oneOverX = reciprocal(x);

    f32_4x xSign = s32_4x_less(x, zero_s32_4x());
    f32_4x ySign = s32_4x_less(y, zero_s32_4x());

    f32_4x xZeroResult = select4x(x, ySign, oneOverX);

    //s32 e = (yu & MATS_F32_EXP_MASK >> MATS_F32_EXP_SHIFT;
    //s32 m = (1 << (150 - e));
    f32_4x e = s32_4x_srl(yu, MATS_F32_EXP_SHIFT);
    f32_4x shiftCount = s32_4x_sub(S32_4x(MATS_F32_EXP_BIAS + MATS_F32_EXP_SHIFT), e);
    // TODO(michiel): Ieuw
    f32_4x yMask;
    yMask.u[0] = 1 << shiftCount.u[0];
    yMask.u[1] = 1 << shiftCount.u[1];
    yMask.u[2] = 1 << shiftCount.u[2];
    yMask.u[3] = 1 << shiftCount.u[3];
    yMask = s32_4x_equal(s32_4x_and(yMask, yu), zero_s32_4x());
    yMask = s32_4x_and_not(xSign, yMask);
    f32_4x modSign = S32_4x(1 << (EXP2F_TABLE_BITS + 11));
    signBias = s32_4x_and(modSign, yMask);

    f64_2x logXLo, logXHi;
    pow32_log2_4x(xu, &logXLo, &logXHi);
    f64_2x yLogXLo = F64_2x(y) * logXLo;
    f64_2x yLogXHi = F64_2x(y, true) * logXHi;

    f32_4x result = pow32_exp2_4x(yLogXLo, yLogXHi, signBias);

    result = select4x(result, xIsZero, xZeroResult);
    result = select4x(result, yIsZero, one4x);

    return result;
}

internal f32_4x
pow10_32_4x(f32_4x x)
{
    return pow32_4x(F32_4x(10.0f), x);
}
