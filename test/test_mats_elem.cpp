#include <errno.h>
#include <time.h>
#include <math.h>
#include <x86intrin.h>
#include <xmmintrin.h>

#define STR_FMT(x)   safe_truncate_to_s32(x.size), (char *)x.data
#include "../libberdip/src/common.h"
#include "../libberdip/src/multilane.h"

#define IEEE_754_2008_SNAN 1
#include "../mats/mats_common.h"
#include "../mats/mats_constants.h"

#define sqrt32   sqrt32_nonsse
#define hypot32  hypot32_nonsse
#define exp32    exp32_nonsse
#define exp2_32  exp2_32_nonsse
#define pow2_32  pow2_32_nonsse
#define log32    log32_nonsse
#define log2_32  log2_32_nonsse
#define log10_32 log10_32_nonsse
#define MATS_USE_SSE2 0
#define MATS_USE_SSE4 0
#include "../mats/mats_defines.h"
#include "../mats/mats_elem.h"
#undef sqrt32
#undef hypot32
#undef exp32
#undef exp2_32
#undef pow2_32
#undef log32
#undef log2_32
#undef log10_32

#define sqrt32   sqrt32_sse
#define hypot32  hypot32_sse
#define exp32    exp32_not_used
#define exp2_32  exp2_32_not_used
#define log32    log32_not_used
#define log2_32  log2_32_not_used
#define pow2_32  pow2_32_not_used
#define log10_32 log10_32_not_used
#undef  MATS_USE_SSE2
#undef  MATS_USE_SSE4
#define MATS_USE_SSE2 1
#define MATS_USE_SSE4 1
#include "../mats/mats_defines.h"
#include "../mats/mats_elem.h"
#undef sqrt32
#undef hypot32
#undef exp32
#undef exp2_32
#undef log32
#undef log2_32
#undef pow2_32
#undef log10_32

internal f32
hypot_ieee754(f32 x, f32 y)
{
	s32 ha = (s32)u32f32(x).u & 0x7FFFFFFF;
    s32 hb = (s32)u32f32(y).u & 0x7FFFFFFF;

    if (hb > ha) {
        s32 temp = ha;
        ha = hb;
        hb = temp;
    }

    f32 a = u32f32((u32)ha).f;
    f32 b = u32f32((u32)hb).f;

    if ((ha - hb) > 0x0F000000) {
        return a + b; /* x/y > 2**30 */
    }

    s32 k = 0;
	if (ha > 0x58800000)
    {   /* a > 2**50 */
        if (!FLT_UWORD_IS_FINITE(ha))
        {
            f32 w = a + b;
            if (FLT_UWORD_IS_INFINITE(ha)) {
                w = a;
            } else if (FLT_UWORD_IS_INFINITE(hb)) {
                w = b;
            }
            return w;
        }
        /* scale a and b by 2**-68 */
        ha -= 0x22000000;
        hb -= 0x22000000;
        k  += 68;
        a = u32f32((u32)ha).f;
        b = u32f32((u32)hb).f;
	}

    if (hb < 0x26800000L) {
        /* b < 2**-50 */
        if (FLT_UWORD_IS_ZERO(hb)) {
            return a;
        } else if (FLT_UWORD_IS_SUBNORMAL(hb)) {
            f32 t = u32f32((u32)0x7E800000).f; /* t = 2^126 */
            a *= t;
            b *= t;
            k -= 126;
        } else { /* scale a and b by 2^68 */
            ha += 0x22000000; /* a *= 2^68 */
            hb += 0x22000000; /* b *= 2^68 */
            k  -= 68;
            a   = u32f32((u32)ha).f;
            b   = u32f32((u32)hb).f;
        }
	}

    /* medium size a and b */
    f32 w = a - b;
    if (w > b)
    {
        f32 t1 = u32f32((u32)ha & 0xFFFFF000).f;
        f32 t2 = a - t1;
        w = sqrt32_sse(t1 * t1 - (b * (-b) - t2 * (a + t1)));
    }
    else
    {
        a = a + a;
        f32 y1 = u32f32((u32)hb & 0xFFFFF000).f;
        f32 y2 = b - y1;
        f32 t1 = u32f32((u32)(ha + 0x00800000) & 0xFFFFF000).f;
        f32 t2 = a - t1;
        w = sqrt32_sse(t1 * y1 - (w * (-w) - (t1 * y2 + t2 * b)));
    }

    if (k != 0)
    {
        f32 t1 = u32f32((u32)0x3F800000 + (k << 23)).f;
        w = w * t1;
    }
    return w;
}

internal f32
exp32_sse(f32 x)
{
    WideMath xw;
    xw.m = _mm_set_ss(x);
    WideMath topMask;
    topMask.mi = _mm_set1_epi32(0x7FF00000);
    WideMath f88top;
    f88top.mi = _mm_and_si128(_mm_castps_si128(_mm_set_ss(88.0f)),
                              topMask.mi);
    WideMath overflow;
    overflow.m = _mm_set_ss(0x1p97f);
    overflow.m = _mm_mul_ss(overflow.m, overflow.m);
    WideMath underflow;
    underflow.m = _mm_set_ss(0x1p-95f);
    underflow.m = _mm_mul_ss(underflow.m, underflow.m);

    WideMath errors;
    WideMath abstop;
    abstop.mi = _mm_and_si128(xw.mi, topMask.mi);
    WideMath errorMask;
    errorMask.mi = _mm_cmplt_epi32(f88top.mi, abstop.mi);

    WideMath remMask = errorMask;
    WideMath zeroMask;
    zeroMask.mi = _mm_and_si128(_mm_cmpeq_epi32(xw.mi, _mm_castps_si128(_mm_set_ss(-F32_INF))),
                                remMask.mi);
    remMask.mi = _mm_andnot_si128(zeroMask.mi, remMask.mi);

    WideMath infOutMask;
    infOutMask.mi = _mm_and_si128(_mm_cmplt_epi32(_mm_and_si128(_mm_castps_si128(_mm_set_ss(F32_INF)),
                                                                topMask.mi), abstop.mi), remMask.mi);
    errors.m = _mm_and_ps(infOutMask.m, _mm_add_ps(xw.m, xw.m));
    remMask.mi = _mm_andnot_si128(infOutMask.mi, remMask.mi);

    WideMath overflowMask;
    overflowMask.mi = _mm_and_si128(_mm_castps_si128(_mm_cmpgt_ps(xw.m, _mm_set_ss(0x1.62E42Ep6f))),
                                    remMask.mi);
    errors.m = _mm_or_ps(_mm_and_ps(overflowMask.m, overflow.m), errors.m);
    remMask.mi = _mm_andnot_si128(overflowMask.mi, remMask.mi);

    WideMath underflowMask;
    underflowMask.mi = _mm_and_si128(_mm_castps_si128(_mm_cmplt_ps(xw.m, _mm_set_ss(-0x1.9FE368p6f))),
                                     remMask.mi);
    errors.m = _mm_or_ps(_mm_and_ps(underflowMask.m, underflow.m), errors.m);

    WideMath xd;
    xd.md = _mm_cvtps_pd(xw.m);
    WideMath z;
    z.md = _mm_mul_pd(_mm_set_sd(gExp2F32_InvLn2Scaled), xd.md);
    WideMath kd;
    kd.md = _mm_round_pd(z.md, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    u64 ki = (u64)_mm_cvtsd_si64(kd.md);

    u64 t = gExp2F32_Table[ki % (1 << EXP2F_TABLE_BITS)];
    t += ki << (52 - EXP2F_TABLE_BITS);

    WideMath r;
    r.md = _mm_sub_pd(z.md, kd.md);
    WideMath s;
    s.mi = _mm_set1_epi64x(t);

#define EXP2F_N  ((f64)(1 << EXP2F_TABLE_BITS))
    WideMath poly0;
    poly0.md = _mm_set_sd(0x1.c6af84b912394p-5 / EXP2F_N / EXP2F_N / EXP2F_N);
    WideMath poly1;
    poly1.md = _mm_set_sd(0x1.ebfce50fac4f3p-3 / EXP2F_N / EXP2F_N);
    WideMath poly2;
    poly2.md = _mm_set_sd(0x1.62e42ff0c52d6p-1 / EXP2F_N);
#undef EXP2F_N
    z.md = _mm_add_pd(_mm_mul_pd(poly0.md, r.md), poly1.md);

    WideMath r2;
    r2.md = _mm_mul_pd(r.md, r.md);

    WideMath y;
    y.md = _mm_add_pd(_mm_mul_pd(poly2.md, r.md), _mm_set_sd(1.0));
    y.md = _mm_add_pd(y.md, _mm_mul_pd(z.md, r2.md));
    y.md = _mm_mul_pd(y.md, s.md);

    WideMath result;
    result.m = _mm_or_ps(_mm_andnot_ps(errorMask.m, _mm_cvtpd_ps(y.md)), errors.m);
    return result.e[0];
}

#include "../mats/mats_elem_ext.h"
#include "../mats/mats_elem4x.h"

#include "test_common.cpp"

internal f32_4x
sqrtf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = sqrtf(x.e[0]);
    result.e[1] = sqrtf(x.e[1]);
    result.e[2] = sqrtf(x.e[2]);
    result.e[3] = sqrtf(x.e[3]);
    return result;
}

internal f32_4x
hypotf_4x(f32_4x x, f32_4x y)
{
    f32_4x result;
    result.e[0] = hypotf(x.e[0], y.e[0]);
    result.e[1] = hypotf(x.e[1], y.e[1]);
    result.e[2] = hypotf(x.e[2], y.e[2]);
    result.e[3] = hypotf(x.e[3], y.e[3]);
    return result;
}

internal f32_4x
expf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = expf(x.e[0]);
    result.e[1] = expf(x.e[1]);
    result.e[2] = expf(x.e[2]);
    result.e[3] = expf(x.e[3]);
    return result;
}

internal f32_4x
exp2f_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = exp2f(x.e[0]);
    result.e[1] = exp2f(x.e[1]);
    result.e[2] = exp2f(x.e[2]);
    result.e[3] = exp2f(x.e[3]);
    return result;
}

internal f32_4x
logf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = logf(x.e[0]);
    result.e[1] = logf(x.e[1]);
    result.e[2] = logf(x.e[2]);
    result.e[3] = logf(x.e[3]);
    return result;
}

internal f32_4x
log2f_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = log2f(x.e[0]);
    result.e[1] = log2f(x.e[1]);
    result.e[2] = log2f(x.e[2]);
    result.e[3] = log2f(x.e[3]);
    return result;
}

internal f32_4x
log10f_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = log10f(x.e[0]);
    result.e[1] = log10f(x.e[1]);
    result.e[2] = log10f(x.e[2]);
    result.e[3] = log10f(x.e[3]);
    return result;
}

internal f32_4x
powf_4x(f32_4x a, f32_4x b)
{
    f32_4x result;
    result.e[0] = powf(a.e[0], b.e[0]);
    result.e[1] = powf(a.e[1], b.e[1]);
    result.e[2] = powf(a.e[2], b.e[2]);
    result.e[3] = powf(a.e[3], b.e[3]);
    return result;
}

internal f32_4x
pow32_4x_temp(f32_4x a, f32_4x b)
{
    f32_4x result;
    result.e[0] = pow32(a.e[0], b.e[0]);
    result.e[1] = pow32(a.e[1], b.e[1]);
    result.e[2] = pow32(a.e[2], b.e[2]);
    result.e[3] = pow32(a.e[3], b.e[3]);
    return result;
}

internal f32_4x
log10_32_fast_4x_t(f32_4x x)
{
    return log10_32_fast_4x(x);
}

enum DoTestFlag
{
    DoTest_Sqrt       = 0x00000001,
    DoTest_Hypot      = 0x00000002,
    DoTest_Exp        = 0x00000004,
    DoTest_Exp2       = 0x00000008,
    DoTest_Log        = 0x00000010,
    DoTest_Log2       = 0x00000020,
    DoTest_Log10      = 0x00000040,
    DoTest_Pow        = 0x00000080,
    DoTest_FuncMask   = 0x000000FF,

    DoTest_NoFpBehave = 0x00100000,
    DoTest_SpecMask   = 0x0FF00000,

    DoTest_Comp       = 0x10000000,
    DoTest_Speed      = 0x20000000,
    DoTest_Wide       = 0x40000000,
    DoTest_TypeMask   = 0xF0000000,
};

s32 main(s32 argc, char **argv)
{
    u32 doTests = U32_MAX;

    if (argc > 1)
    {
        u32 tests = DoTest_NoFpBehave;
        for (u32 index = 1; index < (u32)argc; ++index)
        {
            if (argv[index][0] == '-') {
                char *str = argv[index];
                ++str;
                while (*str) {
                    if (*str == 'c') {
                        tests |= DoTest_Comp;
                    } else if (*str == 's') {
                        tests |= DoTest_Speed;
                    } else if (*str == 'w') {
                        tests |= DoTest_Wide;
                    } else {
                        i_expect(*str == 'f');
                        tests &= ~DoTest_NoFpBehave;
                    }
                    ++str;
                }
            } else if (strings_are_equal("sqrt", argv[index])) {
                tests |= DoTest_Sqrt;
            } else if (strings_are_equal("hypot", argv[index])) {
                tests |= DoTest_Hypot;
            } else if (strings_are_equal("log", argv[index])) {
                tests |= DoTest_Log;
            } else if (strings_are_equal("log2", argv[index])) {
                tests |= DoTest_Log2;
            } else if (strings_are_equal("log10", argv[index])) {
                tests |= DoTest_Log10;
            } else if (strings_are_equal("exp", argv[index])) {
                tests |= DoTest_Exp;
            } else if (strings_are_equal("exp2", argv[index])) {
                tests |= DoTest_Exp2;
            } else if (strings_are_equal("pow", argv[index])) {
                tests |= DoTest_Pow;
            }
        }

        if (tests & DoTest_FuncMask) {
            doTests = (doTests & ~DoTest_FuncMask) | tests;
        }
        if (tests & DoTest_SpecMask) {
            doTests = (doTests & ~DoTest_SpecMask) | tests;
        }
        if (tests & DoTest_TypeMask) {
            doTests = (doTests & ~DoTest_TypeMask) | tests;
        }
    }

    if ((doTests & DoTest_NoFpBehave) == 0)
    {
        set_default_fp_behavior();
    }

    u32 tests = 50000000;
    f32 minVal = 1.0e-9f;
    f32 maxVal = 123.8f;

#if 0
    u32 testsA = 100;
    f32 minValA = -1.0f;
    f32 maxValA =  2.0f;
    u32 testsB = 100;
    f32 minValB = -5.0f;
    f32 maxValB =  5.0f;
#else
    u32 testsA = 10000;
    f32 minValA = -10.0f;
    f32 maxValA =  10.0f;
    u32 testsB = 4000;
    f32 minValB = -13.0f;
    f32 maxValB =  14.0f;
#endif

#if 0
    f32 pRes = powf(2.34f, 2.123f);
    f32_4x mRes = pow32_4x(F32_4x(2.34f), F32_4x(2.123f));
    fprintf(stdout, "P: %a %a | %g %g\n", pRes, mRes.e[0], pRes, mRes.e[0]);
    f64 plRes = pow32_log2(u32f32(2.34f).u);
    f64_2x mlResLo, mlResHi;
    pow32_log2_4x(F32_4x(2.34f), &mlResLo, &mlResHi);
    fprintf(stdout, "PL: %a %a | %g %g\n", plRes, mlResLo.e[0], plRes, mlResLo.e[0]);
    f32 peRes = pow32_exp2(plRes, 0);
    f32_4x meRes = pow32_exp2_4x(mlResLo, mlResHi, zero_f32_4x());
    fprintf(stdout, "PE: %a %a | %g %g\n", peRes, meRes.e[0], peRes, meRes.e[0]);
#endif

    fprintf(stdout, "\n");

    if (doTests & DoTest_Comp)
    {
        fprintf(stdout, "Correctness\n\n");

        if (doTests & DoTest_Sqrt)
        {
            fprintf(stdout, "Square root\n");
            f32 stdSec = call_comp_x(stdlib, sqrt, f, 0.0f);
            call_comp_x(mats, sqrt, 32_nonsse, stdSec);
            call_comp_x(matsse, sqrt, 32_sse, stdSec);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Hypot)
        {
            fprintf(stdout, "Hypothenuse\n");
            f32 stdSec = call_comp_x2(stdlib, hypot, f, 0.0f);
            call_comp_x2(mats, hypot, 32_nonsse, stdSec);
            call_comp_x2(mats, hypot, 32_sse, stdSec);
            call_comp_x2(ieee754, hypot, _ieee754, stdSec);
            call_comp_x2_4x(matsse, hypot, 32_4x, stdSec);
            call_comp_x2_4x(fatsse, hypot, 32_fast_4x, stdSec);
            fprintf(stdout, "\n");
        }

        minVal = -12.0f;
        maxVal = 12.0f;

        if (doTests & DoTest_Exp)
        {
            fprintf(stdout, "Exp\n");
            f32 stdSecExp32 = call_comp_x(stdlib, exp, f, 0.0f);
            call_comp_x(mats, exp, 32_nonsse, stdSecExp32);
            call_comp_x(matsse, exp, 32_sse, stdSecExp32);
            call_comp_x_4x(matsse, exp, 32_4x, stdSecExp32);
            call_comp_x_4x(fatsse, exp, 32_fast_4x, stdSecExp32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Exp2)
        {
            fprintf(stdout, "Exp2\n");
            f32 stdSecExp232 = call_comp_x(stdlib, exp2, f, 0.0f);
            call_comp_x(mats, exp2, _32_nonsse, stdSecExp232);
            call_comp_x_4x(matsse, exp2, _32_4x, stdSecExp232);
            call_comp_x_4x(fatsse, exp2, _32_fast_4x, stdSecExp232);
            fprintf(stdout, "\n");
        }

        minVal = 1.0e-9f;
        maxVal = 123.8f;

        if (doTests & DoTest_Log)
        {
            fprintf(stdout, "Log\n");
            f32 stdSecLog32 = call_comp_x(stdlib, log, f, 0.0f);
            call_comp_x(mats, log, 32_nonsse, stdSecLog32);
            call_comp_x_4x(matsse, log, 32_4x, stdSecLog32);
            call_comp_x_4x(fatsse, log, 32_fast_4x, stdSecLog32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Log2)
        {
            fprintf(stdout, "Log2\n");
            f32 stdSecLog232 = call_comp_x(stdlib, log2, f, 0.0f);
            call_comp_x(mats, log2, _32_nonsse, stdSecLog232);
            call_comp_x_4x(matsse, log2, _32_4x, stdSecLog232);
            call_comp_x_4x(fatsse, log2, _32_fast_4x, stdSecLog232);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Log10)
        {
            fprintf(stdout, "Log10\n");
            f32 stdSec = call_comp_x(stdlib, log10, f, 0.0f);
            call_comp_x(mats, log10, _32_nonsse, stdSec);
            call_comp_x_4x(matsse, log10, _32_4x, stdSec);
            call_comp_x_4x(fatsse, log10, _32_fast_4x_t, stdSec);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Pow)
        {
            fprintf(stdout, "Pow\n");
            f32 stdSecPow32 = call_comp_x2(stdlib, pow, f, 0.0f);
            call_comp_x2(mats, pow, 32, stdSecPow32);
            call_comp_x2_4x(matsse, pow, 32_4x, stdSecPow32);
            call_comp_x2_4x(mats, pow, 32_4x_temp, stdSecPow32);
            fprintf(stdout, "\n");
        }
    }

    minVal = 1.0e-9f;
    maxVal = 123.8f;

    if (doTests & DoTest_Speed)
    {
        fprintf(stdout, "Speed\n\n");

        if (doTests & DoTest_Sqrt)
        {
            fprintf(stdout, "Square root\n");
            f32 stdSec = call_spd(stdlib, sqrt, f, 0.0f);
            call_spd(mats, sqrt, 32_nonsse, stdSec);
            call_spd(mats, sqrt, 32_sse, stdSec);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Hypot)
        {
            fprintf(stdout, "Hypothenuse\n");
            f32 stdSec = call_spd2(stdlib, hypot, f, 0.0f);
            call_spd2(mats, hypot, 32_nonsse, stdSec);
            call_spd2(mats, hypot, 32_sse, stdSec);
            call_spd2(ieee754, hypot, _ieee754, stdSec);
            fprintf(stdout, "\n");
        }

        minVal = -12.0f;
        maxVal = 12.0f;

        if (doTests & DoTest_Exp)
        {
            fprintf(stdout, "Exp\n");
            f32 stdSec = call_spd(stdlib, exp, f, 0.0f);
            call_spd(mats, exp, 32_nonsse, stdSec);
            call_spd(matsse, exp, 32_sse, stdSec);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Exp2)
        {
            fprintf(stdout, "Exp2\n");
            f32 stdSec = call_spd(stdlib, exp2, f, 0.0f);
            call_spd(mats, exp2, _32_nonsse, stdSec);
            fprintf(stdout, "\n");
        }

        minVal = 1.0e-9f;
        maxVal = 123.8f;

        if (doTests & DoTest_Log)
        {
            fprintf(stdout, "Log\n");
            f32 stdSec = call_spd(stdlib, log, f, 0.0f);
            call_spd(mats, log, 32_nonsse, stdSec);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Log2)
        {
            fprintf(stdout, "Log2\n");
            f32 stdSec = call_spd(stdlib, log2, f, 0.0f);
            call_spd(mats, log2, _32_nonsse, stdSec);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Log10)
        {
            fprintf(stdout, "Log10\n");
            f32 stdSec = call_spd(stdlib, log10, f, 0.0f);
            call_spd(mats, log10, _32_nonsse, stdSec);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Pow)
        {
            fprintf(stdout, "Pow\n");
            f32 spdSecPow32 = call_spd2(stdlib, pow, f, 0.0f);
            call_spd2(mats, pow, 32, spdSecPow32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Wide)
        {
            tests = 5000000;

            minVal = 1.0e-9f;
            maxVal = 123.8f;

            if (doTests & DoTest_Sqrt)
            {
                fprintf(stdout, "Sqrt wide\n");
                f32 stdSec = call_spd_4x(stdlib, sqrt, f_4x, 0.0f);
                call_spd_4x(mats, sqrt, 32_4x, stdSec);
                fprintf(stdout, "\n");
            }

            if (doTests & DoTest_Hypot)
            {
                fprintf(stdout, "Hypot wide\n");
                f32 stdSec = call_spd2_4x(stdlib, hypot, f_4x, 0.0f);
                call_spd2_4x(mats, hypot, 32_4x, stdSec);
                call_spd2_4x(fats, hypot, 32_fast_4x, stdSec);
                fprintf(stdout, "\n");
            }

            minVal = -12.0f;
            maxVal = 12.0f;

            if (doTests & DoTest_Exp)
            {
                fprintf(stdout, "Exp wide\n");
                f32 stdSec = call_spd_4x(stdlib, exp, f_4x, 0.0f);
                call_spd_4x(mats, exp, 32_4x, stdSec);
                call_spd_4x(fats, exp, 32_fast_4x, stdSec);
                fprintf(stdout, "\n");
            }

            if (doTests & DoTest_Exp2)
            {
                fprintf(stdout, "Exp2 wide\n");
                f32 stdSec = call_spd_4x(stdlib, exp2, f_4x, 0.0f);
                call_spd_4x(mats, exp2, _32_4x, stdSec);
                call_spd_4x(fats, exp2, _32_fast_4x, stdSec);
                fprintf(stdout, "\n");
            }

            minVal = 1.0e-9f;
            maxVal = 123.8f;

            if (doTests & DoTest_Log)
            {
                fprintf(stdout, "Log wide\n");
                f32 stdSec = call_spd_4x(stdlib, log, f_4x, 0.0f);
                call_spd_4x(mats, log, 32_4x, stdSec);
                call_spd_4x(fats, log, 32_fast_4x, stdSec);
                fprintf(stdout, "\n");
            }

            if (doTests & DoTest_Log2)
            {
                fprintf(stdout, "Log2 wide\n");
                f32 stdSec = call_spd_4x(stdlib, log2, f_4x, 0.0f);
                call_spd_4x(mats, log2, _32_4x, stdSec);
                call_spd_4x(fats, log2, _32_fast_4x, stdSec);
                fprintf(stdout, "\n");
            }

            if (doTests & DoTest_Log10)
            {
                fprintf(stdout, "Log10 wide\n");
                f32 stdSec = call_spd_4x(stdlib, log10, f_4x, 0.0f);
                call_spd_4x(mats, log10, _32_4x, stdSec);
                call_spd_4x(fats, log10, _32_fast_4x_t, stdSec);
                fprintf(stdout, "\n");
            }

            if (doTests & DoTest_Pow)
            {
                fprintf(stdout, "Pow wide\n");
                f32 spdSecPow32 = call_spd2_4x(stdlib, pow, f_4x, 0.0f);
                call_spd2_4x(mats, pow, 32_4x, spdSecPow32);
                call_spd2_4x(tats, pow, 32_4x_temp, spdSecPow32);
                fprintf(stdout, "\n");
            }

        }
    }

    return 0;
}
