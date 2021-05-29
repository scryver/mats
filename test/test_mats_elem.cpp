#include <errno.h>
#include <time.h>
#include <math.h>
#include <x86intrin.h>
#include <xmmintrin.h>

#define STR_FMT(x)   safe_truncate_to_s32(x.size), (char *)x.data
#include "../libberdip/src/common.h"
#include "../libberdip/src/multilane.h"

#define MATS_F32_ABS_MASK  0x7FFFFFFF
#define MATS_F64_ABS_MASK  0x7FFFFFFFFFFFFFFFULL
#define IEEE_754_2008_SNAN 1
#include "../mats/mats_common.h"
#include "../mats/mats_constants.h"
#include "../mats/mats4x.h"

#define sqrt32   sqrt32_nonsse
#define hypot32  hypot32_nonsse
#define exp32    exp32_nonsse
#define exp2_32  exp2_32_nonsse
#define pow2_32  pow2_32_nonsse
#define log32    log32_nonsse
#define log2_32  log2_32_nonsse
#define log10_32 log10_32_nonsse
#define expm1_32 expm1_32_nonsse
#define log1p32  log1p32_nonsse
#define log1p_fast32  log1p_fast32_nonsse
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
#undef expm1_32
#undef log1p32
#undef log1p_fast32

#define sqrt32   sqrt32_sse
#define hypot32  hypot32_sse
#define exp32    exp32_not_used
#define exp2_32  exp2_32_not_used
#define log32    log32_not_used
#define log2_32  log2_32_not_used
#define pow2_32  pow2_32_not_used
#define log10_32 log10_32_not_used
#define expm1_32 expm1_32_not_used
#define log1p32  log1p32_not_used
#define log1p_fast32  log1p_fast32_not_used
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
#undef expm1_32
#undef log1p32
#undef log1p_fast32

internal f32
hypot_ieee754(f32 x, f32 y)
{
	s32 ha = MATS_S32_FROM_F32(x) & MATS_F32_ABS_MASK;
    s32 hb = MATS_S32_FROM_F32(y) & MATS_F32_ABS_MASK;

    if (hb > ha) {
        s32 temp = ha;
        ha = hb;
        hb = temp;
    }

    f32 a = MATS_F32_FROM_S32(ha);
    f32 b = MATS_F32_FROM_S32(hb);

    if ((ha - hb) > 0x0F000000) {
        return a + b; /* x/y > 2**30 */
    }

    s32 k = 0;
	if (ha > 0x58800000)
    {   /* a > 2**50 */
        if (!MATS_F32_UWORD_IS_FINITE(ha))
        {
            f32 w = a + b;
            if (MATS_F32_UWORD_IS_INFINITE(ha)) {
                w = a;
            } else if (MATS_F32_UWORD_IS_INFINITE(hb)) {
                w = b;
            }
            return w;
        }
        /* scale a and b by 2**-68 */
        ha -= 0x22000000;
        hb -= 0x22000000;
        k  += 68;
        a = MATS_F32_FROM_S32(ha);
        b = MATS_F32_FROM_S32(hb);
	}

    if (hb < 0x26800000L) {
        /* b < 2**-50 */
        if (MATS_F32_UWORD_IS_ZERO(hb)) {
            return a;
        } else if (MATS_F32_UWORD_IS_SUBNORMAL(hb)) {
            f32 t = MATS_F32_FROM_U32(0x7E800000); /* t = 2^126 */
            a *= t;
            b *= t;
            k -= 126;
        } else { /* scale a and b by 2^68 */
            ha += 0x22000000; /* a *= 2^68 */
            hb += 0x22000000; /* b *= 2^68 */
            k  -= 68;
            a   = MATS_F32_FROM_S32(ha);
            b   = MATS_F32_FROM_S32(hb);
        }
	}

    /* medium size a and b */
    f32 w = a - b;
    if (w > b)
    {
        f32 t1 = MATS_F32_FROM_U32(ha & 0xFFFFF000);
        f32 t2 = a - t1;
        w = sqrt32_sse(t1 * t1 - (b * (-b) - t2 * (a + t1)));
    }
    else
    {
        a = a + a;
        f32 y1 = MATS_F32_FROM_U32(hb & 0xFFFFF000);
        f32 y2 = b - y1;
        f32 t1 = MATS_F32_FROM_U32((ha + 0x00800000) & 0xFFFFF000);
        f32 t2 = a - t1;
        w = sqrt32_sse(t1 * y1 - (w * (-w) - (t1 * y2 + t2 * b)));
    }

    if (k != 0)
    {
        f32 t1 = MATS_F32_FROM_U32((u32)0x3F800000 + (k << MATS_F32_EXP_SHIFT));
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

WIDE_FUNC_FROM_F32(sqrtf)
WIDE_FUNC_FROM_F32_F32(hypotf)
WIDE_FUNC_FROM_F32(expf)
WIDE_FUNC_FROM_F32(exp2f)
WIDE_FUNC_FROM_F32(expm1f)
WIDE_FUNC_FROM_F32(logf)
WIDE_FUNC_FROM_F32(log2f)
WIDE_FUNC_FROM_F32(log10f)
WIDE_FUNC_FROM_F32(log1pf)
WIDE_FUNC_FROM_F32_F32(powf)

#define pow32_temp pow32
WIDE_FUNC_FROM_F32_F32(pow32_temp)

internal f32_4x
log10_32_fast_4x_t(f32_4x x)
{
    return log10_32_fast_4x(x);
}

#define expm1_temp32 expm1_32_nonsse
WIDE_FUNC_FROM_F32(expm1_temp32)

enum DoTestFlag
{
    DoTest_sqrt       = 0x00000001,
    DoTest_hypot      = 0x00000002,
    DoTest_exp        = 0x00000004,
    DoTest_exp2       = 0x00000008,
    DoTest_log        = 0x00000010,
    DoTest_log2       = 0x00000020,
    DoTest_log10      = 0x00000040,
    DoTest_expm1      = 0x00000080,
    DoTest_log1p      = 0x00000100,
    DoTest_pow        = 0x00000200,
    DoTest_FuncMask   = 0x00000FFF,

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
                tests |= DoTest_sqrt;
            } else if (strings_are_equal("hypot", argv[index])) {
                tests |= DoTest_hypot;
            } else if (strings_are_equal("exp", argv[index])) {
                tests |= DoTest_exp;
            } else if (strings_are_equal("exp2", argv[index])) {
                tests |= DoTest_exp2;
            } else if (strings_are_equal("log", argv[index])) {
                tests |= DoTest_log;
            } else if (strings_are_equal("log2", argv[index])) {
                tests |= DoTest_log2;
            } else if (strings_are_equal("log10", argv[index])) {
                tests |= DoTest_log10;
            } else if (strings_are_equal("expm1", argv[index])) {
                tests |= DoTest_expm1;
            } else if (strings_are_equal("log1p", argv[index])) {
                tests |= DoTest_log1p;
            } else if (strings_are_equal("pow", argv[index])) {
                tests |= DoTest_pow;
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
    f64 plRes = pow32_log2(MATS_U32_FROM_F32(2.34f));
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

        BEGIN_TEST(doTests, sqrt, call_comp);
        call_comp(mats, sqrt, 32_nonsse, stdSec);
        call_comp(matsse, sqrt, 32_sse, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, hypot, call_comp2);
        call_comp2(mats, hypot, 32_nonsse, stdSec);
        call_comp2(mats, hypot, 32_sse, stdSec);
        call_comp2(ieee754, hypot, _ieee754, stdSec);
        call_comp2_4x(matsse, hypot, 32_4x, stdSec);
        call_comp2_4x(fatsse, hypot, 32_fast_4x, stdSec);
        END_TEST();

        minVal = -12.0f;
        maxVal = 12.0f;

        BEGIN_TEST(doTests, exp, call_comp);
        call_comp(mats, exp, 32_nonsse, stdSec);
        call_comp(matsse, exp, 32_sse, stdSec);
        call_comp_4x(matsse, exp, 32_4x, stdSec);
        call_comp_4x(fatsse, exp, 32_fast_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, exp2, call_comp);
        call_comp(mats, exp2, _32_nonsse, stdSec);
        call_comp_4x(matsse, exp2, _32_4x, stdSec);
        call_comp_4x(fatsse, exp2, _32_fast_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, expm1, call_comp);
        call_comp(mats, expm1, _32_nonsse, stdSec);
        call_comp_4x(matsse, expm1, _32_4x, stdSec);
        //call_comp_4x(fatsse, expm1, _32_fast_4x_t, stdSec);
        END_TEST();

        minVal = 1.0e-9f;
        maxVal = 123.8f;

        BEGIN_TEST(doTests, log, call_comp);
        call_comp(mats, log, 32_nonsse, stdSec);
        call_comp_4x(matsse, log, 32_4x, stdSec);
        call_comp_4x(fatsse, log, 32_fast_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, log2, call_comp);
        call_comp(mats, log2, _32_nonsse, stdSec);
        call_comp_4x(matsse, log2, _32_4x, stdSec);
        call_comp_4x(fatsse, log2, _32_fast_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, log10, call_comp);
        call_comp(mats, log10, _32_nonsse, stdSec);
        call_comp_4x(matsse, log10, _32_4x, stdSec);
        call_comp_4x(fatsse, log10, _32_fast_4x_t, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, log1p, call_comp);
        call_comp(mats, log1p, 32_nonsse, stdSec);
        call_comp(mats, log1p, _fast32_nonsse, stdSec);
        call_comp_4x(mats4, log1p, 32_4x, stdSec);
        //call_comp_4x(fatsse, log1p, _32_fast_4x_t, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, pow, call_comp2);
        call_comp2(mats, pow, 32, stdSec);
        call_comp2_4x(matsse, pow, 32_4x, stdSec);
        call_comp2_4x(mats, pow, 32_temp_4x, stdSec);
        END_TEST();
    }

    minVal = 1.0e-9f;
    maxVal = 123.8f;

    if (doTests & DoTest_Speed)
    {
        fprintf(stdout, "Speed\n\n");

        BEGIN_TEST(doTests, sqrt, call_spd);
        call_spd(mats, sqrt, 32_nonsse, stdSec);
        call_spd(matsse, sqrt, 32_sse, stdSec);
        call_spd_4x(mats4, sqrt, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, hypot, call_spd2);
        call_spd2(mats, hypot, 32_nonsse, stdSec);
        call_spd2(matsse, hypot, 32_sse, stdSec);
        call_spd2(ieee754, hypot, _ieee754, stdSec);
        call_spd2_4x(mats4, hypot, 32_4x, stdSec);
        call_spd2_4x(mats4, hypot, 32_fast_4x, stdSec);
        END_TEST();

        minVal = -12.0f;
        maxVal = 12.0f;

        BEGIN_TEST(doTests, exp, call_spd);
        call_spd(mats, exp, 32_nonsse, stdSec);
        call_spd(matsse, exp, 32_sse, stdSec);
        call_spd_4x(mats4, exp, 32_4x, stdSec);
        call_spd_4x(mats4, exp, 32_fast_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, exp2, call_spd);
        call_spd(mats, exp2, _32_nonsse, stdSec);
        call_spd_4x(mats4, exp2, _32_4x, stdSec);
        call_spd_4x(mats4, exp2, _32_fast_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, expm1, call_spd);
        call_spd(mats, expm1, _32_nonsse, stdSec);
        call_spd_4x(mats, expm1, _32_4x, stdSec);
        END_TEST();

        minVal = 1.0e-9f;
        maxVal = 123.8f;

        BEGIN_TEST(doTests, log, call_spd);
        call_spd(mats, log, 32_nonsse, stdSec);
        call_spd_4x(mats4, log, 32_4x, stdSec);
        call_spd_4x(mats4, log, 32_fast_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, log2, call_spd);
        call_spd(mats, log2, _32_nonsse, stdSec);
        call_spd_4x(mats4, log2, _32_4x, stdSec);
        call_spd_4x(mats4, log2, _32_fast_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, log10, call_spd);
        call_spd(mats, log10, _32_nonsse, stdSec);
        call_spd_4x(mats4, log10, _32_4x, stdSec);
        call_spd_4x(mats4, log10, _32_fast_4x_t, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, log1p, call_spd);
        call_spd(mats, log1p, 32_nonsse, stdSec);
        call_spd(mats, log1p, _fast32_nonsse, stdSec);
        call_spd_4x(mats4, log1p, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, pow, call_spd2);
        call_spd2(mats, pow, 32, stdSec);
        call_spd2_4x(mats4, pow, 32_4x, stdSec);
        END_TEST();
    }

    if (doTests & DoTest_Wide)
    {
        tests = 5000000;

        minVal = 1.0e-9f;
        maxVal = 123.8f;

        BEGIN_TEST_WIDE(doTests, sqrt, call_spd_4x);
        call_spd_4x(mats, sqrt, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, hypot, call_spd2_4x);
        call_spd2_4x(mats, hypot, 32_4x, stdSec);
        call_spd2_4x(fats, hypot, 32_fast_4x, stdSec);
        END_TEST();

        minVal = -12.0f;
        maxVal = 12.0f;

        BEGIN_TEST_WIDE(doTests, exp, call_spd_4x);
        call_spd_4x(mats, exp, 32_4x, stdSec);
        call_spd_4x(fats, exp, 32_fast_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, exp2, call_spd_4x);
        call_spd_4x(mats, exp2, _32_4x, stdSec);
        call_spd_4x(fats, exp2, _32_fast_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, expm1, call_spd_4x);
        call_spd_4x(mats, expm1, _32_4x, stdSec);
        call_spd_4x(matst, expm1, _temp32_4x, stdSec);
        END_TEST();

        minVal = 1.0e-9f;
        maxVal = 123.8f;

        BEGIN_TEST_WIDE(doTests, log, call_spd_4x);
        call_spd_4x(mats, log, 32_4x, stdSec);
        call_spd_4x(fats, log, 32_fast_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, log2, call_spd_4x);
        call_spd_4x(mats, log2, _32_4x, stdSec);
        call_spd_4x(fats, log2, _32_fast_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, log10, call_spd_4x);
        call_spd_4x(mats, log10, _32_4x, stdSec);
        call_spd_4x(fats, log10, _32_fast_4x_t, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, log1p, call_spd_4x);
        call_spd_4x(mats, log1p, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, pow, call_spd2_4x);
        call_spd2_4x(mats, pow, 32_4x, stdSec);
        call_spd2_4x(tats, pow, 32_temp_4x, stdSec);
        END_TEST();
    }

    return 0;
}
