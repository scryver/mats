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
#define exp32    exp32_nonsse
#define exp2_32  exp2_32_nonsse
#define log32    log32_nonsse
#define log2_32  log2_32_nonsse
#define MATS_USE_SSE2 0
#include "../mats/mats_defines.h"
#include "../mats/mats_elem.h"
#undef sqrt32
#undef exp32
#undef exp2_32
#undef log32
#undef log2_32

#define sqrt32   sqrt32_sse
#define exp32    exp32_not_used
#define exp2_32  exp2_32_not_used
#define log32    log32_not_used
#define log2_32  log2_32_not_used
#undef  MATS_USE_SSE2
#undef  MATS_USE_SSE4
#define MATS_USE_SSE2 1
#define MATS_USE_SSE4 1
#include "../mats/mats_defines.h"
#include "../mats/mats_elem.h"
#undef sqrt32
#undef exp32
#undef exp2_32
#undef log32
#undef log2_32

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

enum DoTestFlag
{
    DoTest_Sqrt       = 0x00000001,
    DoTest_Log        = 0x00000002,
    DoTest_Log2       = 0x00000004,
    DoTest_Exp        = 0x00000008,
    DoTest_Exp2       = 0x00000010,
    DoTest_Pow        = 0x00000020,
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
            } else if (strings_are_equal("log", argv[index])) {
                tests |= DoTest_Log;
            } else if (strings_are_equal("log2", argv[index])) {
                tests |= DoTest_Log2;
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

    u32 testsA = 10000;
    f32 minValA = 1.0e-1f;
    f32 maxValA = 12.8f;
    u32 testsB = 10000;
    f32 minValB = -10.0f;
    f32 maxValB =  10.0f;

    fprintf(stdout, "\n");

    if (doTests & DoTest_Comp)
    {
        fprintf(stdout, "Correctness\n\n");

        if (doTests & DoTest_Sqrt)
        {
            fprintf(stdout, "Square root\n");
            f32 stdSecSqrt32 = run_comp_f32_x(string("stdlib sqrt"), 15, "sqrt", tests, minVal, maxVal, sqrt, sqrtf, 0.0f);
            run_comp_f32_x(string("mats sqrt"), 15, "sqrt", tests, minVal, maxVal, sqrt, sqrt32_nonsse, stdSecSqrt32);
            run_comp_f32_x(string("mats sqrt sse"), 15, "sqrt", tests, minVal, maxVal, sqrt, sqrt32_sse, stdSecSqrt32);
            fprintf(stdout, "\n");
        }

        minVal = -12.0f;
        maxVal = 12.0f;

        if (doTests & DoTest_Exp)
        {
            fprintf(stdout, "Exp\n");
            f32 stdSecExp32 = run_comp_f32_x(string("stdlib exp"), 15, "exp", tests, minVal, maxVal, exp, expf, 0.0f);
            run_comp_f32_x(string("mats exp"), 15, "exp", tests, minVal, maxVal, exp, exp32_nonsse, stdSecExp32);
            run_comp_f32_x(string("mats exp sse"), 15, "exp", tests, minVal, maxVal, exp, exp32_sse, stdSecExp32);
            run_comp_f32_4x_x(string("mats exp 4x"), 15, "exp", tests, minVal, maxVal, exp, exp32_4x, stdSecExp32);
            run_comp_f32_4x_x(string("mats fast 4x"), 15, "exp", tests, minVal, maxVal, exp, exp32_fast_4x, stdSecExp32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Exp2)
        {
            fprintf(stdout, "Exp2\n");
            f32 stdSecExp232 = run_comp_f32_x(string("stdlib exp2"), 15, "exp2", tests, minVal, maxVal, exp2, exp2f, 0.0f);
            run_comp_f32_x(string("mats exp2"), 15, "exp2", tests, minVal, maxVal, exp2, exp2_32_nonsse, stdSecExp232);
            run_comp_f32_4x_x(string("mats exp2 4x"), 15, "exp2", tests, minVal, maxVal, exp2, exp2_32_4x, stdSecExp232);
            run_comp_f32_4x_x(string("mats fast 4x"), 15, "exp2", tests, minVal, maxVal, exp2, exp2_32_fast_4x, stdSecExp232);
            fprintf(stdout, "\n");
        }

        minVal = 1.0e-9f;
        maxVal = 123.8f;

        if (doTests & DoTest_Log)
        {
            fprintf(stdout, "Log\n");
            f32 stdSecLog32 = run_comp_f32_x(string("stdlib log"), 15, "log", tests, minVal, maxVal, log, logf, 0.0f);
            run_comp_f32_x(string("mats log"), 15, "log", tests, minVal, maxVal, log, log32_nonsse, stdSecLog32);
            run_comp_f32_4x_x(string("mats log 4x"), 15, "log", tests, minVal, maxVal, log, log32_4x, stdSecLog32);
            run_comp_f32_4x_x(string("mats fast 4x"), 15, "log", tests, minVal, maxVal, log, log32_fast_4x, stdSecLog32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Log2)
        {
            fprintf(stdout, "Log2\n");
            f32 stdSecLog232 = run_comp_f32_x(string("stdlib log2"), 15, "log2", tests, minVal, maxVal, log2, log2f, 0.0f);
            run_comp_f32_x(string("mats log2"), 15, "log2", tests, minVal, maxVal, log2, log2_32_nonsse, stdSecLog232);
            run_comp_f32_4x_x(string("mats log2 4x"), 15, "log2", tests, minVal, maxVal, log2, log2_32_4x, stdSecLog232);
            run_comp_f32_4x_x(string("mats fast 4x"), 15, "log2", tests, minVal, maxVal, log2, log2_32_fast_4x, stdSecLog232);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Pow)
        {
            fprintf(stdout, "Pow\n");
            f32 stdSecPow32 = call_comp_x2(stdlib, pow, f, 0.0f);
            call_comp_x2(mats, pow, 32, stdSecPow32);
            //call_comp_x2(mats, pow, 32_temp, stdSecPow32);
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
            f32 spdSecSqrt32 = run_speed_f32(string("stdlib sqrt"), 15, "sqrt", tests, minVal, maxVal, sqrtf, 0.0f);
            run_speed_f32(string("mats sqrt"), 15, "sqrt", tests, minVal, maxVal, sqrt32_nonsse, spdSecSqrt32);
            run_speed_f32(string("mats sqrt sse"), 15, "sqrt", tests, minVal, maxVal, sqrt32_sse, spdSecSqrt32);
            fprintf(stdout, "\n");
        }

        minVal = -12.0f;
        maxVal = 12.0f;

        if (doTests & DoTest_Exp)
        {
            fprintf(stdout, "Exp\n");
            f32 spdSecExp32 = run_speed_f32(string("stdlib exp"), 15, "exp", tests, minVal, maxVal, expf, 0.0f);
            run_speed_f32(string("mats exp"), 15, "exp", tests, minVal, maxVal, exp32_nonsse, spdSecExp32);
            run_speed_f32(string("mats exp sse"), 15, "exp", tests, minVal, maxVal, exp32_sse, spdSecExp32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Exp2)
        {
            fprintf(stdout, "Exp2\n");
            f32 spdSecExp232 = run_speed_f32(string("stdlib exp2"), 15, "exp2", tests, minVal, maxVal, exp2f, 0.0f);
            run_speed_f32(string("mats exp2"), 15, "exp2", tests, minVal, maxVal, exp2_32_nonsse, spdSecExp232);
            fprintf(stdout, "\n");
        }

        minVal = 1.0e-9f;
        maxVal = 123.8f;

        if (doTests & DoTest_Log)
        {
            fprintf(stdout, "Log\n");
            f32 spdSecLog32 = run_speed_f32(string("stdlib log"), 15, "log", tests, minVal, maxVal, logf, 0.0f);
            run_speed_f32(string("mats log"), 15, "log", tests, minVal, maxVal, log32_nonsse, spdSecLog32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Log2)
        {
            fprintf(stdout, "Log2\n");
            f32 spdSecLog232 = run_speed_f32(string("stdlib log2"), 15, "log2", tests, minVal, maxVal, log2f, 0.0f);
            run_speed_f32(string("mats log2"), 15, "log2", tests, minVal, maxVal, log2_32_nonsse, spdSecLog232);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Pow)
        {
            fprintf(stdout, "Pow\n");
            f32 spdSecPow32 = call_spd2(stdlib, pow, f, 0.0f);
            call_spd2(mats, pow, 32, spdSecPow32);
            //call_spd2(mats, pow, 32_temp, spdSecPow32);
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
                f32 spdSecSqrt4x = run_speed_f32_4x(string("stdlib sqrt 4x"), 15, "sqrt_4x", tests, minVal, maxVal, sqrtf_4x, 0.0f);
                run_speed_f32_4x(string("mats sqrt 4x"), 15, "sqrt_4x", tests, minVal, maxVal, sqrt32_4x, spdSecSqrt4x);
                fprintf(stdout, "\n");
            }

            minVal = -12.0f;
            maxVal = 12.0f;

            if (doTests & DoTest_Exp)
            {
                fprintf(stdout, "Exp wide\n");
                f32 spdSecExp4x = run_speed_f32_4x(string("stdlib exp 4x"), 15, "exp_4x", tests, minVal, maxVal, expf_4x, 0.0f);
                run_speed_f32_4x(string("mats exp 4x"), 15, "exp_4x", tests, minVal, maxVal, exp32_4x, spdSecExp4x);
                run_speed_f32_4x(string("mats fast 4x"), 15, "exp_4x", tests, minVal, maxVal, exp32_fast_4x, spdSecExp4x);
                fprintf(stdout, "\n");
            }

            if (doTests & DoTest_Exp2)
            {
                fprintf(stdout, "Exp2 wide\n");
                f32 spdSecExp4x = run_speed_f32_4x(string("stdlib exp2 4x"), 15, "exp2_4x", tests, minVal, maxVal, exp2f_4x, 0.0f);
                run_speed_f32_4x(string("mats exp2 4x"), 15, "exp2_4x", tests, minVal, maxVal, exp2_32_4x, spdSecExp4x);
                run_speed_f32_4x(string("mats fast 4x"), 15, "exp2_4x", tests, minVal, maxVal, exp2_32_fast_4x, spdSecExp4x);
                fprintf(stdout, "\n");
            }

            minVal = 1.0e-9f;
            maxVal = 123.8f;

            if (doTests & DoTest_Log)
            {
                fprintf(stdout, "Log wide\n");
                f32 spdSecLog32 = run_speed_f32_4x(string("stdlib log 4x"), 15, "log", tests, minVal, maxVal, logf_4x, 0.0f);
                run_speed_f32_4x(string("mats log 4x"), 15, "log", tests, minVal, maxVal, log32_4x, spdSecLog32);
                run_speed_f32_4x(string("mats fast 4x"), 15, "log", tests, minVal, maxVal, log32_fast_4x, spdSecLog32);
                fprintf(stdout, "\n");
            }

            if (doTests & DoTest_Log2)
            {
                fprintf(stdout, "Log2 wide\n");
                f32 spdSecLog232 = run_speed_f32_4x(string("stdlib log2 4x"), 15, "log2", tests, minVal, maxVal, log2f_4x, 0.0f);
                run_speed_f32_4x(string("mats log2 4x"), 15, "log2", tests, minVal, maxVal, log2_32_4x, spdSecLog232);
                run_speed_f32_4x(string("mats fast 4x"), 15, "log2", tests, minVal, maxVal, log2_32_fast_4x, spdSecLog232);
            }
        }
    }

    return 0;
}
