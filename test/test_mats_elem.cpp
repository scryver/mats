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
#define MATH_FUNC_F32_4x_FROM_F32_4x(name)     f32_4x name(f32_4x x)
typedef MATH_FUNC_F32_4x_FROM_F32_4x(MathFuncF32_4xFromF32_4x);

global volatile f32_4x gRunSpeedSum4x;
internal f32
run_speed_f32_4x(String name, u32 maxNameSize, char *func, u32 tests, f32 minVal, f32 maxVal,
                 MathFuncF32_4xFromF32_4x *testFunc, f32 secondsBase)
{
    gRunSpeedSum4x.m = _mm_setzero_ps();
    f32 oneOverTests = 1.0f / (f32)tests;
    f32 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f32_4x inputVal = F32_4x((f32)index * scale + minVal, 0.0f, F32_MAX, -F32_MAX);
        f32_4x testRes0 = testFunc(inputVal);
        f32_4x testRes1 = testFunc(testRes0);
        f32_4x testRes2 = testFunc(testRes1);
        f32_4x testRes3 = testFunc(testRes2);
        gRunSpeedSum4x.m = _mm_add_ps(gRunSpeedSum4x.m, testRes3.m);
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    if (secondsBase == 0.0f) {
        print_speed_info(name, maxNameSize, func, 16*tests, seconds, gRunSpeedSum4x.e[0], 1.0f);
    } else {
        print_speed_info(name, maxNameSize, func, 16*tests, seconds, gRunSpeedSum4x.e[0], seconds / secondsBase);
    }
    return seconds;
}

internal f32_4x
expf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = expf(x.e[0]);
    result.e[1] = expf(x.e[1]);
    result.e[2] = expf(x.e[2]);
    result.e[3] = expf(x.e[3]);
    result.e[0] = expf(result.e[0]);
    result.e[1] = expf(result.e[1]);
    result.e[2] = expf(result.e[2]);
    result.e[3] = expf(result.e[3]);
    result.e[0] = expf(result.e[0]);
    result.e[1] = expf(result.e[1]);
    result.e[2] = expf(result.e[2]);
    result.e[3] = expf(result.e[3]);
    result.e[0] = expf(result.e[0]);
    result.e[1] = expf(result.e[1]);
    result.e[2] = expf(result.e[2]);
    result.e[3] = expf(result.e[3]);
    return result;
}

s32 main(s32 argc, char **argv)
{
    u32 tests = 50000000;
    f32 minVal = 1.0e-9f;
    f32 maxVal = 123.8f;

    fprintf(stdout, "Pow test\n");
    fprintf(stdout, "1.2^-2.3: %8a %8a | %g %g\n", powf(1.2f, -2.3f), pow32(1.2f, -2.3f), powf(1.2f, -2.3f), pow32(1.2f, -2.3f));
    fprintf(stdout, "1.2^2.3 : %8a %8a | %g %g\n", powf(1.2f, 2.3f), pow32(1.2f, 2.3f), powf(1.2f, 2.3f), pow32(1.2f, 2.3f));

    fprintf(stdout, "\nCorrectness\n\n");

    fprintf(stdout, "Square root\n");
    f32 stdSecSqrt32 = run_comp_f32_x(string("stdlib sqrt"), 15, "sqrt", tests, minVal, maxVal, sqrt, sqrtf, 0.0f);
    run_comp_f32_x(string("mats sqrt"), 15, "sqrt", tests, minVal, maxVal, sqrt, sqrt32_nonsse, stdSecSqrt32);
    run_comp_f32_x(string("mats sqrt sse"), 15, "sqrt", tests, minVal, maxVal, sqrt, sqrt32_sse, stdSecSqrt32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Log\n");
    f32 stdSecLog32 = run_comp_f32_x(string("stdlib log"), 15, "log", tests, minVal, maxVal, log, logf, 0.0f);
    run_comp_f32_x(string("mats log"), 15, "log", tests, minVal, maxVal, log, log32_nonsse, stdSecLog32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Log2\n");
    f32 stdSecLog232 = run_comp_f32_x(string("stdlib log2"), 15, "log2", tests, minVal, maxVal, log2, log2f, 0.0f);
    run_comp_f32_x(string("mats log2"), 15, "log2", tests, minVal, maxVal, log2, log2_32_nonsse, stdSecLog232);
    fprintf(stdout, "\n");

    minVal = -12.0f;
    maxVal = 12.0f;
    fprintf(stdout, "Exp\n");
    f32 stdSecExp32 = run_comp_f32_x(string("stdlib exp"), 15, "exp", tests, minVal, maxVal, exp, expf, 0.0f);
    run_comp_f32_x(string("mats exp"), 15, "exp", tests, minVal, maxVal, exp, exp32_nonsse, stdSecExp32);
    run_comp_f32_x(string("mats exp sse"), 15, "exp", tests, minVal, maxVal, exp, exp32_sse, stdSecExp32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Exp2\n");
    f32 stdSecExp232 = run_comp_f32_x(string("stdlib exp2"), 15, "exp2", tests, minVal, maxVal, exp2, exp2f, 0.0f);
    run_comp_f32_x(string("mats exp2"), 15, "exp2", tests, minVal, maxVal, exp2, exp2_32_nonsse, stdSecExp232);
    fprintf(stdout, "\n");

    fprintf(stdout, "Speed\n\n");

    minVal = 1.0e-9f;
    maxVal = 123.8f;
    fprintf(stdout, "Square root\n");
    f32 spdSecSqrt32 = run_speed_f32(string("stdlib sqrt"), 15, "sqrt", tests, minVal, maxVal, sqrtf, 0.0f);
    run_speed_f32(string("mats sqrt"), 15, "sqrt", tests, minVal, maxVal, sqrt32_nonsse, spdSecSqrt32);
    run_speed_f32(string("mats sqrt sse"), 15, "sqrt", tests, minVal, maxVal, sqrt32_sse, spdSecSqrt32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Log\n");
    f32 spdSecLog32 = run_speed_f32(string("stdlib log"), 15, "log", tests, minVal, maxVal, logf, 0.0f);
    run_speed_f32(string("mats log"), 15, "log", tests, minVal, maxVal, log32_nonsse, spdSecLog32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Log2\n");
    f32 spdSecLog232 = run_speed_f32(string("stdlib log2"), 15, "log2", tests, minVal, maxVal, log2f, 0.0f);
    run_speed_f32(string("mats log2"), 15, "log2", tests, minVal, maxVal, log2_32_nonsse, spdSecLog232);
    fprintf(stdout, "\n");

    minVal = -12.0f;
    maxVal = 12.0f;
    fprintf(stdout, "Exp\n");
    f32 spdSecExp32 = run_speed_f32(string("stdlib exp"), 15, "exp", tests, minVal, maxVal, expf, 0.0f);
    run_speed_f32(string("mats exp"), 15, "exp", tests, minVal, maxVal, exp32_nonsse, spdSecExp32);
    run_speed_f32(string("mats exp sse"), 15, "exp", tests, minVal, maxVal, exp32_sse, spdSecExp32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Exp2\n");
    f32 spdSecExp232 = run_speed_f32(string("stdlib exp2"), 15, "exp2", tests, minVal, maxVal, exp2f, 0.0f);
    run_speed_f32(string("mats exp2"), 15, "exp2", tests, minVal, maxVal, exp2_32_nonsse, spdSecExp232);
    fprintf(stdout, "\n");

    tests = 5000000;
#if 1
    fprintf(stdout, "Exp wide\n");
    f32 spdSecExp4x = run_speed_f32_4x(string("stdlib exp 4x"), 15, "exp_4x", tests, minVal, maxVal, expf_4x, 0.0f);
    run_speed_f32_4x(string("mats exp 4x"), 15, "exp_4x", tests, minVal, maxVal, exp32_4x, spdSecExp4x);
    fprintf(stdout, "\n");
#endif

    return 0;
}
