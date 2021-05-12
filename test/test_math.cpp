#include <errno.h>
#include <time.h>
#include <math.h>
#include <x86intrin.h>
#include <xmmintrin.h>

#define STR_FMT(x)   safe_truncate_to_s32(x.size), (char *)x.data
#include "../libberdip/src/common.h"
#include "../libberdip/src/multilane.h"

#include "../src/sun93.h"

#if 0
internal f32
absolute(f32 x)
{
    U32F32 xf = u32f32(x);
    xf.u &= ~F32_SIGN_MASK;
    return xf.f;
}
#endif

#include "../mats/mats_defines.h"
#include "../mats/mats_common.h"
#include "../mats/mats_constants.h"

//#include "../src/constants.cpp"
#include "../src/sun93_sincos.cpp"
#include "../src/arm_sincos.cpp"
#include "../src/arm_exp_sincos.cpp"
#include "../src/sincos_4x.cpp"
#include "../src/sun93_elem.cpp"

#include "../mats/mats_elem.h"
#include "../mats/mats_trig.h"

#include "test_common.cpp"

internal void
print_info(String name, u32 maxNameSize, u32 calls, u32 tests, f32 seconds, f64 sum, f64 errSum, f32 secondsRatio)
{
    f64 testCalls = (f64)((u64)calls * tests);
    persist char *spaces = "                                                                ";
    fprintf(stdout, "%.*s%.*s: %f, %f sec, %10.0f sin/sec, avg %f usec, %f, %f, err %e\n",
            STR_FMT(name), (u32)(maxNameSize - name.size), spaces, secondsRatio,
            seconds, testCalls / seconds, seconds / (testCalls / 1000000.0), sum, sum * (1.0 / testCalls), errSum * (1.0 / testCalls));
}

#define TRIG_FUNC_F32(name)    f32 name(f32 angle)
typedef TRIG_FUNC_F32(TrigFuncF32);
#define TRIG_FUNC_F32_4X(name) f32_4x name(f32_4x angles)
typedef TRIG_FUNC_F32_4X(TrigFuncF32_4x);
#define SINCOS_F32_4x(name)    SinCos4x name(f32_4x angles)
typedef SINCOS_F32_4x(SinCosF32_4x);

f32
run_test_f32(String name, u32 maxNameSize, u32 tests, f32 scaleFactor,
             TrigFuncF32 *sin_func, TrigFuncF32 *cos_func,
             f32 secondsBase)
{
    f64 sum = 0.0;
    f64 errSum = 0.0;
    f32 oneOverTests = 1.0f / (f32)tests;
    f32 scale = scaleFactor * oneOverTests;
    f32 halfFactor = scaleFactor * 0.5f;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f32 ratio0 = (f32)index * scale - halfFactor;
        f32 ratio1 = (f32)index * 0.5f * scale - halfFactor;
        f32 ratio2 = (f32)index * 0.4f * scale - halfFactor;
        f32 ratio3 = (f32)index * 0.3f * scale - halfFactor;
        f32 sine0 = sin_func(ratio0);
        f32 sine1 = sin_func(ratio1);
        f32 sine2 = sin_func(ratio2);
        f32 sine3 = sin_func(ratio3);
        f32 cosine0 = cos_func(ratio0);
        f32 cosine1 = cos_func(ratio1);
        f32 cosine2 = cos_func(ratio2);
        f32 cosine3 = cos_func(ratio3);
        sum += sine0;
        sum += cosine0;
        errSum += 1.0f - (square(sine0) + square(cosine0));
        errSum += 1.0f - (square(sine1) + square(cosine1));
        errSum += 1.0f - (square(sine2) + square(cosine2));
        errSum += 1.0f - (square(sine3) + square(cosine3));
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 calls = 4;

    if (secondsBase == 0.0f) {
        print_info(name, maxNameSize, calls, tests, seconds, sum, errSum, 1.0f);
    } else {
        print_info(name, maxNameSize, calls, tests, seconds, sum, errSum, seconds / secondsBase);
    }
    return seconds;
}

f32
run_test_f32_4x(String name, u32 maxNameSize, u32 tests, f32 scaleFactor,
                TrigFuncF32_4x *sin_func, TrigFuncF32_4x *cos_func,
                f32 secondsBase)
{
    f64 sum = 0.0;
    f64 errSum = 0.0;
    f32 oneOverTests = 1.0f / (f32)tests;
    f32_4x scale;
    scale.m = _mm_set1_ps(scaleFactor * oneOverTests);
    f32_4x halfFactor;
    halfFactor.m = _mm_set1_ps(scaleFactor * 0.5f);
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f32_4x ratio1;
        ratio1.m = _mm_setr_ps((f32)index, (f32)index * 0.5f, (f32)index * 0.4f, (f32)index * 0.3f);
        ratio1.m = _mm_sub_ps(_mm_mul_ps(ratio1.m, scale.m), halfFactor.m);

        f32_4x sin1 = sin_func(ratio1);
        f32_4x cos1 = cos_func(ratio1);
        sum += sin1.e[0];
        sum += cos1.e[0];

        f32_4x err;
        err.m = _mm_sub_ps(_mm_set1_ps(1.0f),
                           _mm_add_ps(_mm_mul_ps(sin1.m, sin1.m),
                                      _mm_mul_ps(cos1.m, cos1.m)));

        errSum += err.e[0];
        errSum += err.e[1];
        errSum += err.e[2];
        errSum += err.e[3];
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 calls = 4;

    if (secondsBase == 0.0f) {
        print_info(name, maxNameSize, calls, tests, seconds, sum, errSum, 1.0f);
    } else {
        print_info(name, maxNameSize, calls, tests, seconds, sum, errSum, seconds / secondsBase);
    }
    return seconds;
}

internal f32
run_test_sincos_4x(String name, u32 maxNameSize, u32 tests, f32 scaleFactor,
                   SinCosF32_4x *sincos_func, f32 secondsBase)
{
    f64 sum = 0.0;
    f64 errSum = 0.0;
    f32 oneOverTests = 1.0f / (f32)tests;
    f32_4x scale;
    scale.m = _mm_set1_ps(scaleFactor * oneOverTests);
    f32_4x halfFactor;
    halfFactor.m = _mm_set1_ps(scaleFactor * 0.5f);
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f32_4x ratio1;
        ratio1.m = _mm_setr_ps((f32)index, (f32)index * 0.5f, (f32)index * 0.4f, (f32)index * 0.3f);
        ratio1.m = _mm_sub_ps(_mm_mul_ps(ratio1.m, scale.m), halfFactor.m);

        SinCos4x sinCos1 = sincos_func(ratio1);

        sum += sinCos1.sin.e[0];
        sum += sinCos1.cos.e[0];

        f32_4x err;
        err.m = _mm_sub_ps(_mm_set1_ps(1.0f),
                           _mm_add_ps(_mm_mul_ps(sinCos1.sin.m, sinCos1.sin.m),
                                      _mm_mul_ps(sinCos1.cos.m, sinCos1.cos.m)));

        errSum += err.e[0];
        errSum += err.e[1];
        errSum += err.e[2];
        errSum += err.e[3];
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 calls = 4;

    if (secondsBase == 0.0f) {
        print_info(name, maxNameSize, calls, tests, seconds, sum, errSum, 1.0f);
    } else {
        print_info(name, maxNameSize, calls, tests, seconds, sum, errSum, seconds / secondsBase);
    }
    return seconds;
}

#include "test_sincos.cpp"

s32 main(s32 argc, char **argv)
{
    fprintf(stdout, "Sin/cos test\n");

    fprintf(stdout, "0: %g\n", sun93_cosf(0.0f));
    fprintf(stdout, "1: %g\n", sun93_cosf(1.0f));
    fprintf(stdout, "2: %g\n", sun93_cosf(2.0f));
    fprintf(stdout, "PI: %g\n", sun93_cosf(F32_PI));
    fprintf(stdout, "PI/2: %g\n", sun93_cosf(F32_PI*0.5f));

    fprintf(stdout, "0: %g\n", sun93_sinf(0.0f));
    fprintf(stdout, "1: %g\n", sun93_sinf(1.0f));
    fprintf(stdout, "2: %g\n", sun93_sinf(2.0f));
    fprintf(stdout, "PI: %g\n", sun93_sinf(F32_PI));
    fprintf(stdout, "PI/2: %g\n", sun93_sinf(F32_PI*0.5f));

#define SUM(x) (sun93_sinf(x)*sun93_sinf(x) + sun93_cosf(x)*sun93_cosf(x))
    fprintf(stdout, "0: %g\n", SUM(0.0f));
    fprintf(stdout, "1: %g\n", SUM(1.0f));
    fprintf(stdout, "2: %g\n", SUM(2.0f));
    fprintf(stdout, "PI: %g\n", SUM(F32_PI));
    fprintf(stdout, "PI/2: %g\n", SUM(F32_PI*0.5f));
#undef SUM

    fprintf(stdout, "ARM\n");
    fprintf(stdout, "0: %g\n", arm_cosf(0.0f));
    fprintf(stdout, "1: %g\n", arm_cosf(1.0f));
    fprintf(stdout, "2: %g\n", arm_cosf(2.0f));
    fprintf(stdout, "PI: %g\n", arm_cosf(F32_PI));
    fprintf(stdout, "PI/2: %g\n", arm_cosf(F32_PI*0.5f));

    fprintf(stdout, "0: %g\n", arm_sinf(0.0f));
    fprintf(stdout, "1: %g\n", arm_sinf(1.0f));
    fprintf(stdout, "2: %g\n", arm_sinf(2.0f));
    fprintf(stdout, "PI: %g\n", arm_sinf(F32_PI));
    fprintf(stdout, "PI/2: %g\n", arm_sinf(F32_PI*0.5f));

#define SUM(x) (arm_sinf(x)*arm_sinf(x) + arm_cosf(x)*arm_cosf(x))
    fprintf(stdout, "0: %g\n", SUM(0.0f));
    fprintf(stdout, "1: %g\n", SUM(1.0f));
    fprintf(stdout, "2: %g\n", SUM(2.0f));
    fprintf(stdout, "PI: %g\n", SUM(F32_PI));
    fprintf(stdout, "PI/2: %g\n", SUM(F32_PI*0.5f));
#undef SUM

    fprintf(stdout, "ARM Expanded\n");
    fprintf(stdout, "0: %g\n", arm_exp_cosf(0.0f));
    fprintf(stdout, "1: %g\n", arm_exp_cosf(1.0f));
    fprintf(stdout, "2: %g\n", arm_exp_cosf(2.0f));
    fprintf(stdout, "PI: %g\n", arm_exp_cosf(F32_PI));
    fprintf(stdout, "PI/2: %g\n", arm_exp_cosf(F32_PI*0.5f));

    fprintf(stdout, "0: %g\n", arm_exp_sinf(0.0f));
    fprintf(stdout, "1: %g\n", arm_exp_sinf(1.0f));
    fprintf(stdout, "2: %g\n", arm_exp_sinf(2.0f));
    fprintf(stdout, "PI: %g\n", arm_exp_sinf(F32_PI));
    fprintf(stdout, "PI/2: %g\n", arm_exp_sinf(F32_PI*0.5f));

#define SUM(x) (arm_exp_sinf(x)*arm_exp_sinf(x) + arm_exp_cosf(x)*arm_exp_cosf(x))
    fprintf(stdout, "0: %g\n", SUM(0.0f));
    fprintf(stdout, "1: %g\n", SUM(1.0f));
    fprintf(stdout, "2: %g\n", SUM(2.0f));
    fprintf(stdout, "PI: %g\n", SUM(F32_PI));
    fprintf(stdout, "PI/2: %g\n", SUM(F32_PI*0.5f));
#undef SUM

    f32_4x zeroes;
    f32_4x ones;
    f32_4x twos;
    f32_4x pies;
    f32_4x halfPies;
    zeroes.m = _mm_setzero_ps();
    ones.m = _mm_set1_ps(1.0f);
    twos.m = _mm_set1_ps(2.0f);
    pies.m = _mm_set1_ps(F32_PI);
    halfPies.m = _mm_set1_ps(0.5f * F32_PI);
    fprintf(stdout, "ARMish 4x\n");
    fprintf(stdout, "0: %g\n", arm_cosf_4x(zeroes).e[0]);
    fprintf(stdout, "1: %g\n", arm_cosf_4x(ones).e[0]);
    fprintf(stdout, "2: %g\n", arm_cosf_4x(twos).e[0]);
    fprintf(stdout, "PI: %g\n", arm_cosf_4x(pies).e[0]);
    fprintf(stdout, "PI/2: %g\n", arm_cosf_4x(halfPies).e[0]);

    fprintf(stdout, "0: %g\n", arm_sinf_4x(zeroes).e[0]);
    fprintf(stdout, "1: %g\n", arm_sinf_4x(ones).e[0]);
    fprintf(stdout, "2: %g\n", arm_sinf_4x(twos).e[0]);
    fprintf(stdout, "PI: %g\n", arm_sinf_4x(pies).e[0]);
    fprintf(stdout, "PI/2: %g\n", arm_sinf_4x(halfPies).e[0]);

#define SUM(x) (arm_sinf_4x(x).e[0]*arm_sinf_4x(x).e[0] + arm_cosf_4x(x).e[0]*arm_cosf_4x(x).e[0])
    fprintf(stdout, "0: %g\n", SUM(zeroes));
    fprintf(stdout, "1: %g\n", SUM(ones));
    fprintf(stdout, "2: %g\n", SUM(twos));
    fprintf(stdout, "PI: %g\n", SUM(pies));
    fprintf(stdout, "PI/2: %g\n", SUM(halfPies));
#undef SUM

    fprintf(stdout, "Sqrt 100, 10, 1, 0.1: %g %g %g %g\n", sqrtf(100.0f), sqrtf(10.0f), sqrtf(1.0f), sqrtf(0.1f));
    fprintf(stdout, "Sqrt 100, 10, 1, 0.1: %g %g %g %g\n", sun93_sqrtf(100.0f), sun93_sqrtf(10.0f), sun93_sqrtf(1.0f), sun93_sqrtf(0.1f));

    test_cosf();
    test_sinf();

    u32 tests = 20000000;

#if 0
    f32 scaler = 2.0f*F32_TAU;
    f32 secondsStd = run_test_f32(string("stdlib f32"), 14, tests, scaler, sinf, cosf, 0.0f);
    run_test_f32(string("sun f32"), 14, tests, scaler, sun93_sinf, sun93_cosf, secondsStd);
    run_test_f32(string("arm f32"), 14, tests, scaler, arm_sinf, arm_cosf, secondsStd);
    run_test_f32(string("arm exp f32"), 14, tests, scaler, arm_exp_sinf, arm_exp_cosf, secondsStd);
    run_test_f32_4x(string("arm f32 4x"), 14, tests, scaler, arm_sinf_4x, arm_cosf_4x, secondsStd);
    run_test_sincos_4x(string("sincos 4x"), 14, tests, scaler, arm_sincosf_4x, secondsStd);
#endif

    f32 stdSecSin = run_comp_f32(string("stdlib sinf"), 14, "sin", tests, -F32_TAU, F32_TAU, sinf, sinf, 0.0f);
    run_comp_f32(string("trig sin"), 14, "sin", tests, -F32_TAU, F32_TAU, sinf, sin32, stdSecSin);
    run_comp_f32(string("sun sinf"), 14, "sin", tests, -F32_TAU, F32_TAU, sinf, sun93_sinf, stdSecSin);
    run_comp_f32(string("arm sinf"), 14, "sin", tests, -F32_TAU, F32_TAU, sinf, arm_sinf, stdSecSin);
    run_comp_f32(string("aexp sinf"), 14, "sin", tests, -F32_TAU, F32_TAU, sinf, arm_exp_sinf, stdSecSin);

    stdSecSin = run_comp_f32(string("stdlib sinf"), 14, "sin", tests, -0.1f, 0.1f, sinf, sinf, 0.0f);
    run_comp_f32(string("trig sin"), 14, "sin", tests, -0.1f, 0.1f, sinf, sin32, stdSecSin);
    run_comp_f32(string("sun sinf"), 14, "sin", tests, -0.1f, 0.1f, sinf, sun93_sinf, stdSecSin);
    run_comp_f32(string("arm sinf"), 14, "sin", tests, -0.1f, 0.1f, sinf, arm_sinf, stdSecSin);
    run_comp_f32(string("aexp sinf"), 14, "sin", tests, -0.1f, 0.1f, sinf, arm_exp_sinf, stdSecSin);

    f32 stdSecSqrt = run_comp_f32(string("stdlib sqrtf"), 14, "sqrt", tests, 1.0e-5f, 100.0f, sqrtf, sqrtf, 0.0f);
    run_comp_f32(string("sun sqrtf"), 14, "sqrt", tests, 1.0e-5f, 100.0f, sqrtf, sun93_sqrtf, stdSecSqrt);
    run_comp_f32(string("sse sqrtf"), 14, "sqrt", tests, 1.0e-5f, 100.0f, sqrtf, sqrt32, stdSecSqrt);

    f32 stdSecAtan = run_comp_f32(string("stdlib atanf"), 14, "atan", tests, -1.0f, 1.0f, atanf, atanf, 0.0f);
    run_comp_f32(string("sun atanf"), 14, "atan", tests, -1.0f, 1.0f, atanf, sun93_atanf, stdSecAtan);

    return 0;
}
