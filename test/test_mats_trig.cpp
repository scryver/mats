#include <errno.h>
#include <time.h>
#include <math.h>
#include <x86intrin.h>
#include <xmmintrin.h>

#define NO_INTRINSICS 1
#define STR_FMT(x)   safe_truncate_to_s32(x.size), (char *)x.data
#include "../libberdip/src/common.h"
#include "../libberdip/src/multilane.h"
#include "../libberdip/src/trigonometry_v2.h"

#define MATS_USE_SSE4 1
#define MATS_USE_SSE2 1
#define IEEE_754_2008_SNAN 1
#include "../mats/mats_common.h"
#include "../mats/mats_constants.h"
#include "../mats/mats_defines.h"
#include "../mats/mats_elem.h"
#include "../mats/mats_trig.h"

#include "../src/arm_sincos.cpp"

#include "test_common.cpp"

internal f64
sincos_match64(f64 val)
{
    return 1.0;
}

internal f32
sincosf_match(f32 val)
{
    return square(sinf(val)) + square(cosf(val));
}

internal f32
sincos32_match(f32 val)
{
    return square(sin32(val)) + square(cos32(val));
}

internal f32
sincos32_prec_match(f32 val)
{
    return square(sin32_prec(val)) + square(cos32_prec(val));
}

internal f32
arm_sincosf_match(f32 val)
{
    return square(arm_sinf(val)) + square(arm_cosf(val));
}

internal f32
sincos_pi_match(f32 val)
{
    return square(sin_pi(val)) + square(cos_pi(val));
}

s32 main(s32 argc, char **argv)
{
    u32 tests = 50000000;
    f32 minVal = -F32_TAU; // -10.0f;
    f32 maxVal = F32_TAU; // 10.0f;

    fprintf(stdout, "\nCorrectness\n\n");

    fprintf(stdout, "Sin\n");
    f32 stdSecSin32 = run_comp_f32_x(string("stdlib sin"), 15, "sin", tests, minVal, maxVal, sin, sinf, 0.0f);
    run_comp_f32_x(string("mats sin"), 15, "sin", tests, minVal, maxVal, sin, sin32, stdSecSin32);
    run_comp_f32_x(string("mats sin prec"), 15, "sin", tests, minVal, maxVal, sin, sin32_prec, stdSecSin32);
    run_comp_f32_x(string("arm sin"), 15, "sin", tests, minVal, maxVal, sin, arm_sinf, stdSecSin32);
    run_comp_f32_x(string("dip sin"), 15, "sin", tests, minVal, maxVal, sin, sin_pi, stdSecSin32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Cos\n");
    f32 stdSecCos32 = run_comp_f32_x(string("stdlib cos"), 15, "cos", tests, minVal, maxVal, cos, cosf, 0.0f);
    run_comp_f32_x(string("mats cos"), 15, "cos", tests, minVal, maxVal, cos, cos32, stdSecCos32);
    run_comp_f32_x(string("mats cos prec"), 15, "cos", tests, minVal, maxVal, cos, cos32_prec, stdSecCos32);
    run_comp_f32_x(string("arm cos"), 15, "cos", tests, minVal, maxVal, cos, arm_cosf, stdSecCos32);
    run_comp_f32_x(string("dip cos"), 15, "cos", tests, minVal, maxVal, cos, cos_pi, stdSecCos32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Sin Cos matching\n");
    f32 stdSCMatchSin32 = run_comp_f32_x(string("stdlib match"), 15, "match", tests, minVal, maxVal, sincos_match64, sincosf_match, 0.0f);
    run_comp_f32_x(string("mats match"), 15, "match", tests, minVal, maxVal, sincos_match64, sincos32_match, stdSCMatchSin32);
    run_comp_f32_x(string("mats match prc"), 15, "match", tests, minVal, maxVal, sincos_match64, sincos32_prec_match, stdSCMatchSin32);
    run_comp_f32_x(string("arm match"), 15, "match", tests, minVal, maxVal, sincos_match64, arm_sincosf_match, stdSCMatchSin32);
    run_comp_f32_x(string("dip match"), 15, "match", tests, minVal, maxVal, sincos_match64, sincos_pi_match, stdSCMatchSin32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Tan\n");
    f32 stdSecTan32 = run_comp_f32_x(string("stdlib tan"), 15, "tan", tests, minVal, maxVal, tan, tanf, 0.0f);
    run_comp_f32_x(string("mats tan"), 15, "tan", tests, minVal, maxVal, tan, tan32, stdSecTan32);
    fprintf(stdout, "\n");

    minVal = -100.0f;
    maxVal = 100.0f;

    fprintf(stdout, "ATan\n");
    f32 stdSecATan32 = run_comp_f32_x(string("stdlib atan"), 15, "atan", tests, minVal, maxVal, atan, atanf, 0.0f);
    run_comp_f32_x(string("mats atan"), 15, "atan", tests, minVal, maxVal, atan, atan32, stdSecATan32);
    fprintf(stdout, "\n");

    minVal = -1.0f;
    maxVal = 1.0f;

    fprintf(stdout, "ASin\n");
    f32 stdSecASin32 = run_comp_f32_x(string("stdlib asin"), 15, "asin", tests, minVal, maxVal, asin, asinf, 0.0f);
    run_comp_f32_x(string("mats asin"), 15, "asin", tests, minVal, maxVal, asin, asin32, stdSecASin32);
    fprintf(stdout, "\n");

    fprintf(stdout, "ACos\n");
    f32 stdSecACos32 = run_comp_f32_x(string("stdlib acos"), 15, "acos", tests, minVal, maxVal, acos, acosf, 0.0f);
    run_comp_f32_x(string("mats acos"), 15, "acos", tests, minVal, maxVal, acos, acos32, stdSecACos32);
    fprintf(stdout, "\n");

#if 0
    fprintf(stdout, "SinCos\n");
    f32 stdSecSinCos32 = run_comp_f32(string("stdlib sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, sincosf, 0.0f);
    run_comp_f32(string("mats sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, sincos32, stdSecSinCos32);
    fprintf(stdout, "\n");
#endif

    fprintf(stdout, "Speed\n\n");
    minVal = -F32_TAU; // -10.0f;
    maxVal = F32_TAU; // 10.0f;

    fprintf(stdout, "Sin\n");
    f32 spdSecSin32 = run_speed_f32(string("stdlib sin"), 15, "sin", tests, minVal, maxVal, sinf, 0.0f);
    run_speed_f32(string("mats sin"), 15, "sin", tests, minVal, maxVal, sin32, spdSecSin32);
    run_speed_f32(string("mats sin prec"), 15, "sin", tests, minVal, maxVal, sin32_prec, spdSecSin32);
    run_speed_f32(string("arm sin"), 15, "sin", tests, minVal, maxVal, arm_sinf, spdSecSin32);
    run_speed_f32(string("dip sin"), 15, "sin", tests, minVal, maxVal, sin_pi, spdSecSin32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Cos\n");
    f32 spdSecCos32 = run_speed_f32(string("stdlib cos"), 15, "cos", tests, minVal, maxVal, cosf, 0.0f);
    run_speed_f32(string("mats cos"), 15, "cos", tests, minVal, maxVal, cos32, spdSecCos32);
    run_speed_f32(string("mats cos prec"), 15, "cos", tests, minVal, maxVal, cos32_prec, spdSecCos32);
    run_speed_f32(string("arm cos"), 15, "cos", tests, minVal, maxVal, arm_cosf, spdSecCos32);
    run_speed_f32(string("dip cos"), 15, "cos", tests, minVal, maxVal, cos_pi, spdSecCos32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Tan\n");
    f32 spdSecTan32 = run_speed_f32(string("stdlib tan"), 15, "tan", tests, minVal, maxVal, tanf, 0.0f);
    run_speed_f32(string("mats tan"), 15, "tan", tests, minVal, maxVal, tan32, spdSecTan32);
    fprintf(stdout, "\n");

    minVal = -100.0f;
    maxVal = 100.0f;

    fprintf(stdout, "ATan\n");
    f32 spdSecATan32 = run_speed_f32(string("stdlib atan"), 15, "atan", tests, minVal, maxVal, atanf, 0.0f);
    run_speed_f32(string("mats atan"), 15, "atan", tests, minVal, maxVal, atan32, spdSecATan32);
    fprintf(stdout, "\n");

    minVal = -1.0f;
    maxVal = 1.0f;

    fprintf(stdout, "ASin\n");
    f32 spdSecASin32 = run_speed_f32(string("stdlib asin"), 15, "asin", tests, minVal, maxVal, asinf, 0.0f);
    run_speed_f32(string("mats asin"), 15, "asin", tests, minVal, maxVal, asin32, spdSecASin32);
    fprintf(stdout, "\n");

    fprintf(stdout, "ACos\n");
    f32 spdSecACos32 = run_speed_f32(string("stdlib acos"), 15, "acos", tests, minVal, maxVal, acosf, 0.0f);
    run_speed_f32(string("mats acos"), 15, "acos", tests, minVal, maxVal, acos32, spdSecACos32);
    fprintf(stdout, "\n");

#if 0
    fprintf(stdout, "SinCos\n");
    f32 spdSecSinCos32 = run_speed_f32(string("stdlib sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, 0.0f);
    run_speed_f32(string("mats sincos"), 15, "sincos", tests, minVal, maxVal, sincos32, spdSecSinCos32);
    fprintf(stdout, "\n");
#endif

    return 0;
}
