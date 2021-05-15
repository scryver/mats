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
#include "../mats/mats_elem_ext.h"
#include "../mats/mats_elem4x.h"
#include "../mats/mats_trig.h"
#include "../mats/mats_trig4x.h"

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
arm_sincosf_match(f32 val)
{
    return square(arm_sinf(val)) + square(arm_cosf(val));
}

internal f32
sincos_pi_match(f32 val)
{
    return square(sin_pi(val)) + square(cos_pi(val));
}

internal f32_4x
sinf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = sinf(x.e[0]);
    result.e[1] = sinf(x.e[1]);
    result.e[2] = sinf(x.e[2]);
    result.e[3] = sinf(x.e[3]);
    return result;
}

internal f32_4x
cosf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = cosf(x.e[0]);
    result.e[1] = cosf(x.e[1]);
    result.e[2] = cosf(x.e[2]);
    result.e[3] = cosf(x.e[3]);
    return result;
}

internal f32_4x
tanf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = tanf(x.e[0]);
    result.e[1] = tanf(x.e[1]);
    result.e[2] = tanf(x.e[2]);
    result.e[3] = tanf(x.e[3]);
    return result;
}

internal f32_4x
acosf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = acosf(x.e[0]);
    result.e[1] = acosf(x.e[1]);
    result.e[2] = acosf(x.e[2]);
    result.e[3] = acosf(x.e[3]);
    return result;
}

internal f32_4x
acos32_temp_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = acos32(x.e[0]);
    result.e[1] = acos32(x.e[1]);
    result.e[2] = acos32(x.e[2]);
    result.e[3] = acos32(x.e[3]);
    return result;
}

internal f32_4x
asinf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = asinf(x.e[0]);
    result.e[1] = asinf(x.e[1]);
    result.e[2] = asinf(x.e[2]);
    result.e[3] = asinf(x.e[3]);
    return result;
}

internal f32_4x
asin32_temp_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = asin32(x.e[0]);
    result.e[1] = asin32(x.e[1]);
    result.e[2] = asin32(x.e[2]);
    result.e[3] = asin32(x.e[3]);
    return result;
}

internal f32_4x
atanf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = atanf(x.e[0]);
    result.e[1] = atanf(x.e[1]);
    result.e[2] = atanf(x.e[2]);
    result.e[3] = atanf(x.e[3]);
    return result;
}

internal f32_4x
atan32_temp_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = atan32(x.e[0]);
    result.e[1] = atan32(x.e[1]);
    result.e[2] = atan32(x.e[2]);
    result.e[3] = atan32(x.e[3]);
    return result;
}

enum DoTestFlag
{
    DoTest_Cos        = 0x00000001,
    DoTest_Sin        = 0x00000002,
    DoTest_Tan        = 0x00000004,
    DoTest_ACos       = 0x00000008,
    DoTest_ASin       = 0x00000010,
    DoTest_ATan       = 0x00000020,
    DoTest_ATan2      = 0x00000040,
    DoTest_CosH       = 0x00000080,
    DoTest_SinH       = 0x00000100,
    DoTest_TanH       = 0x00000200,
    DoTest_ACosH      = 0x00000400,
    DoTest_ASinH      = 0x00000800,
    DoTest_ATanH      = 0x00001000,
    DoTest_FuncMask   = 0x0000FFFF,

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
            } else if (strings_are_equal("cos", argv[index])) {
                tests |= DoTest_Cos;
            } else if (strings_are_equal("sin", argv[index])) {
                tests |= DoTest_Sin;
            } else if (strings_are_equal("tan", argv[index])) {
                tests |= DoTest_Tan;
            } else if (strings_are_equal("acos", argv[index])) {
                tests |= DoTest_ACos;
            } else if (strings_are_equal("asin", argv[index])) {
                tests |= DoTest_ASin;
            } else if (strings_are_equal("atan", argv[index])) {
                tests |= DoTest_ATan;
            } else if (strings_are_equal("atan2", argv[index])) {
                tests |= DoTest_ATan2;
            } else if (strings_are_equal("cosh", argv[index])) {
                tests |= DoTest_CosH;
            } else if (strings_are_equal("sinh", argv[index])) {
                tests |= DoTest_SinH;
            } else if (strings_are_equal("tanh", argv[index])) {
                tests |= DoTest_TanH;
            } else if (strings_are_equal("acosh", argv[index])) {
                tests |= DoTest_ACosH;
            } else if (strings_are_equal("asinh", argv[index])) {
                tests |= DoTest_ASinH;
            } else if (strings_are_equal("atanh", argv[index])) {
                tests |= DoTest_ATanH;
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
    f32 minVal = -100.0f;
    f32 maxVal = 100.0f;

    u32 testsA = 10000;
    u32 testsB = 5000;
    f32 minValA = -100.0f;
    f32 maxValA = 100.0f;
    f32 minValB = -100.0f;
    f32 maxValB = 100.0f;

    if (doTests & DoTest_Comp)
    {
        fprintf(stdout, "\nCorrectness\n\n");

        if (doTests & DoTest_Cos)
        {
            fprintf(stdout, "Cos\n");
            f32 stdSecCos32 = run_comp_f32_x(string("stdlib cos"), 15, "cos", tests, minVal, maxVal, cos, cosf, 0.0f);
            run_comp_f32_x(string("mats cos"), 15, "cos", tests, minVal, maxVal, cos, cos32, stdSecCos32);
            run_comp_f32_x(string("arm cos"), 15, "cos", tests, minVal, maxVal, cos, arm_cosf, stdSecCos32);
            run_comp_f32_x(string("dip cos"), 15, "cos", tests, minVal, maxVal, cos, cos_pi, stdSecCos32);
            run_comp_f32_4x_x(string("mats cos 4x"), 15, "cos_4x", tests, minVal, maxVal, cos, cos32_4x, stdSecCos32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Sin)
        {
            fprintf(stdout, "Sin\n");
            f32 stdSecSin32 = run_comp_f32_x(string("stdlib sin"), 15, "sin", tests, minVal, maxVal, sin, sinf, 0.0f);
            run_comp_f32_x(string("mats sin"), 15, "sin", tests, minVal, maxVal, sin, sin32, stdSecSin32);
            run_comp_f32_x(string("arm sin"), 15, "sin", tests, minVal, maxVal, sin, arm_sinf, stdSecSin32);
            run_comp_f32_x(string("dip sin"), 15, "sin", tests, minVal, maxVal, sin, sin_pi, stdSecSin32);
            run_comp_f32_4x_x(string("mats sin 4x"), 15, "sin_4x", tests, minVal, maxVal, sin, sin32_4x, stdSecSin32);
            fprintf(stdout, "\n");
        }

        if ((doTests & (DoTest_Sin | DoTest_Cos)) == (DoTest_Sin | DoTest_Cos))
        {
            fprintf(stdout, "Sin Cos matching\n");
            f32 stdSCMatchSin32 = run_comp_f32_x(string("stdlib match"), 15, "match", tests, minVal, maxVal, sincos_match64, sincosf_match, 0.0f);
            run_comp_f32_x(string("mats match"), 15, "match", tests, minVal, maxVal, sincos_match64, sincos32_match, stdSCMatchSin32);
            run_comp_f32_x(string("arm match"), 15, "match", tests, minVal, maxVal, sincos_match64, arm_sincosf_match, stdSCMatchSin32);
            run_comp_f32_x(string("dip match"), 15, "match", tests, minVal, maxVal, sincos_match64, sincos_pi_match, stdSCMatchSin32);
            // TODO(michiel): Add 4x
            fprintf(stdout, "\n");

#if 0
            fprintf(stdout, "SinCos\n");
            f32 stdSecSinCos32 = run_comp_f32(string("stdlib sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, sincosf, 0.0f);
            run_comp_f32(string("mats sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, sincos32, stdSecSinCos32);
            fprintf(stdout, "\n");
#endif
        }

        if (doTests & DoTest_Tan)
        {
            fprintf(stdout, "Tan\n");
            f32 stdSecTan32 = run_comp_f32_x(string("stdlib tan"), 15, "tan", tests, minVal, maxVal, tan, tanf, 0.0f);
            run_comp_f32_x(string("mats tan"), 15, "tan", tests, minVal, maxVal, tan, tan32, stdSecTan32);
            run_comp_f32_4x_x(string("mats tan 4x"), 15, "tan_4x", tests, minVal, maxVal, tan, tan32_4x, stdSecTan32);
            fprintf(stdout, "\n");
        }

        minVal = -1.2f;
        maxVal =  1.2f;

        if (doTests & DoTest_ACos)
        {
            fprintf(stdout, "ACos\n");
            f32 stdSecACos32 = run_comp_f32_x(string("stdlib acos"), 15, "acos", tests, minVal, maxVal, acos, acosf, 0.0f);
            run_comp_f32_x(string("mats acos"), 15, "acos", tests, minVal, maxVal, acos, acos32, stdSecACos32);
            run_comp_f32_4x_x(string("mats acos 4x"), 15, "acos_4x", tests, minVal, maxVal, acos, acos32_4x, stdSecACos32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_ASin)
        {
            fprintf(stdout, "ASin\n");
            f32 stdSecASin32 = run_comp_f32_x(string("stdlib asin"), 15, "asin", tests, minVal, maxVal, asin, asinf, 0.0f);
            run_comp_f32_x(string("mats asin"), 15, "asin", tests, minVal, maxVal, asin, asin32, stdSecASin32);
            run_comp_f32_4x_x(string("mats asin 4x"), 15, "asin_4x", tests, minVal, maxVal, asin, asin32_4x, stdSecASin32);
            fprintf(stdout, "\n");
        }

        minVal = -100.0f;
        maxVal = 100.0f;

        if (doTests & DoTest_ATan)
        {
            fprintf(stdout, "ATan\n");
            f32 stdSecATan32 = run_comp_f32_x(string("stdlib atan"), 15, "atan", tests, minVal, maxVal, atan, atanf, 0.0f);
            run_comp_f32_x(string("mats atan"), 15, "atan", tests, minVal, maxVal, atan, atan32, stdSecATan32);
            run_comp_f32_4x_x(string("mats atan 4x"), 15, "atan_4x", tests, minVal, maxVal, atan, atan32_4x, stdSecATan32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_ATan2)
        {
            fprintf(stdout, "ATan2\n");
            f32 stdSecATan232 = call_comp_x2(stdlib, atan2, f, 0.0f);
            call_comp_x2(mats, atan2, _32, stdSecATan232);
            fprintf(stdout, "\n");
        }

        minVal = -22.0f;
        maxVal = 22.0f;

        if (doTests & DoTest_CosH)
        {
            fprintf(stdout, "CosH\n");
            f32 stdSec = run_comp_f32_x(string("stdlib cosh"), 15, "cosh", tests, minVal, maxVal, cosh, coshf, 0.0f);
            run_comp_f32_x(string("mats cosh"), 15, "cosh", tests, minVal, maxVal, cosh, cosh32, stdSec);
            //run_comp_f32_4x_x(string("mats cosh 4x"), 15, "cosh_4x", tests, minVal, maxVal, cosh, cosh32_4x, stdSec);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_SinH)
        {
            fprintf(stdout, "SinH\n");
            f32 stdSec = run_comp_f32_x(string("stdlib sinh"), 15, "sinh", tests, minVal, maxVal, sinh, sinhf, 0.0f);
            run_comp_f32_x(string("mats sinh"), 15, "sinh", tests, minVal, maxVal, sinh, sinh32, stdSec);
            //run_comp_f32_4x_x(string("mats sinh 4x"), 15, "sinh_4x", tests, minVal, maxVal, sinh, sinh32_4x, stdSec);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_TanH)
        {
            fprintf(stdout, "TanH\n");
            f32 stdSec = run_comp_f32_x(string("stdlib tanh"), 15, "tanh", tests, minVal, maxVal, tanh, tanhf, 0.0f);
            run_comp_f32_x(string("mats tanh"), 15, "tanh", tests, minVal, maxVal, tanh, tanh32, stdSec);
            //run_comp_f32_4x_x(string("mats tanh 4x"), 15, "tanh_4x", tests, minVal, maxVal, tanh, tanh32_4x, stdSec);
            fprintf(stdout, "\n");
        }

        minVal = 1.0f;
        maxVal = 100.0f;

        if (doTests & DoTest_ACosH)
        {
            fprintf(stdout, "ACosH\n");
            f32 stdSec = run_comp_f32_x(string("stdlib acosh"), 15, "acosh", tests, minVal, maxVal, acosh, acoshf, 0.0f);
            run_comp_f32_x(string("mats acosh"), 15, "acosh", tests, minVal, maxVal, acosh, acosh32, stdSec);
            //run_comp_f32_4x_x(string("mats acosh 4x"), 15, "acosh_4x", tests, minVal, maxVal, acosh, acosh32_4x, stdSec);
            fprintf(stdout, "\n");
        }

        minVal = -100.0f;
        maxVal = 100.0f;

        if (doTests & DoTest_ASinH)
        {
            fprintf(stdout, "ASinH\n");
            f32 stdSec = run_comp_f32_x(string("stdlib asinh"), 15, "asinh", tests, minVal, maxVal, asinh, asinhf, 0.0f);
            run_comp_f32_x(string("mats asinh"), 15, "asinh", tests, minVal, maxVal, asinh, asinh32, stdSec);
            //run_comp_f32_4x_x(string("mats asinh 4x"), 15, "asinh_4x", tests, minVal, maxVal, asinh, asinh32_4x, stdSec);
            fprintf(stdout, "\n");
        }

        minVal = -1.0f;
        maxVal = 1.0f;

        if (doTests & DoTest_ATanH)
        {
            fprintf(stdout, "ATanH\n");
            f32 stdSec = run_comp_f32_x(string("stdlib atanh"), 15, "atanh", tests, minVal, maxVal, atanh, atanhf, 0.0f);
            run_comp_f32_x(string("mats atanh"), 15, "atanh", tests, minVal, maxVal, atanh, atanh32, stdSec);
            //run_comp_f32_4x_x(string("mats atanh 4x"), 15, "atanh_4x", tests, minVal, maxVal, atanh, atanh32_4x, stdSec);
            fprintf(stdout, "\n");
        }
    }

    if (doTests & DoTest_Speed)
    {
        fprintf(stdout, "Speed\n\n");
        minVal = -100.0f;
        maxVal = 100.0f;

        if (doTests & DoTest_Cos)
        {
            fprintf(stdout, "Cos\n");
            f32 spdSecCos32 = run_speed_f32(string("stdlib cos"), 15, "cos", tests, minVal, maxVal, cosf, 0.0f);
            run_speed_f32(string("mats cos"), 15, "cos", tests, minVal, maxVal, cos32, spdSecCos32);
            run_speed_f32(string("arm cos"), 15, "cos", tests, minVal, maxVal, arm_cosf, spdSecCos32);
            run_speed_f32(string("dip cos"), 15, "cos", tests, minVal, maxVal, cos_pi, spdSecCos32);
            run_speed_f32_4x(string("mats cos 4x"), 15, "cos", tests, minVal, maxVal, cos32_4x, spdSecCos32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Sin)
        {
            fprintf(stdout, "Sin\n");
            f32 spdSecSin32 = run_speed_f32(string("stdlib sin"), 15, "sin", tests, minVal, maxVal, sinf, 0.0f);
            run_speed_f32(string("mats sin"), 15, "sin", tests, minVal, maxVal, sin32, spdSecSin32);
            run_speed_f32(string("arm sin"), 15, "sin", tests, minVal, maxVal, arm_sinf, spdSecSin32);
            run_speed_f32(string("dip sin"), 15, "sin", tests, minVal, maxVal, sin_pi, spdSecSin32);
            run_speed_f32_4x(string("mats sin 4x"), 15, "sin", tests, minVal, maxVal, sin32_4x, spdSecSin32);
            fprintf(stdout, "\n");
        }

#if 0
        if ((doTests & (DoTest_Sin | DoTest_Cos)) == (DoTest_Sin | DoTest_Cos))
        {
            fprintf(stdout, "SinCos\n");
            f32 spdSecSinCos32 = run_speed_f32(string("stdlib sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, 0.0f);
            run_speed_f32(string("mats sincos"), 15, "sincos", tests, minVal, maxVal, sincos32, spdSecSinCos32);
            fprintf(stdout, "\n");
        }
#endif

        if (doTests & DoTest_Tan)
        {
            fprintf(stdout, "Tan\n");
            f32 spdSecTan32 = run_speed_f32(string("stdlib tan"), 15, "tan", tests, minVal, maxVal, tanf, 0.0f);
            run_speed_f32(string("mats tan"), 15, "tan", tests, minVal, maxVal, tan32, spdSecTan32);
            run_speed_f32_4x(string("mats tan 4x"), 15, "tan", tests, minVal, maxVal, tan32_4x, spdSecTan32);
            fprintf(stdout, "\n");
        }

        minVal = -1.2f;
        maxVal = 1.2f;

        if (doTests & DoTest_ACos)
        {
            fprintf(stdout, "ACos\n");
            f32 spdSecACos32 = run_speed_f32(string("stdlib acos"), 15, "acos", tests, minVal, maxVal, acosf, 0.0f);
            run_speed_f32(string("mats acos"), 15, "acos", tests, minVal, maxVal, acos32, spdSecACos32);
            run_speed_f32_4x(string("mats acos 4x"), 15, "acos", tests, minVal, maxVal, acos32_4x, spdSecACos32);
            run_speed_f32_4x(string("mats tcos 4x"), 15, "acos", tests, minVal, maxVal, acos32_temp_4x, spdSecACos32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_ASin)
        {
            fprintf(stdout, "ASin\n");
            f32 spdSecASin32 = run_speed_f32(string("stdlib asin"), 15, "asin", tests, minVal, maxVal, asinf, 0.0f);
            run_speed_f32(string("mats asin"), 15, "asin", tests, minVal, maxVal, asin32, spdSecASin32);
            run_speed_f32_4x(string("mats asin 4x"), 15, "asin", tests, minVal, maxVal, asin32_4x, spdSecASin32);
            run_speed_f32_4x(string("mats tsin 4x"), 15, "asin", tests, minVal, maxVal, asin32_temp_4x, spdSecASin32);
            fprintf(stdout, "\n");
        }

        minVal = -100.0f;
        maxVal = 100.0f;

        if (doTests & DoTest_ATan)
        {
            fprintf(stdout, "ATan\n");
            f32 spdSecATan32 = run_speed_f32(string("stdlib atan"), 15, "atan", tests, minVal, maxVal, atanf, 0.0f);
            run_speed_f32(string("mats atan"), 15, "atan", tests, minVal, maxVal, atan32, spdSecATan32);
            run_speed_f32_4x(string("mats atan 4x"), 15, "atan", tests, minVal, maxVal, atan32_4x, spdSecATan32);
            run_speed_f32_4x(string("mats ttan 4x"), 15, "atan", tests, minVal, maxVal, atan32_temp_4x, spdSecATan32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_ATan2)
        {
            fprintf(stdout, "ATan2\n");
            f32 spdSecATan232 = call_spd2(stdlib, atan2, f, 0.0f);
            call_spd2(mats, atan2, _32, spdSecATan232);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_CosH)
        {
            fprintf(stdout, "CosH\n");
            f32 spdSecCosH32 = run_speed_f32(string("stdlib cosh"), 15, "cosh", tests, minVal, maxVal, coshf, 0.0f);
            run_speed_f32(string("mats cosh"), 15, "cosh", tests, minVal, maxVal, cosh32, spdSecCosH32);
            //run_speed_f32_4x(string("mats cosh 4x"), 15, "cosh", tests, minVal, maxVal, cosh32_4x, spdSecCosH32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_SinH)
        {
            fprintf(stdout, "SinH\n");
            f32 spdSecSinH32 = run_speed_f32(string("stdlib sinh"), 15, "sinh", tests, minVal, maxVal, sinhf, 0.0f);
            run_speed_f32(string("mats sinh"), 15, "sinh", tests, minVal, maxVal, sinh32, spdSecSinH32);
            //run_speed_f32_4x(string("mats sinh 4x"), 15, "sinh", tests, minVal, maxVal, sinh32_4x, spdSecSinH32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_TanH)
        {
            fprintf(stdout, "TanH\n");
            f32 spdSecTanH32 = run_speed_f32(string("stdlib tanh"), 15, "tanh", tests, minVal, maxVal, tanhf, 0.0f);
            run_speed_f32(string("mats tanh"), 15, "tanh", tests, minVal, maxVal, tanh32, spdSecTanH32);
            //run_speed_f32_4x(string("mats tanh 4x"), 15, "tanh", tests, minVal, maxVal, tanh32_4x, spdSecTanH32);
            fprintf(stdout, "\n");
        }

        minVal = 1.0f;
        maxVal = 100.0f;

        if (doTests & DoTest_ACosH)
        {
            fprintf(stdout, "ACosH\n");
            f32 stdSec = run_speed_f32(string("stdlib acosh"), 15, "acosh", tests, minVal, maxVal, acoshf, 0.0f);
            run_speed_f32(string("mats acosh"), 15, "acosh", tests, minVal, maxVal, acosh32, stdSec);
            fprintf(stdout, "\n");
        }

        minVal = -100.0f;
        maxVal = 100.0f;

        if (doTests & DoTest_ASinH)
        {
            fprintf(stdout, "ASinH\n");
            f32 stdSec = run_speed_f32(string("stdlib asinh"), 15, "asinh", tests, minVal, maxVal, asinhf, 0.0f);
            run_speed_f32(string("mats asinh"), 15, "asinh", tests, minVal, maxVal, asinh32, stdSec);
            fprintf(stdout, "\n");
        }

        minVal = -1.0f;
        maxVal = 1.0f;

        if (doTests & DoTest_ATanH)
        {
            fprintf(stdout, "ATanH\n");
            f32 stdSec = run_speed_f32(string("stdlib atanh"), 15, "atanh", tests, minVal, maxVal, atanhf, 0.0f);
            run_speed_f32(string("mats atanh"), 15, "atanh", tests, minVal, maxVal, atanh32, stdSec);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Wide)
        {
            tests = 5000000;
            minVal = -100.0f;
            maxVal = 100.0f;

            if (doTests & DoTest_Cos)
            {
                fprintf(stdout, "Cos wide\n");
                f32 spdSecCos4x = run_speed_f32_4x(string("stdlib cos 4x"), 15, "cos_4x", tests, minVal, maxVal, cosf_4x, 0.0f);
                run_speed_f32_4x(string("mats cos 4x"), 15, "cos_4x", tests, minVal, maxVal, cos32_4x, spdSecCos4x);
                fprintf(stdout, "\n");
            }

            if (doTests & DoTest_Sin)
            {
                fprintf(stdout, "Sin wide\n");
                f32 spdSecSin4x = run_speed_f32_4x(string("stdlib sin 4x"), 15, "sin_4x", tests, minVal, maxVal, sinf_4x, 0.0f);
                run_speed_f32_4x(string("mats sin 4x"), 15, "sin_4x", tests, minVal, maxVal, sin32_4x, spdSecSin4x);
                fprintf(stdout, "\n");
            }

            if (doTests & DoTest_Tan)
            {
                fprintf(stdout, "Tan wide\n");
                f32 spdSecTan4x = run_speed_f32_4x(string("stdlib tan 4x"), 15, "tan_4x", tests, minVal, maxVal, tanf_4x, 0.0f);
                run_speed_f32_4x(string("mats tan 4x"), 15, "tan_4x", tests, minVal, maxVal, tan32_4x, spdSecTan4x);
                fprintf(stdout, "\n");
            }

            minVal = -1.2f;
            maxVal = 1.2f;

            if (doTests & DoTest_ACos)
            {
                fprintf(stdout, "ACos wide\n");
                f32 spdSecACos4x = run_speed_f32_4x(string("stdlib acos 4x"), 15, "acos_4x", tests, minVal, maxVal, acosf_4x, 0.0f);
                run_speed_f32_4x(string("mats acos 4x"), 15, "acos_4x", tests, minVal, maxVal, acos32_4x, spdSecACos4x);
                run_speed_f32_4x(string("mats tcos 4x"), 15, "acos_4x", tests, minVal, maxVal, acos32_temp_4x, spdSecACos4x);
                fprintf(stdout, "\n");
            }

            if (doTests & DoTest_ASin)
            {
                fprintf(stdout, "ASin wide\n");
                f32 spdSecASin4x = run_speed_f32_4x(string("stdlib asin 4x"), 15, "asin_4x", tests, minVal, maxVal, asinf_4x, 0.0f);
                run_speed_f32_4x(string("mats asin 4x"), 15, "asin_4x", tests, minVal, maxVal, asin32_4x, spdSecASin4x);
                run_speed_f32_4x(string("mats tsin 4x"), 15, "asin_4x", tests, minVal, maxVal, asin32_temp_4x, spdSecASin4x);
                fprintf(stdout, "\n");
            }

            minVal = -100.0f;
            maxVal = 100.0f;

            if (doTests & DoTest_ATan)
            {
                fprintf(stdout, "ATan wide\n");
                f32 spdSecATan4x = run_speed_f32_4x(string("stdlib atan 4x"), 15, "atan_4x", tests, minVal, maxVal, atanf_4x, 0.0f);
                run_speed_f32_4x(string("mats atan 4x"), 15, "atan_4x", tests, minVal, maxVal, atan32_4x, spdSecATan4x);
                run_speed_f32_4x(string("mats ttan 4x"), 15, "atan_4x", tests, minVal, maxVal, atan32_temp_4x, spdSecATan4x);
                fprintf(stdout, "\n");
            }
        }
    }

    return 0;
}
