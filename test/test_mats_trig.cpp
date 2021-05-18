#include <errno.h>
#include <time.h>
#include <math.h>
#include <x86intrin.h>
#include <xmmintrin.h>

#define cosf_arm  arm_cosf
#define sinf_arm  arm_sinf

#define NO_INTRINSICS 1
#define STR_FMT(x)   safe_truncate_to_s32(x.size), (char *)x.data
#include "../libberdip/src/common.h"
#include "../libberdip/src/multilane.h"
#include "../libberdip/src/trigonometry_v2.h"

#define MATS_USE_SSE4 1
#define MATS_USE_SSE2 1
#define IEEE_754_2008_SNAN 1
#include "../mats/mats_defines.h"
#include "../mats/mats_common.h"
#include "../mats/mats_constants.h"
#include "../mats/mats_elem.h"
#include "../mats/mats_elem_ext.h"
#include "../mats/mats_elem4x.h"
#include "../mats/mats_trig.h"
#include "../mats/mats_trig4x.h"

#include "../src/arm_sincos.cpp"

#include "test_common.cpp"

internal f64
sincos_match(f64 val)
{
    return 1.0;
}

internal f32
sincos_matchf(f32 val)
{
    return square(sinf(val)) + square(cosf(val));
}

internal f32
sincos_match32(f32 val)
{
    return square(sin32(val)) + square(cos32(val));
}

internal f32
sincos_matchf_arm(f32 val)
{
    return square(arm_sinf(val)) + square(arm_cosf(val));
}

internal f32
sincos_match_pi(f32 val)
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
    DoTest_cos        = 0x00000001,
    DoTest_sin        = 0x00000002,
    DoTest_tan        = 0x00000004,
    DoTest_acos       = 0x00000008,
    DoTest_asin       = 0x00000010,
    DoTest_atan       = 0x00000020,
    DoTest_atan2      = 0x00000040,
    DoTest_cosh       = 0x00000080,
    DoTest_sinh       = 0x00000100,
    DoTest_tanh       = 0x00000200,
    DoTest_acosh      = 0x00000400,
    DoTest_asinh      = 0x00000800,
    DoTest_atanh      = 0x00001000,
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
                tests |= DoTest_cos;
            } else if (strings_are_equal("sin", argv[index])) {
                tests |= DoTest_sin;
            } else if (strings_are_equal("tan", argv[index])) {
                tests |= DoTest_tan;
            } else if (strings_are_equal("acos", argv[index])) {
                tests |= DoTest_acos;
            } else if (strings_are_equal("asin", argv[index])) {
                tests |= DoTest_asin;
            } else if (strings_are_equal("atan", argv[index])) {
                tests |= DoTest_atan;
            } else if (strings_are_equal("atan2", argv[index])) {
                tests |= DoTest_atan2;
            } else if (strings_are_equal("cosh", argv[index])) {
                tests |= DoTest_cosh;
            } else if (strings_are_equal("sinh", argv[index])) {
                tests |= DoTest_sinh;
            } else if (strings_are_equal("tanh", argv[index])) {
                tests |= DoTest_tanh;
            } else if (strings_are_equal("acosh", argv[index])) {
                tests |= DoTest_acosh;
            } else if (strings_are_equal("asinh", argv[index])) {
                tests |= DoTest_asinh;
            } else if (strings_are_equal("atanh", argv[index])) {
                tests |= DoTest_atanh;
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

        BEGIN_TEST(doTests, cos, call_comp_x);
        call_comp_x(mats, cos, 32, stdSec);
        call_comp_x(arm, cos, f_arm, stdSec);
        call_comp_x(dip, cos, _pi, stdSec);
        call_comp_x_4x(mats, cos, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, sin, call_comp_x);
        call_comp_x(mats, sin, 32, stdSec);
        call_comp_x(arm, sin, f_arm, stdSec);
        call_comp_x(dip, sin, _pi, stdSec);
        call_comp_x_4x(mats, sin, 32_4x, stdSec);
        END_TEST();

        if ((doTests & (DoTest_Sin | DoTest_Cos)) == (DoTest_Sin | DoTest_Cos))
        {
            fprintf(stdout, "Sin Cos matching\n");
            f32 stdSec = call_comp_x(stdlib, sincos_match, f, 0.0f);
            call_comp_x(mats, sincos_match, 32, stdSec);
            call_comp_x(arm, sincos_match, f_arm, stdSec);
            call_comp_x(dip, sincos_match, _pi, stdSec);
            // TODO(michiel): Add 4x
            fprintf(stdout, "\n");

#if 0
            fprintf(stdout, "SinCos\n");
            f32 stdSecSinCos32 = run_comp_f32(string("stdlib sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, sincosf, 0.0f);
            run_comp_f32(string("mats sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, sincos32, stdSecSinCos32);
            fprintf(stdout, "\n");
#endif
        }

        BEGIN_TEST(doTests, tan, call_comp_x);
        call_comp_x(mats, tan, 32, stdSec);
        call_comp_x_4x(mats, tan, 32_4x, stdSec);
        END_TEST();

        minVal = -1.2f;
        maxVal =  1.2f;

        BEGIN_TEST(doTests, acos, call_comp_x);
        call_comp_x(mats, acos, 32, stdSec);
        call_comp_x_4x(mats, acos, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, asin, call_comp_x);
        call_comp_x(mats, asin, 32, stdSec);
        call_comp_x_4x(mats, asin, 32_4x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, atan, call_comp_x);
        call_comp_x(mats, atan, 32, stdSec);
        call_comp_x_4x(mats, atan, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, atan2, call_comp_x2);
        call_comp_x2(mats, atan2, _32, stdSec);
        END_TEST();

        minVal = -22.0f;
        maxVal = 22.0f;

        BEGIN_TEST(doTests, cosh, call_comp_x);
        call_comp_x(mats, cosh, 32, stdSec);
        //call_comp_x_4x(mats, cosh, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, sinh, call_comp_x);
        call_comp_x(mats, sinh, 32, stdSec);
        //call_comp_x_4x(mats, sinh, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, tanh, call_comp_x);
        call_comp_x(mats, tanh, 32, stdSec);
        //call_comp_x_4x(mats, tanh, 32_4x, stdSec);
        END_TEST();

        minVal = 1.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, acosh, call_comp_x);
        call_comp_x(mats, acosh, 32, stdSec);
        //call_comp_x_4x(mats, acosh, 32_4x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, asinh, call_comp_x);
        call_comp_x(mats, asinh, 32, stdSec);
        //call_comp_x_4x(mats, asinh, 32_4x, stdSec);
        END_TEST();

        minVal = -1.0f;
        maxVal = 1.0f;

        BEGIN_TEST(doTests, atanh, call_comp_x);
        call_comp_x(mats, atanh, 32, stdSec);
        //call_comp_x_4x(mats, atanh, 32_4x, stdSec);
        END_TEST();

    }

    if (doTests & DoTest_Speed)
    {
        fprintf(stdout, "Speed\n\n");
        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, cos, call_spd);
        call_spd(mats, cos, 32, stdSec);
        call_spd(arm, cos, f_arm, stdSec);
        call_spd(dip, cos, _pi, stdSec);
        call_spd_4x(mats, cos, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, sin, call_spd);
        call_spd(mats, sin, 32, stdSec);
        call_spd(arm, sin, f_arm, stdSec);
        call_spd(dip, sin, _pi, stdSec);
        call_spd_4x(mats, sin, 32_4x, stdSec);
        END_TEST();

#if 0
        if ((doTests & (DoTest_Sin | DoTest_Cos)) == (DoTest_Sin | DoTest_Cos))
        {
            fprintf(stdout, "SinCos\n");
            f32 spdSecSinCos32 = run_speed_f32(string("stdlib sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, 0.0f);
            run_speed_f32(string("mats sincos"), 15, "sincos", tests, minVal, maxVal, sincos32, spdSecSinCos32);
            fprintf(stdout, "\n");
        }
#endif

        BEGIN_TEST(doTests, tan, call_spd);
        call_spd(mats, tan, 32, stdSec);
        call_spd_4x(mats, tan, 32_4x, stdSec);
        END_TEST();

        minVal = -1.2f;
        maxVal = 1.2f;

        BEGIN_TEST(doTests, acos, call_spd);
        call_spd(mats, acos, 32, stdSec);
        call_spd_4x(mats, acos, 32_4x, stdSec);
        call_spd_4x(tats, acos, 32_temp_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, asin, call_spd);
        call_spd(mats, asin, 32, stdSec);
        call_spd_4x(mats, asin, 32_4x, stdSec);
        call_spd_4x(tats, asin, 32_temp_4x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, atan, call_spd);
        call_spd(mats, atan, 32, stdSec);
        call_spd_4x(mats, atan, 32_4x, stdSec);
        call_spd_4x(tats, atan, 32_temp_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, atan2, call_spd2);
        call_spd2(mats, atan2, _32, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, cosh, call_spd);
        call_spd(mats, cosh, 32, stdSec);
        //call_spd_4x(mats, cosh, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, sinh, call_spd);
        call_spd(mats, sinh, 32, stdSec);
        //call_spd_4x(mats, sinh, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, tanh, call_spd);
        call_spd(mats, tanh, 32, stdSec);
        //call_spd_4x(mats, tanh, 32_4x, stdSec);
        END_TEST();

        minVal = 1.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, acosh, call_spd);
        call_spd(mats, acosh, 32, stdSec);
        //call_spd_4x(mats, acosh, 32_4x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, asinh, call_spd);
        call_spd(mats, asinh, 32, stdSec);
        //call_spd_4x(mats, asinh, 32_4x, stdSec);
        END_TEST();

        minVal = -1.0f;
        maxVal = 1.0f;

        BEGIN_TEST(doTests, atanh, call_spd);
        call_spd(mats, atanh, 32, stdSec);
        //call_spd_4x(mats, atanh, 32_4x, stdSec);
        END_TEST();

        if (doTests & DoTest_Wide)
        {
            tests = 5000000;
            minVal = -100.0f;
            maxVal = 100.0f;

            BEGIN_TEST_WIDE(doTests, cos, call_spd_4x);
            call_spd_4x(mats, cos, 32_4x, stdSec);
            END_TEST();

            BEGIN_TEST_WIDE(doTests, sin, call_spd_4x);
            call_spd_4x(mats, sin, 32_4x, stdSec);
            END_TEST();

            BEGIN_TEST_WIDE(doTests, tan, call_spd_4x);
            call_spd_4x(mats, tan, 32_4x, stdSec);
            END_TEST();

            minVal = -1.2f;
            maxVal = 1.2f;

            BEGIN_TEST_WIDE(doTests, acos, call_spd_4x);
            call_spd_4x(mats, acos, 32_4x, stdSec);
            call_spd_4x(tats, acos, 32_temp_4x, stdSec);
            END_TEST();

            BEGIN_TEST_WIDE(doTests, asin, call_spd_4x);
            call_spd_4x(mats, asin, 32_4x, stdSec);
            call_spd_4x(tats, asin, 32_temp_4x, stdSec);
            END_TEST();

            minVal = -100.0f;
            maxVal = 100.0f;

            BEGIN_TEST_WIDE(doTests, atan, call_spd_4x);
            call_spd_4x(mats, atan, 32_4x, stdSec);
            call_spd_4x(tats, atan, 32_temp_4x, stdSec);
            END_TEST();
        }
    }

    return 0;
}
