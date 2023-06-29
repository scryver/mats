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
#include "../mats/mats_constants.h"
#include "../mats/mats_common.h"
#include "../mats/mats_rounding.h"
#include "../mats/mats_elem.h"
#include "../mats/mats_elem_ext.h"
#include "../mats/mats_trig.h"
#include "../mats/mats_wide.h"

#include "../src/arm_sincos.cpp"

//#include "musl_temp.cpp"
//#include "newlib_temp.cpp"

#include "test_common.cpp"

func void
sincos_wrap64(f64 val, f64 *s, f64 *c)
{
    SinCos64 result = sincos64(val);
    *s = result.sin;
    *c = result.cos;
}

func void
sincos_wrap32(f32 val, f32 *s, f32 *c)
{
    SinCos32 result = sincos32(val);
    *s = result.sin;
    *c = result.cos;
}

func f64
sincos_match(f64 val)
{
    return square(sin(val)) + square(cos(val));
}

func f32
sincos_matchf(f32 val)
{
    return square(sinf(val)) + square(cosf(val));
}

func f32
sincos_match32(f32 val)
{
    return square(sin32(val)) + square(cos32(val));
}

func f32
sincos_match_cs_32(f32 val)
{
    SinCos32 cs = sincos32(val);
    return square(cs.sin) + square(cs.cos);
}

func f32
sincos_matchf_arm(f32 val)
{
    return square(arm_sinf(val)) + square(arm_cosf(val));
}

func f32
sincos_match_pi(f32 val)
{
    return square(sin_pi(val)) + square(cos_pi(val));
}

func f64
sincos_match64(f64 val)
{
    return square(sin64(val)) + square(cos64(val));
}

func f64
sincos_match_cs_64(f64 val)
{
    SinCos64 cs = sincos64(val);
    return square(cs.sin) + square(cs.cos);
}

func f32_4x
sinf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = sinf(x.e[0]);
    result.e[1] = sinf(x.e[1]);
    result.e[2] = sinf(x.e[2]);
    result.e[3] = sinf(x.e[3]);
    return result;
}

func f32_4x
cosf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = cosf(x.e[0]);
    result.e[1] = cosf(x.e[1]);
    result.e[2] = cosf(x.e[2]);
    result.e[3] = cosf(x.e[3]);
    return result;
}

func f32_4x
tanf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = tanf(x.e[0]);
    result.e[1] = tanf(x.e[1]);
    result.e[2] = tanf(x.e[2]);
    result.e[3] = tanf(x.e[3]);
    return result;
}

func f32_4x
acosf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = acosf(x.e[0]);
    result.e[1] = acosf(x.e[1]);
    result.e[2] = acosf(x.e[2]);
    result.e[3] = acosf(x.e[3]);
    return result;
}

func f32_4x
acos32_temp_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = acos32(x.e[0]);
    result.e[1] = acos32(x.e[1]);
    result.e[2] = acos32(x.e[2]);
    result.e[3] = acos32(x.e[3]);
    return result;
}

func f32_4x
asinf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = asinf(x.e[0]);
    result.e[1] = asinf(x.e[1]);
    result.e[2] = asinf(x.e[2]);
    result.e[3] = asinf(x.e[3]);
    return result;
}

func f32_4x
asin32_temp_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = asin32(x.e[0]);
    result.e[1] = asin32(x.e[1]);
    result.e[2] = asin32(x.e[2]);
    result.e[3] = asin32(x.e[3]);
    return result;
}

func f32_4x
atanf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = atanf(x.e[0]);
    result.e[1] = atanf(x.e[1]);
    result.e[2] = atanf(x.e[2]);
    result.e[3] = atanf(x.e[3]);
    return result;
}

func f32_4x
atan2f_4x(f32_4x y, f32_4x x)
{
    f32_4x result;
    result.e[0] = atan2f(y.e[0], x.e[0]);
    result.e[1] = atan2f(y.e[1], x.e[1]);
    result.e[2] = atan2f(y.e[2], x.e[2]);
    result.e[3] = atan2f(y.e[3], x.e[3]);
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

        BEGIN_TEST(doTests, cos, call_comp);
        call_comp(mats, cos, 32, stdSec);
        call_comp(arm, cos, f_arm, stdSec);
        call_comp(dip, cos, _pi, stdSec);
        call_comp_4x(mats4, cos, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, sin, call_comp);
        call_comp(mats, sin, 32, stdSec);
        call_comp(arm, sin, f_arm, stdSec);
        call_comp(dip, sin, _pi, stdSec);
        call_comp_4x(mats4, sin, 32_4x, stdSec);
        END_TEST();

        if ((doTests & (DoTest_sin | DoTest_cos)) == (DoTest_sin | DoTest_cos))
        {
            fprintf(stdout, "sincos32\n");
            f32 stdSec = call_comp_r2(stdlib, sincos, f, 0.0f);
            call_comp_r2(mats, sincos, _wrap32, stdSec);
            fprintf(stdout, "\n");

            fprintf(stdout, "Sin Cos matching\n");
            stdSec = call_comp(stdlib, sincos_match, f, 0.0f);
            call_comp(mats, sincos_match, 32, stdSec);
            call_comp(mats, sincos_match, _cs_32, stdSec);
            call_comp(arm, sincos_match, f_arm, stdSec);
            call_comp(dip, sincos_match, _pi, stdSec);
            // TODO(michiel): Add 4x
            fprintf(stdout, "\n");

#if 0
            fprintf(stdout, "SinCos\n");
            f32 stdSecSinCos32 = run_comp_f32(string("stdlib sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, sincosf, 0.0f);
            run_comp_f32(string("mats sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, sincos32, stdSecSinCos32);
            fprintf(stdout, "\n");
#endif
        }

        BEGIN_TEST(doTests, tan, call_comp);
        call_comp(mats, tan, 32, stdSec);
        call_comp_4x(mats4, tan, 32_4x, stdSec);
        END_TEST();

        minVal = -1.2f;
        maxVal =  1.2f;

        BEGIN_TEST(doTests, acos, call_comp);
        call_comp(mats, acos, 32, stdSec);
        call_comp_4x(mats4, acos, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, asin, call_comp);
        call_comp(mats, asin, 32, stdSec);
        call_comp_4x(mats4, asin, 32_4x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, atan, call_comp);
        call_comp(mats, atan, 32, stdSec);
        call_comp_4x(mats4, atan, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, atan2, call_comp2);
        call_comp2(mats, atan2, _32, stdSec);
        call_comp2_4x(mats4, atan2, _32_4x, stdSec);
        END_TEST();

        minVal = -22.0f;
        maxVal = 22.0f;

        BEGIN_TEST(doTests, cosh, call_comp);
        call_comp(mats, cosh, 32, stdSec);
        //call_comp_4x(mats4, cosh, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, sinh, call_comp);
        call_comp(mats, sinh, 32, stdSec);
        //call_comp_4x(mats4, sinh, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, tanh, call_comp);
        call_comp(mats, tanh, 32, stdSec);
        //call_comp_4x(mats4, tanh, 32_4x, stdSec);
        END_TEST();

        minVal = 1.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, acosh, call_comp);
        call_comp(mats, acosh, 32, stdSec);
        //call_comp_4x(mats4, acosh, 32_4x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, asinh, call_comp);
        call_comp(mats, asinh, 32, stdSec);
        //call_comp_4x(mats4, asinh, 32_4x, stdSec);
        END_TEST();

        minVal = -1.0f;
        maxVal = 1.0f;

        BEGIN_TEST(doTests, atanh, call_comp);
        call_comp(mats, atanh, 32, stdSec);
        //call_comp_4x(mats4, atanh, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, cos, call_comp64);
        call_comp64(mats, cos, 64, stdSec);
        //call_comp64_2x(mats4, cos, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, sin, call_comp64);
        call_comp64(mats, sin, 64, stdSec);
        //call_comp64_2x(mats4, sin, 64_2x, stdSec);
        END_TEST();

        if ((doTests & (DoTest_sin | DoTest_cos)) == (DoTest_sin | DoTest_cos))
        {
            fprintf(stdout, "sincos64\n");
            f32 stdSec = call_comp64_r2(stdlib, sincos, , 0.0f);
            call_comp64_r2(mats, sincos, _wrap64, stdSec);
            fprintf(stdout, "\n");

            fprintf(stdout, "Sin Cos matching\n");
            stdSec = call_comp64(stdlib, sincos_match, , 0.0f);
            call_comp64(mats, sincos_match, 64, stdSec);
            call_comp64(mats, sincos_match, _cs_64, stdSec);
            // TODO(michiel): Add 4x
            fprintf(stdout, "\n");

#if 0
            fprintf(stdout, "SinCos\n");
            f64 stdSecSinCos64 = run_comp_f64(string("stdlib sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, sincosf, 0.0f);
            run_comp_f64(string("mats sincos"), 15, "sincos", tests, minVal, maxVal, sincosf, sincos64, stdSecSinCos64);
            fprintf(stdout, "\n");
#endif
        }

        BEGIN_TEST64(doTests, tan, call_comp64);
        call_comp64(mats, tan, 64, stdSec);
        //call_comp64(musl, tan, _musl, stdSec);
        //call_comp64(newlib, tan, _newlib, stdSec);
        //call_comp64_2x(mats4, tan, 64_2x, stdSec);
        END_TEST();

        minVal = -1.2f;
        maxVal =  1.2f;

        BEGIN_TEST64(doTests, acos, call_comp64);
        call_comp64(mats, acos, 64, stdSec);
        //call_comp64_2x(mats4, acos, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, asin, call_comp64);
        call_comp64(mats, asin, 64, stdSec);
        //call_comp64_2x(mats4, asin, 64_2x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST64(doTests, atan, call_comp64);
        call_comp64(mats, atan, 64, stdSec);
        //call_comp64_2x(mats4, atan, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, atan2, call_comp64_2);
        call_comp64_2(mats, atan2, _64, stdSec);
        //call_comp64_2_2x(mats4, atan2, _64_2x, stdSec);
        END_TEST();

        minVal = -22.0f;
        maxVal = 22.0f;

        BEGIN_TEST64(doTests, cosh, call_comp64);
        call_comp64(mats, cosh, 64, stdSec);
        //call_comp64_2x(mats4, cosh, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, sinh, call_comp64);
        call_comp64(mats, sinh, 64, stdSec);
        //call_comp64_2x(mats4, sinh, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, tanh, call_comp64);
        call_comp64(mats, tanh, 64, stdSec);
        //call_comp64_2x(mats4, tanh, 64_2x, stdSec);
        END_TEST();

        minVal = 1.0f;
        maxVal = 100.0f;

        BEGIN_TEST64(doTests, acosh, call_comp64);
        call_comp64(mats, acosh, 64, stdSec);
        //call_comp64_2x(mats4, acosh, 64_2x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST64(doTests, asinh, call_comp64);
        call_comp64(mats, asinh, 64, stdSec);
        //call_comp64_2x(mats4, asinh, 64_2x, stdSec);
        END_TEST();

        minVal = -1.0f;
        maxVal = 1.0f;

        BEGIN_TEST64(doTests, atanh, call_comp64);
        call_comp64(mats, atanh, 64, stdSec);
        //call_comp64_2x(mats4, atanh, 64_2x, stdSec);
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
        call_spd_4x(mats4, cos, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, sin, call_spd);
        call_spd(mats, sin, 32, stdSec);
        call_spd(arm, sin, f_arm, stdSec);
        call_spd(dip, sin, _pi, stdSec);
        call_spd_4x(mats4, sin, 32_4x, stdSec);
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
        call_spd_4x(mats4, tan, 32_4x, stdSec);
        END_TEST();

        minVal = -1.2f;
        maxVal = 1.2f;

        BEGIN_TEST(doTests, acos, call_spd);
        call_spd(mats, acos, 32, stdSec);
        call_spd_4x(mats4, acos, 32_4x, stdSec);
        call_spd_4x(tats4, acos, 32_temp_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, asin, call_spd);
        call_spd(mats, asin, 32, stdSec);
        call_spd_4x(mats4, asin, 32_4x, stdSec);
        call_spd_4x(tats4, asin, 32_temp_4x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, atan, call_spd);
        call_spd(mats, atan, 32, stdSec);
        call_spd_4x(mats4, atan, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, atan2, call_spd2);
        call_spd2(mats, atan2, _32, stdSec);
        call_spd2_4x(mats4, atan2, _32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, cosh, call_spd);
        call_spd(mats, cosh, 32, stdSec);
        //call_spd_4x(mats4, cosh, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, sinh, call_spd);
        call_spd(mats, sinh, 32, stdSec);
        //call_spd_4x(mats4, sinh, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, tanh, call_spd);
        call_spd(mats, tanh, 32, stdSec);
        //call_spd_4x(mats4, tanh, 32_4x, stdSec);
        END_TEST();

        minVal = 1.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, acosh, call_spd);
        call_spd(mats, acosh, 32, stdSec);
        //call_spd_4x(mats4, acosh, 32_4x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST(doTests, asinh, call_spd);
        call_spd(mats, asinh, 32, stdSec);
        //call_spd_4x(mats4, asinh, 32_4x, stdSec);
        END_TEST();

        minVal = -1.0f;
        maxVal = 1.0f;

        BEGIN_TEST(doTests, atanh, call_spd);
        call_spd(mats, atanh, 32, stdSec);
        //call_spd_4x(mats4, atanh, 32_4x, stdSec);
        END_TEST();

        //
        // NOTE(michiel): 64-bit
        //

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST64(doTests, cos, call_spd64);
        call_spd64(mats, cos, 64, stdSec);
        //call_spd64_2x(matsse, cos, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, sin, call_spd64);
        call_spd64(mats, sin, 64, stdSec);
        //call_spd64_2x(matsse, sin, 64_2x, stdSec);
        END_TEST();

#if 0
        if ((doTests & (DoTest_Sin | DoTest_Cos)) == (DoTest_Sin | DoTest_Cos))
        {
            fprintf(stdout, "SinCos\n");
            f64 spdSecSinCos64 = run_speed_f64(string("stdlib sincos"), 15, "sincos", tests, minVal, maxVal, sincos, 0.0f);
            run_speed_f64(string("mats sincos"), 15, "sincos", tests, minVal, maxVal, sincos64, spdSecSinCos64);
            fprintf(stdout, "\n");
        }
#endif

        BEGIN_TEST64(doTests, tan, call_spd64);
        call_spd64(mats, tan, 64, stdSec);
        //call_spd64(musl, tan, _musl, stdSec);
        //call_spd64(newlib, tan, _newlib, stdSec);
        //call_spd64_2x(matsse, tan, 64_2x, stdSec);
        END_TEST();

        minVal = -1.2f;
        maxVal = 1.2f;

        BEGIN_TEST64(doTests, acos, call_spd64);
        call_spd64(mats, acos, 64, stdSec);
        //call_spd64_2x(mats4, acos, 64_4x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, asin, call_spd64);
        call_spd64(mats, asin, 64, stdSec);
        //call_spd64_2x(mats4, asin, 64_2x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST64(doTests, atan, call_spd64);
        call_spd64(mats, atan, 64, stdSec);
        //call_spd64_2x(mats4, atan, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, atan2, call_spd64_2);
        call_spd64_2(mats, atan2, _64, stdSec);
        //call_spd64_2_2x(mats4, atan2, _64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, cosh, call_spd64);
        call_spd64(mats, cosh, 64, stdSec);
        //call_spd64_2x(mats4, cosh, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, sinh, call_spd64);
        call_spd64(mats, sinh, 64, stdSec);
        //call_spd64_2x(mats4, sinh, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, tanh, call_spd64);
        call_spd64(mats, tanh, 64, stdSec);
        //call_spd64_2x(mats4, tanh, 64_2x, stdSec);
        END_TEST();

        minVal = 1.0f;
        maxVal = 100.0f;

        BEGIN_TEST64(doTests, acosh, call_spd64);
        call_spd64(mats, acosh, 64, stdSec);
        //call_spd64_2x(mats4, acosh, 64_2x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST64(doTests, asinh, call_spd64);
        call_spd64(mats, asinh, 64, stdSec);
        //call_spd64_2x(mats4, asinh, 64_2x, stdSec);
        END_TEST();

        minVal = -1.0f;
        maxVal = 1.0f;

        BEGIN_TEST64(doTests, atanh, call_spd64);
        call_spd64(mats, atanh, 64, stdSec);
        //call_spd64_2x(mats4, atanh, 64_2x, stdSec);
        END_TEST();
    }

    if (doTests & DoTest_Wide)
    {
        tests = 5000000;
        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST_WIDE(doTests, cos, call_spd_4x);
        call_spd_4x(mats4, cos, 32_4x, stdSec);
        call_spd_4x(dip4, cos, _f32_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, sin, call_spd_4x);
        call_spd_4x(mats4, sin, 32_4x, stdSec);
        call_spd_4x(dip4, sin, _f32_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, tan, call_spd_4x);
        call_spd_4x(mats4, tan, 32_4x, stdSec);
        END_TEST();

        minVal = -1.2f;
        maxVal = 1.2f;

        BEGIN_TEST_WIDE(doTests, acos, call_spd_4x);
        call_spd_4x(mats4, acos, 32_4x, stdSec);
        call_spd_4x(tats4, acos, 32_temp_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, asin, call_spd_4x);
        call_spd_4x(mats4, asin, 32_4x, stdSec);
        call_spd_4x(tats4, asin, 32_temp_4x, stdSec);
        END_TEST();

        minVal = -100.0f;
        maxVal = 100.0f;

        BEGIN_TEST_WIDE(doTests, atan, call_spd_4x);
        call_spd_4x(mats4, atan, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, atan2, call_spd2_4x);
        call_spd2_4x(mats4, atan2, _32_4x, stdSec);
        END_TEST();
    }

    return 0;
}
