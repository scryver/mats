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
atanf_4x(f32_4x x)
{
    f32_4x result;
    result.e[0] = atanf(x.e[0]);
    result.e[1] = atanf(x.e[1]);
    result.e[2] = atanf(x.e[2]);
    result.e[3] = atanf(x.e[3]);
    return result;
}

internal void
set_default_fp_behavior(void)
{
#define FLUSH_TO_ZERO_BIT (1 << 15)
#define ROUNDING_CONTROL_BITS (3 << 13)
#define PRECISION_MASK (1 << 12)
#define UNDERFLOW_MASK (1 << 11)
#define OVERFLOW_MASK (1 << 10)
#define DBZ_MASK (1 << 9)
#define DENORMAL_OP_MASK (1 << 8)
#define INVALID_OP_MASK (1 << 7)
#define DENORMALS_ARE_ZERO (1 << 6)
    u32 fpControlMask = (FLUSH_TO_ZERO_BIT |
                         // ROUNDING_CONTROL_BITS |
                         PRECISION_MASK |
                         UNDERFLOW_MASK |
                         OVERFLOW_MASK |
                         DBZ_MASK |
                         DENORMAL_OP_MASK |
                         INVALID_OP_MASK |
                         DENORMALS_ARE_ZERO);
    u32 desiredBits = fpControlMask;
    u32 oldControlBits = _mm_getcsr();
    u32 newControlBits = (oldControlBits & ~fpControlMask) | desiredBits;
    _mm_setcsr(newControlBits);
}

enum DoTestFlag
{
    DoTest_Cos  = 0x00000001,
    DoTest_Sin  = 0x00000002,
    DoTest_Tan  = 0x00000004,
    DoTest_ACos = 0x00000008,
    DoTest_ASin = 0x00000010,
    DoTest_ATan = 0x00000020,

    DoTest_FuncMask = 0x000000FF,

    DoTest_Comp  = 0x10000000,
    DoTest_Speed = 0x20000000,
    DoTest_Wide  = 0x40000000,

    DoTest_TypeMask = 0xF0000000,
};

s32 main(s32 argc, char **argv)
{
    //set_default_fp_behavior();

    u32 doTests = U32_MAX;

    if (argc > 1)
    {
        u32 tests = 0;
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
                    } else {
                        i_expect(*str == 'w');
                        tests |= DoTest_Wide;
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
            }
        }

        if (tests & DoTest_FuncMask) {
            doTests = (doTests & ~DoTest_FuncMask) | tests;
        }
        if (tests & DoTest_TypeMask) {
            doTests = (doTests & ~DoTest_TypeMask) | tests;
        }
    }

    // TODO(michiel): Arguments: select function and speed vs acc
    u32 tests = 50000000;
    f32 minVal = -F32_TAU; // -10.0f;
    f32 maxVal = F32_TAU; // 10.0f;

    if (doTests & DoTest_Comp)
    {
        fprintf(stdout, "\nCorrectness\n\n");

        if (doTests & DoTest_Cos)
        {
            fprintf(stdout, "Cos\n");
            f32 stdSecCos32 = run_comp_f32_x(string("stdlib cos"), 15, "cos", tests, minVal, maxVal, cos, cosf, 0.0f);
            run_comp_f32_x(string("mats cos"), 15, "cos", tests, minVal, maxVal, cos, cos32, stdSecCos32);
            run_comp_f32_x(string("mats cos prec"), 15, "cos", tests, minVal, maxVal, cos, cos32_prec, stdSecCos32);
            run_comp_f32_x(string("arm cos"), 15, "cos", tests, minVal, maxVal, cos, arm_cosf, stdSecCos32);
            run_comp_f32_x(string("dip cos"), 15, "cos", tests, minVal, maxVal, cos, cos_pi, stdSecCos32);
            run_comp_f32_4x_x(string("mats cos 4x"), 15, "cos_4x", tests, minVal, maxVal, cos, cos32_4x, stdSecCos32);
            run_comp_f32_4x_x(string("mats prec 4x"), 15, "cos_4x", tests, minVal, maxVal, cos, cos32_prec_4x, stdSecCos32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Sin)
        {
            fprintf(stdout, "Sin\n");
            f32 stdSecSin32 = run_comp_f32_x(string("stdlib sin"), 15, "sin", tests, minVal, maxVal, sin, sinf, 0.0f);
            run_comp_f32_x(string("mats sin"), 15, "sin", tests, minVal, maxVal, sin, sin32, stdSecSin32);
            run_comp_f32_x(string("mats sin prec"), 15, "sin", tests, minVal, maxVal, sin, sin32_prec, stdSecSin32);
            run_comp_f32_x(string("arm sin"), 15, "sin", tests, minVal, maxVal, sin, arm_sinf, stdSecSin32);
            run_comp_f32_x(string("dip sin"), 15, "sin", tests, minVal, maxVal, sin, sin_pi, stdSecSin32);
            run_comp_f32_4x_x(string("mats sin 4x"), 15, "sin_4x", tests, minVal, maxVal, sin, sin32_4x, stdSecSin32);
            run_comp_f32_4x_x(string("mats prec 4x"), 15, "sin_4x", tests, minVal, maxVal, sin, sin32_prec_4x, stdSecSin32);
            fprintf(stdout, "\n");
        }

        if ((doTests & (DoTest_Sin | DoTest_Cos)) == (DoTest_Sin | DoTest_Cos))
        {
            fprintf(stdout, "Sin Cos matching\n");
            f32 stdSCMatchSin32 = run_comp_f32_x(string("stdlib match"), 15, "match", tests, minVal, maxVal, sincos_match64, sincosf_match, 0.0f);
            run_comp_f32_x(string("mats match"), 15, "match", tests, minVal, maxVal, sincos_match64, sincos32_match, stdSCMatchSin32);
            run_comp_f32_x(string("mats match prc"), 15, "match", tests, minVal, maxVal, sincos_match64, sincos32_prec_match, stdSCMatchSin32);
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
        maxVal = 1.2f;

        if (doTests & DoTest_ACos)
        {
            fprintf(stdout, "ACos\n");
            f32 stdSecACos32 = run_comp_f32_x(string("stdlib acos"), 15, "acos", tests, minVal, maxVal, acos, acosf, 0.0f);
            run_comp_f32_x(string("mats acos"), 15, "acos", tests, minVal, maxVal, acos, acos32, stdSecACos32);
            run_comp_f32_x(string("mats temp"), 15, "acos", tests, minVal, maxVal, acos, acos32_temp, stdSecACos32);
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
            run_comp_f32_x(string("mats temp"), 15, "atan", tests, minVal, maxVal, atan, atan32_temp, stdSecATan32);
            run_comp_f32_4x_x(string("mats atan 4x"), 15, "atan_4x", tests, minVal, maxVal, atan, atan32_4x, stdSecATan32);
            fprintf(stdout, "\n");
        }
    }

    if (doTests & DoTest_Speed)
    {
        fprintf(stdout, "Speed\n\n");
        minVal = -F32_TAU; // -10.0f;
        maxVal = F32_TAU; // 10.0f;

        if (doTests & DoTest_Cos)
        {
            fprintf(stdout, "Cos\n");
            f32 spdSecCos32 = run_speed_f32(string("stdlib cos"), 15, "cos", tests, minVal, maxVal, cosf, 0.0f);
            run_speed_f32(string("mats cos"), 15, "cos", tests, minVal, maxVal, cos32, spdSecCos32);
            run_speed_f32(string("mats cos prec"), 15, "cos", tests, minVal, maxVal, cos32_prec, spdSecCos32);
            run_speed_f32(string("arm cos"), 15, "cos", tests, minVal, maxVal, arm_cosf, spdSecCos32);
            run_speed_f32(string("dip cos"), 15, "cos", tests, minVal, maxVal, cos_pi, spdSecCos32);
            run_speed_f32_4x(string("mats cos 4x"), 15, "cos", tests, minVal, maxVal, cos32_4x, spdSecCos32);
            run_speed_f32_4x(string("mats prec 4x"), 15, "cos", tests, minVal, maxVal, cos32_prec_4x, spdSecCos32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Sin)
        {
            fprintf(stdout, "Sin\n");
            f32 spdSecSin32 = run_speed_f32(string("stdlib sin"), 15, "sin", tests, minVal, maxVal, sinf, 0.0f);
            run_speed_f32(string("mats sin"), 15, "sin", tests, minVal, maxVal, sin32, spdSecSin32);
            run_speed_f32(string("mats sin prec"), 15, "sin", tests, minVal, maxVal, sin32_prec, spdSecSin32);
            run_speed_f32(string("arm sin"), 15, "sin", tests, minVal, maxVal, arm_sinf, spdSecSin32);
            run_speed_f32(string("dip sin"), 15, "sin", tests, minVal, maxVal, sin_pi, spdSecSin32);
            run_speed_f32_4x(string("mats sin 4x"), 15, "sin", tests, minVal, maxVal, sin32_4x, spdSecSin32);
            run_speed_f32_4x(string("mats prec 4x"), 15, "sin", tests, minVal, maxVal, sin32_prec_4x, spdSecSin32);
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
            run_speed_f32(string("mats temp"), 15, "acos", tests, minVal, maxVal, acos32_temp, spdSecACos32);
            run_speed_f32_4x(string("mats acos 4x"), 15, "acos", tests, minVal, maxVal, acos32_4x, spdSecACos32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_ASin)
        {
            fprintf(stdout, "ASin\n");
            f32 spdSecASin32 = run_speed_f32(string("stdlib asin"), 15, "asin", tests, minVal, maxVal, asinf, 0.0f);
            run_speed_f32(string("mats asin"), 15, "asin", tests, minVal, maxVal, asin32, spdSecASin32);
            run_speed_f32_4x(string("mats asin 4x"), 15, "asin", tests, minVal, maxVal, asin32_4x, spdSecASin32);
            fprintf(stdout, "\n");
        }

        minVal = -100.0f;
        maxVal = 100.0f;

        if (doTests & DoTest_ATan)
        {
            fprintf(stdout, "ATan\n");
            f32 spdSecATan32 = run_speed_f32(string("stdlib atan"), 15, "atan", tests, minVal, maxVal, atanf, 0.0f);
            run_speed_f32(string("mats atan"), 15, "atan", tests, minVal, maxVal, atan32, spdSecATan32);
            run_speed_f32(string("mats temp"), 15, "atan", tests, minVal, maxVal, atan32_temp, spdSecATan32);
            run_speed_f32_4x(string("mats atan 4x"), 15, "atan", tests, minVal, maxVal, atan32_4x, spdSecATan32);
            fprintf(stdout, "\n");
        }

        if (doTests & DoTest_Wide)
        {
            tests = 5000000;
            minVal = -F32_TAU; // -10.0f;
            maxVal = F32_TAU; // 10.0f;

            if (doTests & DoTest_Cos)
            {
                fprintf(stdout, "Cos wide\n");
                f32 spdSecCos4x = run_speed_f32_4x(string("stdlib cos 4x"), 15, "cos_4x", tests, minVal, maxVal, cosf_4x, 0.0f);
                run_speed_f32_4x(string("mats cos 4x"), 15, "cos_4x", tests, minVal, maxVal, cos32_4x, spdSecCos4x);
                run_speed_f32_4x(string("mats prec 4x"), 15, "cos_4x", tests, minVal, maxVal, cos32_prec_4x, spdSecCos4x);
                fprintf(stdout, "\n");
            }

            if (doTests & DoTest_Sin)
            {
                fprintf(stdout, "Sin wide\n");
                f32 spdSecSin4x = run_speed_f32_4x(string("stdlib sin 4x"), 15, "sin_4x", tests, minVal, maxVal, sinf_4x, 0.0f);
                run_speed_f32_4x(string("mats sin 4x"), 15, "sin_4x", tests, minVal, maxVal, sin32_4x, spdSecSin4x);
                run_speed_f32_4x(string("mats prec 4x"), 15, "sin_4x", tests, minVal, maxVal, sin32_prec_4x, spdSecSin4x);
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
                fprintf(stdout, "\n");
            }

            if (doTests & DoTest_ASin)
            {
                fprintf(stdout, "ASin wide\n");
                f32 spdSecASin4x = run_speed_f32_4x(string("stdlib asin 4x"), 15, "asin_4x", tests, minVal, maxVal, asinf_4x, 0.0f);
                run_speed_f32_4x(string("mats asin 4x"), 15, "asin_4x", tests, minVal, maxVal, asin32_4x, spdSecASin4x);
                fprintf(stdout, "\n");
            }

            minVal = -100.0f;
            maxVal = 100.0f;

            if (doTests & DoTest_ATan)
            {
                fprintf(stdout, "ATan wide\n");
                f32 spdSecATan4x = run_speed_f32_4x(string("stdlib atan 4x"), 15, "atan_4x", tests, minVal, maxVal, atanf_4x, 0.0f);
                run_speed_f32_4x(string("mats atan 4x"), 15, "atan_4x", tests, minVal, maxVal, atan32_4x, spdSecATan4x);
                fprintf(stdout, "\n");
            }
        }
    }

    return 0;
}
