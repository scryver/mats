#include <errno.h>
#include <time.h>
#include <math.h>
#include <x86intrin.h>
#include <xmmintrin.h>

#define rint32_sse    round32_sse
#define modulusf      fmodf
#define modulus       fmod

#define STR_FMT(x)   safe_truncate_to_s32(x.size), (char *)x.data
#include "../libberdip/src/common.h"
#include "../libberdip/src/multilane.h"

#define IEEE_754_2008_SNAN 1
#define MATS_F32_ABS_MASK  0x7FFFFFFF
#define MATS_F64_ABS_MASK  0x7FFFFFFFFFFFFFFFULL
#include "../mats/mats_common.h"
#include "../mats/mats_constants.h"

#define floor32     floor32_nonsse
#define ceil32      ceil32_nonsse
#define round32     round32_nonsse
#define trunc32     trunc32_nonsse
#define modulus32   modulus32_nonsse
#define remainder32 remainder32_nonsse
#define floor64     floor64_nonsse
#define ceil64      ceil64_nonsse
#define round64     round64_nonsse
#define trunc64     trunc64_nonsse
#define modulus64   modulus64_nonsse
#define remainder64 remainder64_nonsse
#define MATS_USE_SSE4 0
#include "../mats/mats_defines.h"
#include "../mats/mats_rounding.h"
#undef floor32
#undef ceil32
#undef round32
#undef trunc32
#undef modulus32
#undef remainder32
#undef floor64
#undef ceil64
#undef round64
#undef trunc64
#undef modulus64
#undef remainder64

#define floor32  floor32_sse
#define ceil32   ceil32_sse
#define round32  round32_sse
#define trunc32  trunc32_sse
#define modulus32   modulus32_sse
#define remainder32 remainder32_sse
#define floor64  floor64_sse
#define ceil64   ceil64_sse
#define round64  round64_sse
#define trunc64  trunc64_sse
#define modulus64   modulus64_sse
#define remainder64 remainder64_sse
#undef  MATS_USE_SSE4
#define MATS_USE_SSE4 1
#include "../mats/mats_defines.h"
#include "../mats/mats_rounding.h"
#undef floor32
#undef ceil32
#undef round32
#undef trunc32
#undef modulus32
#undef remainder32
#undef floor64
#undef ceil64
#undef round64
#undef trunc64
#undef modulus64
#undef remainder64

#include "../mats/mats_rounding4x.h"

#include "test_common.cpp"

WIDE_FUNC_FROM_F32(floorf)
WIDE_FUNC_FROM_F32(ceilf)
WIDE_FUNC_FROM_F32(roundf)
WIDE_FUNC_FROM_F32(truncf)
WIDE_FUNC_FROM_F64(floor)
WIDE_FUNC_FROM_F64(ceil)
WIDE_FUNC_FROM_F64(round)
WIDE_FUNC_FROM_F64(trunc);

enum DoTestFlag
{
    DoTest_floor      = 0x00000001,
    DoTest_ceil       = 0x00000002,
    DoTest_round      = 0x00000004,
    DoTest_trunc      = 0x00000008,
    DoTest_modulus    = 0x00000010,
    DoTest_remainder  = 0x00000020,
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
            } else if (strings_are_equal("floor", argv[index])) {
                tests |= DoTest_floor;
            } else if (strings_are_equal("ceil", argv[index])) {
                tests |= DoTest_ceil;
            } else if (strings_are_equal("round", argv[index])) {
                tests |= DoTest_round;
            } else if (strings_are_equal("trunc", argv[index])) {
                tests |= DoTest_trunc;
            } else if (strings_are_equal("mod", argv[index])) {
                tests |= DoTest_modulus;
            } else if (strings_are_equal("rem", argv[index])) {
                tests |= DoTest_remainder;
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

    u32 tests = 200000000;
    f32 minVal = -12.0f;
    f32 maxVal = 123.8f;

    u32 testsA = 20000;
    u32 testsB = 10000;
    f32 minValA = -1200.0f;
    f32 maxValA = 1230.8f;
    f32 minValB = -100.0f;
    f32 maxValB = 100.0f;

    fprintf(stdout, "Test round (note that the sse implementation does round-to-even to break ties)\n");
    fprintf(stdout, "  -10.5 %8a %8a | %g %g\n", roundf(-10.5f), round32_sse(-10.5f), roundf(-10.5f), round32_sse(-10.5f));
    fprintf(stdout, "   -0.5 %8a %8a | %g %g\n", roundf( -0.5f), round32_sse( -0.5f), roundf( -0.5f), round32_sse( -0.5f));
    fprintf(stdout, "    0.5 %8a %8a | %g %g\n", roundf(  0.5f), round32_sse(  0.5f), roundf(  0.5f), round32_sse(  0.5f));

    fprintf(stdout, "\n");

    if (doTests & DoTest_Comp)
    {
        fprintf(stdout, "Correctness\n\n");

        BEGIN_TEST(doTests, floor, call_comp);
        call_comp(mats, floor, 32_nonsse, stdSec);
        call_comp(matsse, floor, 32_sse, stdSec);
        call_comp_4x(mats4, floor, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, ceil, call_comp);
        call_comp(mats, ceil, 32_nonsse, stdSec);
        call_comp(matsse, ceil, 32_sse, stdSec);
        call_comp_4x(mats4, ceil, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, round, call_comp);
        fprintf(stdout, "Round (expect round sse to 'misbehave')\n");
        call_comp(mats, round, 32_nonsse, stdSec);
        call_comp(matsse, round, 32_sse, stdSec);
        call_comp_4x(mats4, round, 32_4x, stdSec);

        fprintf(stdout, "Round sse vs correct call 'rintf'\n");
        f32 stdSecRint32 = call_comp(stdlib, rint, f, 0.0f);
        call_comp(matsse, rint, 32_sse, stdSecRint32 );
        END_TEST();

        BEGIN_TEST(doTests, trunc, call_comp);
        call_comp(mats, trunc, 32_nonsse, stdSec);
        call_comp(matsse, trunc, 32_sse, stdSec);
        call_comp_4x(mats4, trunc, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, modulus, call_comp2);
        call_comp2(mats, modulus, 32_nonsse, stdSec);
        call_comp2(matsse, modulus, 32_nosafe, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, remainder, call_comp2);
        call_comp2(mats, remainder, 32_nonsse, stdSec);
        END_TEST();


        BEGIN_TEST64(doTests, floor, call_comp64);
        call_comp64(mats, floor, 64_nonsse, stdSec);
        call_comp64(matsse, floor, 64_sse, stdSec);
        call_comp64_2x(mats4, floor, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, ceil, call_comp64);
        call_comp64(mats, ceil, 64_nonsse, stdSec);
        call_comp64(matsse, ceil, 64_sse, stdSec);
        call_comp_4x(mats4, ceil, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, round, call_comp64);
        fprintf(stdout, "Round (expect round sse to 'misbehave')\n");
        call_comp64(mats, round, 64_nonsse, stdSec);
        call_comp64(matsse, round, 64_sse, stdSec);
        call_comp64_2x(mats4, round, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, trunc, call_comp64);
        call_comp64(mats, trunc, 64_nonsse, stdSec);
        call_comp64(matsse, trunc, 64_sse, stdSec);
        call_comp64_2x(mats4, trunc, 64_2x, stdSec);
        END_TEST();

    }

    if (doTests & DoTest_Speed)
    {
        fprintf(stdout, "Speed\n\n");

        BEGIN_TEST(doTests, floor, call_spd);
        call_spd(mats, floor, 32_nonsse, stdSec);
        call_spd(matsse, floor, 32_sse, stdSec);
        call_spd_4x(mats4, floor, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, ceil, call_spd);
        call_spd(mats, ceil, 32_nonsse, stdSec);
        call_spd(matsse, ceil, 32_sse, stdSec);
        call_spd_4x(mats4, ceil, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, round, call_spd);
        call_spd(mats, round, 32_nonsse, stdSec);
        call_spd(matsse, round, 32_sse, stdSec);
        call_spd_4x(mats4, round, 32_4x, stdSec);

        fprintf(stdout, "Rint\n");
        f32 spdSecRint32 = call_spd(stdlib, rint, f, 0.0f);
        call_spd(matsse, rint, 32_sse, spdSecRint32);
        END_TEST();

        BEGIN_TEST(doTests, trunc, call_spd);
        call_spd(mats, trunc, 32_nonsse, stdSec);
        call_spd(matsse, trunc, 32_sse, stdSec);
        call_spd_4x(mats4, trunc, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, modulus, call_spd2);
        call_spd2(mats, modulus, 32_nonsse, stdSec);
        call_spd2(fats, modulus, 32_nosafe, stdSec);
        END_TEST();

        BEGIN_TEST(doTests, remainder, call_spd2);
        call_spd2(mats, remainder, 32_nonsse, stdSec);
        END_TEST();


        BEGIN_TEST64(doTests, floor, call_spd64);
        call_spd64(mats, floor, 64_nonsse, stdSec);
        call_spd64(matsse, floor, 64_sse, stdSec);
        call_spd64_2x(mats4, floor, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, ceil, call_spd64);
        call_spd64(mats, ceil, 64_nonsse, stdSec);
        call_spd64(matsse, ceil, 64_sse, stdSec);
        call_spd64_2x(mats4, ceil, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, round, call_spd64);
        call_spd64(mats, round, 64_nonsse, stdSec);
        call_spd64(matsse, round, 64_sse, stdSec);
        call_spd64_2x(mats4, round, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST64(doTests, trunc, call_spd64);
        call_spd64(mats, trunc, 64_nonsse, stdSec);
        call_spd64(matsse, trunc, 64_sse, stdSec);
        call_spd64_2x(mats4, trunc, 64_2x, stdSec);
        END_TEST();

    }

    if (doTests & DoTest_Wide)
    {
        BEGIN_TEST_WIDE(doTests, floor, call_spd_4x);
        call_spd_4x(mats4, floor, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, ceil, call_spd_4x);
        call_spd_4x(mats4, ceil, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, round, call_spd_4x);
        call_spd_4x(mats4, round, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE(doTests, trunc, call_spd_4x);
        call_spd_4x(mats4, trunc, 32_4x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE64(doTests, floor, call_spd64_2x);
        call_spd64_2x(mats4, floor, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE64(doTests, ceil, call_spd64_2x);
        call_spd64_2x(mats4, ceil, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE64(doTests, round, call_spd64_2x);
        call_spd64_2x(mats4, round, 64_2x, stdSec);
        END_TEST();

        BEGIN_TEST_WIDE64(doTests, trunc, call_spd64_2x);
        call_spd64_2x(mats4, trunc, 64_2x, stdSec);
        END_TEST();
    }

    return 0;
}
