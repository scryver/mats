#include <errno.h>
#include <time.h>
#include <math.h>
#include <x86intrin.h>
#include <xmmintrin.h>

#define rint32_sse    round32_sse
#define fmod32_nonsse modulus32_nonsse
#define fmod32_nosafe modulus32_nosafe
#define fmod32_4x     modulus32_4x

#define STR_FMT(x)   safe_truncate_to_s32(x.size), (char *)x.data
#include "../libberdip/src/common.h"
#include "../libberdip/src/multilane.h"

#define IEEE_754_2008_SNAN 1
#define MATS_F32_ABS_MASK  0x7FFFFFFF
#include "../mats/mats_common.h"
#include "../mats/mats_constants.h"

#define floor32     floor32_nonsse
#define ceil32      ceil32_nonsse
#define round32     round32_nonsse
#define trunc32     trunc32_nonsse
#define modulus32   modulus32_nonsse
#define remainder32 remainder32_nonsse
#define MATS_USE_SSE4 0
#include "../mats/mats_defines.h"
#include "../mats/mats_rounding.h"
#undef floor32
#undef ceil32
#undef round32
#undef trunc32
#undef modulus32
#undef remainder32

#define floor32  floor32_sse
#define ceil32   ceil32_sse
#define round32  round32_sse
#define trunc32  trunc32_sse
#define modulus32   modulus32_sse
#define remainder32 remainder32_sse
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

#include "../mats/mats_rounding4x.h"

#include "test_common.cpp"

s32 main(s32 argc, char **argv)
{
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

    fprintf(stdout, "Correctness\n\n");

    fprintf(stdout, "Floor\n");
    f32 stdSecFloor32 = call_comp_x(stdlib, floor, f, 0.0f);
    call_comp_x(mats, floor, 32_nonsse, stdSecFloor32);
    call_comp_x(matsse, floor, 32_sse, stdSecFloor32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Ceil\n");
    f32 stdSecCeil32 = call_comp_x(stdlib, ceil, f, 0.0f);
    call_comp_x(mats, ceil, 32_nonsse, stdSecCeil32);
    call_comp_x(matsse, ceil, 32_sse, stdSecCeil32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Round (expect round sse to 'misbehave')\n");
    f32 stdSecRound32 = call_comp_x(stdlib, round, f, 0.0f);
    call_comp_x(mats, round, 32_nonsse, stdSecRound32);
    call_comp_x(matsse, round, 32_sse, stdSecRound32);
    fprintf(stdout, "Round sse vs correct call 'rintf'\n");
    f32 stdSecRint32 = call_comp_x(stdlib, rint, f, 0.0f);
    call_comp_x(matssse, rint, 32_sse, stdSecRint32 );
    fprintf(stdout, "\n");

    fprintf(stdout, "Truncate\n");
    f32 stdSecTruncate32 = call_comp_x(stdlib, trunc, f, 0.0f);
    call_comp_x(mats, trunc, 32_nonsse, stdSecTruncate32);
    call_comp_x(matsse, trunc, 32_sse, stdSecTruncate32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Mod\n");
    f32 stdSecMod32 = call_comp_x2(stdlib, fmod, f, 0.0f);
    call_comp_x2(mats, fmod, 32_nonsse, stdSecMod32);
    call_comp_x2(mats, fmod, 32_nosafe, stdSecMod32);
    call_comp_x2_4x(mats, fmod, 32_4x, stdSecMod32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Rem\n");
    f32 stdSecRem32 = call_comp_x2(stdlib, remainder, f, 0.0f);
    call_comp_x2(mats, remainder, 32_nonsse, stdSecRem32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Speed\n\n");

    fprintf(stdout, "Floor\n");
    f32 spdSecFloor32 = call_spd(stdlib, floor, f, 0.0f);
    call_spd(mats, floor, 32_nonsse, spdSecFloor32);
    call_spd(matsse, floor, 32_sse, spdSecFloor32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Ceil\n");
    f32 spdSecCeil32 = call_spd(stdlib, ceil, f, 0.0f);
    call_spd(mats, ceil, 32_nonsse, spdSecCeil32);
    call_spd(matsse, ceil, 32_sse, spdSecCeil32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Round\n");
    f32 spdSecRound32 = call_spd(stdlib, round, f, 0.0f);
    call_spd(mats, round, 32_nonsse, spdSecRound32);
    call_spd(matsse, round, 32_sse, spdSecRound32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Rint\n");
    f32 spdSecRint32 = call_spd(stdlib, rint, f, 0.0f);
    call_spd(matsse, rint, 32_sse, spdSecRint32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Truncate\n");
    f32 spdSecTrunc32 = call_spd(stdlib, trunc, f, 0.0f);
    call_spd(mats, trunc, 32_nonsse, spdSecTrunc32);
    call_spd(matsse, trunc, 32_sse, spdSecTrunc32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Mod\n");
    f32 spdSecMod32 = call_spd2(stdlib, fmod, f, 0.0f);
    call_spd2(mats, fmod, 32_nonsse, spdSecMod32);
    call_spd2(mats, fmod, 32_nosafe, spdSecMod32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Rem\n");
    f32 spdSecRem32 = call_spd2(stdlib, remainder, f, 0.0f);
    call_spd2(mats, remainder, 32_nonsse, spdSecRem32);
    fprintf(stdout, "\n");

    return 0;
}
