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

#define floor32  floor32_nonsse
#define ceil32   ceil32_nonsse
#define round32  round32_nonsse
#define trunc32  trunc32_nonsse
#define MATS_USE_SSE4 0
#include "../mats/mats_defines.h"
#include "../mats/mats_rounding.h"
#undef floor32
#undef ceil32
#undef round32
#undef trunc32

#define floor32  floor32_sse
#define ceil32   ceil32_sse
#define round32  round32_sse
#define trunc32  trunc32_sse
#undef  MATS_USE_SSE4
#define MATS_USE_SSE4 1
#include "../mats/mats_defines.h"
#include "../mats/mats_rounding.h"
#undef floor32
#undef ceil32
#undef round32
#undef trunc32

#include "test_common.cpp"

s32 main(s32 argc, char **argv)
{
    u32 tests = 200000000;
    f32 minVal = -12.0f;
    f32 maxVal = 123.8f;

    fprintf(stdout, "Test round (note that the sse implementation does round-to-even to break ties)\n");
    fprintf(stdout, "  -10.5 %8a %8a | %g %g\n", roundf(-10.5f), round32_sse(-10.5f), roundf(-10.5f), round32_sse(-10.5f));
    fprintf(stdout, "   -0.5 %8a %8a | %g %g\n", roundf( -0.5f), round32_sse( -0.5f), roundf( -0.5f), round32_sse( -0.5f));
    fprintf(stdout, "    0.5 %8a %8a | %g %g\n", roundf(  0.5f), round32_sse(  0.5f), roundf(  0.5f), round32_sse(  0.5f));

    fprintf(stdout, "Correctness\n\n");

    fprintf(stdout, "Floor\n");
    f32 stdSecFloor32 = run_comp_f32(string("stdlib floor"), 15, "floor", tests, minVal, maxVal, floorf, floorf, 0.0f);
    run_comp_f32(string("mats floor"), 15, "floor", tests, minVal, maxVal, floorf, floor32_nonsse, stdSecFloor32);
    run_comp_f32(string("mats floor sse"), 15, "floor", tests, minVal, maxVal, floorf, floor32_sse, stdSecFloor32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Ceil\n");
    f32 stdSecCeil32 = run_comp_f32(string("stdlib ceil"), 15, "ceil", tests, minVal, maxVal, ceilf, ceilf, 0.0f);
    run_comp_f32(string("mats ceil"), 15, "ceil", tests, minVal, maxVal, ceilf, ceil32_nonsse, stdSecCeil32);
    run_comp_f32(string("mats ceil sse"), 15, "ceil", tests, minVal, maxVal, ceilf, ceil32_sse, stdSecCeil32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Round (expect round sse to 'misbehave')\n");
    f32 stdSecRound32 = run_comp_f32(string("stdlib round"), 15, "round", tests, minVal, maxVal, roundf, roundf, 0.0f);
    run_comp_f32(string("mats round"), 15, "round", tests, minVal, maxVal, roundf, round32_nonsse, stdSecRound32);
    run_comp_f32(string("mats round sse"), 15, "round", tests, minVal, maxVal, roundf, round32_sse, stdSecRound32);
    fprintf(stdout, "Round sse vs correct call 'rintf'\n");
    run_comp_f32(string("mats round sse"), 15, "round", tests, minVal, maxVal, rintf, round32_sse, stdSecRound32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Truncate\n");
    f32 stdSecTruncate32 = run_comp_f32(string("stdlib trunc"), 15, "trunc", tests, minVal, maxVal, truncf, truncf, 0.0f);
    run_comp_f32(string("mats trunc"), 15, "trunc", tests, minVal, maxVal, truncf, trunc32_nonsse, stdSecTruncate32);
    run_comp_f32(string("mats trunc sse"), 15, "trunc", tests, minVal, maxVal, truncf, trunc32_sse, stdSecTruncate32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Speed\n\n");

    fprintf(stdout, "Floor\n");
    f32 spdSecFloor32 = run_speed_f32(string("stdlib floor"), 15, "floor", tests, minVal, maxVal, floorf, 0.0f);
    run_speed_f32(string("mats floor"), 15, "floor", tests, minVal, maxVal, floor32_nonsse, spdSecFloor32);
    run_speed_f32(string("mats floor sse"), 15, "floor", tests, minVal, maxVal, floor32_sse, spdSecFloor32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Ceil\n");
    f32 spdSecCeil32 = run_speed_f32(string("stdlib ceil"), 15, "ceil", tests, minVal, maxVal, ceilf, 0.0f);
    run_speed_f32(string("mats ceil"), 15, "ceil", tests, minVal, maxVal, ceil32_nonsse, spdSecCeil32);
    run_speed_f32(string("mats ceil sse"), 15, "ceil", tests, minVal, maxVal, ceil32_sse, spdSecCeil32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Round\n");
    f32 spdSecRound32 = run_speed_f32(string("stdlib round"), 15, "round", tests, minVal, maxVal, roundf, 0.0f);
    run_speed_f32(string("mats round"), 15, "round", tests, minVal, maxVal, round32_nonsse, spdSecRound32);
    run_speed_f32(string("mats round sse"), 15, "round", tests, minVal, maxVal, round32_sse, spdSecRound32);
    fprintf(stdout, "\n");

    fprintf(stdout, "Truncate\n");
    f32 spdSecTruncate32 = run_speed_f32(string("stdlib trunc"), 15, "trunc", tests, minVal, maxVal, truncf, 0.0f);
    run_speed_f32(string("mats trunc"), 15, "trunc", tests, minVal, maxVal, trunc32_nonsse, spdSecTruncate32);
    run_speed_f32(string("mats trunc sse"), 15, "trunc", tests, minVal, maxVal, trunc32_sse, spdSecTruncate32);
    fprintf(stdout, "\n");

    return 0;
}
