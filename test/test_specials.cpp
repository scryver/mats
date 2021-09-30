#include <errno.h>
#include <time.h>
#include <math.h>
#include <x86intrin.h>
#include <xmmintrin.h>

#define STR_FMT(x)   safe_truncate_to_s32(x.size), (char *)x.data
#include "../libberdip/src/common.h"
#include "../libberdip/src/multilane.h"

#define MATS_USE_SSE4 1
#include "../mats/mats.h"
#include "../mats/mats4x.h"

#include "test_common.cpp"
//#include "../mats/mats.cpp"

s32 main(s32 argc, char **argv)
{
    f32 pZeroF32 = 0.0f;
    f32 mZeroF32 = -0.0f;
    f32 pInfF32  = F32_INF;
    f32 mInfF32  = -F32_INF;
    f32 pNaNF32  = F32_NAN;
    f32 mNaNF32  = -F32_NAN;

    f32 f32specs[] = {
        pZeroF32, mZeroF32,
        pInfF32, mInfF32,
        pNaNF32, mNaNF32,
    };

    char *spaces = "                                                  ";
    u32 maxNameSize = 10;

    // TODO(michiel): 4x
#define spec_f32(func, name) \
String str##name = string(#name); \
for (u32 index = 0; index < array_count(f32specs); ++index) { \
f32 inputVal = f32specs[index]; \
f32 outputVal = func(inputVal); \
fprintf(stdout, "%.*s%.*s: %8g => %8g\n", (s32)minimum(maxNameSize, str##name.size), (char *)str##name.data, (s32)(maxNameSize - minimum(maxNameSize, str##name.size)), spaces, inputVal, outputVal); \
}

#define spec_f32_f32(func, name) \
String str##name = string(#name); \
for (u32 index = 0; index < array_count(f32specs); ++index) { \
f32 inputVal = f32specs[index]; \
f32 outputVal1 = func(inputVal, 1.2f); \
fprintf(stdout, "%.*s%.*s: (%8g,      1.2) => %8g\n", (s32)minimum(maxNameSize, str##name.size), (char *)str##name.data, (s32)(maxNameSize - minimum(maxNameSize, str##name.size)), spaces, inputVal, outputVal1); \
f32 outputVal2 = func(1.2f, inputVal); \
fprintf(stdout, "%.*s%.*s: (     1.2, %8g) => %8g\n", (s32)minimum(maxNameSize, str##name.size), (char *)str##name.data, (s32)(maxNameSize - minimum(maxNameSize, str##name.size)), spaces, inputVal, outputVal2); \
for (u32 indexB = 0; indexB < array_count(f32specs); ++indexB) { \
f32 inputB = f32specs[indexB]; \
if (indexB != index) { \
f32 outputB = func(inputVal, inputB); \
fprintf(stdout, "%.*s%.*s: (%8g, %8g) => %8g\n", (s32)minimum(maxNameSize, str##name.size), (char *)str##name.data, (s32)(maxNameSize - minimum(maxNameSize, str##name.size)), spaces, inputVal, inputB, outputB); \
} \
} \
}

    spec_f32(absolute32, Abs);
    spec_f32(floor32, Floor);
    spec_f32(ceil32, Ceil);
    spec_f32(round32, Round);
    spec_f32(trunc32, Trunc);
    spec_f32_f32(modulus32, Mod);
    spec_f32_f32(remainder32, Rem);

    spec_f32(sqrt32, Sqrt);
    spec_f32_f32(hypot32, Hypot);
    spec_f32(exp32, Exp);
    spec_f32(exp2_32, Exp2);
    spec_f32(log32, Log);
    spec_f32(log2_32, Log2);
    spec_f32(log10_32, Log10);
    spec_f32(expm1_32, ExpM1);
    spec_f32(log1p32, Log1P);
    spec_f32(log1p_fast32, Log1Pfast);
    spec_f32_f32(pow32, Pow);

    spec_f32(cos32, Cos);
    spec_f32(sin32, Sin);
    spec_f32(tan32, Tan);
    spec_f32(acos32, ACos);
    spec_f32(asin32, ASin);
    spec_f32(atan32, ATan);
    spec_f32_f32(atan2_32, ATan2);
    spec_f32(cosh32, CosH);
    spec_f32(sinh32, SinH);
    spec_f32(tanh32, TanH);
    spec_f32(acosh32, ACosH);
    spec_f32(asinh32, ASinH);
    spec_f32(atanh32, ATanH);

    return 0;
}