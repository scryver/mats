#include <x86intrin.h>
#include <math.h>

#define NO_INTRINSICS 1
#define MATS_USE_SSE4 0
#include "../libberdip/src/common.h"
#include "../libberdip/src/multilane.h"

#include "../mats/mats.h"
#include "../mats/mats_elem4x.h"
#include "../mats/mats_trig4x.h"

internal f32
tan_kernel_logic(f32 x, s32 mod)
{
	s32 ix = (s32)u32f32(x).u;
    u32 hx = ix & 0x7FFFFFFF;
    if (hx < 0x31800000)
    {
        if ((s32)x == 0) {
            if ((hx | (mod + 1)) == 0) {
                return 1.0f;
            } else {
                return (mod == 1) ? 2.0f : 3.0f;
            }
        }
    }
    if (hx >= 0x3F2CA140)
    {
        if (ix < 0) {
            x = -x;
        }
        x *= 2.0f;
    }

    if (hx >= 0x3F2CA140)
    {
        return x;
    }

    if (mod == 1) {
        return 4.0f;
    } else {
        return 5.0f;
    }
}

internal f32_4x
tan_kernel_logic_4x(f32_4x x, f32_4x mod)
{
    f32_4x hx = x & S32_4x(0x7FFFFFFF);

    f32_4x specialMask = s32_4x_less(hx, S32_4x(0x31800000));
    specialMask = s32_4x_and(s32_4x_equal(s32_4x_from_f32_trunc(x), zero_f32_4x()), specialMask);
    f32_4x special = zero_f32_4x();
    {
        f32_4x hxModP1Mask = s32_4x_equal(s32_4x_or(hx, s32_4x_add(mod, S32_4x(1))), zero_f32_4x());
        hxModP1Mask = s32_4x_and(hxModP1Mask, specialMask);
        special = F32_4x(1.0f) & hxModP1Mask;
        f32_4x workMask = and_not(specialMask, hxModP1Mask);
        f32_4x modIsOneMask = s32_4x_equal(mod, S32_4x(1));
#if MATS_USE_SSE4
        f32_4x modIsOne; modIsOne.m = _mm_blendv_ps(F32_4x(3.0f).m, F32_4x(2.0f).m, modIsOneMask.m);
#else
        f32_4x modIsOne = select(F32_4x(3.0f), modIsOneMask, F32_4x(2.0f));
#endif
        special = special | (modIsOne & workMask);
    }

    f32_4x hxGreat = s32_4x_less(S32_4x(0x3F2CA140), hx);
    f32_4x xLessZero = s32_4x_less(x, zero_f32_4x());
    xLessZero = s32_4x_and(hxGreat, s32_4x_and(xLessZero, S32_4x(MATS_F32_SIGN_MASK)));
    x = x ^ xLessZero;
    f32_4x xPre = F32_4x(2.0f) * x;
#if MATS_USE_SSE4
    x.m = _mm_blendv_ps(x.m, xPre.m, hxGreat.m);
#else
    x = select(x, hxGreat, xPre);
#endif
    f32_4x result = F32_4x(4.0f);

    f32_4x greatRes = x;

    f32_4x nonModdedMask = s32_4x_equal(mod, S32_4x(1));
    f32_4x modded = F32_4x(5.0f);

#if MATS_USE_SSE4
    result.m = _mm_blendv_ps(modded.m, result.m, nonModdedMask.m);
    result.m = _mm_blendv_ps(result.m, greatRes.m, hxGreat.m);
    result.m = _mm_blendv_ps(result.m, special.m, specialMask.m);
#else
    result = select(modded, nonModdedMask, result);
    result = select(result, hxGreat, greatRes);
    result = select(result, specialMask, special);
#endif

    return result;
}

int main(int argc, char **argv)
{
    f32 res1 = tan_kernel_logic(0.0f, -1);
    f32 res2 = tan_kernel_logic(0.0f, 1);
    f32 res3 = tan_kernel_logic(0.0f, 0);
    f32 res4 = tan_kernel_logic(0.5f, 1);
    f32 res5 = tan_kernel_logic(0.5f, -1);
    f32 res6 = tan_kernel_logic(5.0f, 1);
    f32 res7 = tan_kernel_logic(-5.0f, 1);
    fprintf(stdout, "Tan kernel(  0, -1): %f\n", res1);
    fprintf(stdout, "Tan kernel(  0,  1): %f\n", res2);
    fprintf(stdout, "Tan kernel(  0,  0): %f\n", res3);
    fprintf(stdout, "Tan kernel(0.5,  1): %f\n", res4);
    fprintf(stdout, "Tan kernel(0.5, -1): %f\n", res5);
    fprintf(stdout, "Tan kernel(  5,  1): %f\n", res6);
    fprintf(stdout, "Tan kernel( -5,  1): %f\n", res7);

    f32_4x res4_1 = tan_kernel_logic_4x(F32_4x(0.0f), S32_4x(-1));
    f32_4x res4_2 = tan_kernel_logic_4x(F32_4x(0.0f), S32_4x(1));
    f32_4x res4_3 = tan_kernel_logic_4x(F32_4x(0.0f), S32_4x(0));
    f32_4x res4_4 = tan_kernel_logic_4x(F32_4x(0.5f), S32_4x(1));
    f32_4x res4_5 = tan_kernel_logic_4x(F32_4x(0.5f), S32_4x(-1));
    f32_4x res4_6 = tan_kernel_logic_4x(F32_4x(5.0f), S32_4x(1));
    f32_4x res4_7 = tan_kernel_logic_4x(F32_4x(-5.0f), S32_4x(1));
    f32_4x res4_8 = tan_kernel_logic_4x(F32_4x(0.0f, 0.0f, 0.0f, 0.5f), S32_4x(-1, 1, 0, 1));
    f32_4x res4_9 = tan_kernel_logic_4x(F32_4x(0.5f, 5.0f, -5.0f, -0.5f), S32_4x(-1, 1, 1, 1));
    fprintf(stdout, "Tan kernel4(  0, -1): %f %f %f %f\n", res4_1.e[0], res4_1.e[1], res4_1.e[2], res4_1.e[3]);
    fprintf(stdout, "Tan kernel4(  0,  1): %f %f %f %f\n", res4_2.e[0], res4_2.e[1], res4_2.e[2], res4_2.e[3]);
    fprintf(stdout, "Tan kernel4(  0,  0): %f %f %f %f\n", res4_3.e[0], res4_3.e[1], res4_3.e[2], res4_3.e[3]);
    fprintf(stdout, "Tan kernel4(0.5,  1): %f %f %f %f\n", res4_4.e[0], res4_4.e[1], res4_4.e[2], res4_4.e[3]);
    fprintf(stdout, "Tan kernel4(0.5, -1): %f %f %f %f\n", res4_5.e[0], res4_5.e[1], res4_5.e[2], res4_5.e[3]);
    fprintf(stdout, "Tan kernel4(  5,  1): %f %f %f %f\n", res4_6.e[0], res4_6.e[1], res4_6.e[2], res4_6.e[3]);
    fprintf(stdout, "Tan kernel4( -5,  1): %f %f %f %f\n", res4_7.e[0], res4_7.e[1], res4_7.e[2], res4_7.e[3]);
    fprintf(stdout, "Tan kernel4: %f %f %f %f\n", res4_8.e[0], res4_8.e[1], res4_8.e[2], res4_8.e[3]);
    fprintf(stdout, "Tan kernel4: %f %f %f %f\n", res4_9.e[0], res4_9.e[1], res4_9.e[2], res4_9.e[3]);

    return 0;
}
