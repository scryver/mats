
#ifndef MATS_USE_SSE4
#define MATS_USE_SSE4  0
#endif

#ifndef MATS_USE_SSE2
#define MATS_USE_SSE2  MATS_USE_SSE4
#endif

#define MATS_USE_SSE  (MATS_USE_SSE2 || MATS_USE_SSE4)

#ifndef IEEE_754_2008_SNAN
#define IEEE_754_2008_SNAN 1
#endif

#ifndef MATS_HAVE_FAST_FMA
#define MATS_HAVE_FAST_FMA 0
#endif

#if MATS_USE_SSE
union WideMath
{
    __m128  m;
    __m128i mi;
    __m128d md;
    f32 e[4];
    u32 u[4];
    f64 ed[2];
    u64 ud[2];
};
#endif

#define MATS_F32_SIGN_MASK    0x80000000
#define MATS_F32_EXP_MASK     0x7F800000
#define MATS_F32_MANT_MASK    0x007FFFFF

#define MATS_F32_SIGN_SHIFT   31
#define MATS_F32_EXP_SHIFT    23

#define MATS_F32_ABS_MASK     0x7FFFFFFF
#define MATS_F32_EXP_BIAS     127
#define MATS_F32_EXP_MAX      127
#define MATS_F32_EXP_MIN     -126

#define MATS_F64_SIGN_MASK    0x8000000000000000ULL
#define MATS_F64_EXP_MASK     0x7FF0000000000000ULL
#define MATS_F64_MANT_MASK    0x000FFFFFFFFFFFFFULL

#define MATS_F64_SIGN_SHIFT   63
#define MATS_F64_EXP_SHIFT    52

#define MATS_F64_ABS_MASK     0x7FFFFFFFFFFFFFFFULL
#define MATS_F64_EXP_BIAS     1023
#define MATS_F64_EXP_MAX      1023
#define MATS_F64_EXP_MIN     -1022
