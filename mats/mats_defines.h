
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
