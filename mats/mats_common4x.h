
internal f32_4x
select4x(f32_4x a, f32_4x mask, f32_4x b)
{
    // NOTE(michiel): if mask return b else a
    f32_4x result;
#if MATS_USE_SSE4
    result.m = _mm_blendv_ps(a.m, b.m, mask.m);
#else
    result.m = _mm_or_ps(_mm_andnot_ps(mask.m, a.m),
                         _mm_and_ps(mask.m, b.m));
#endif
    return result;
}

internal f64_2x
select2x(f64_2x a, f64_2x mask, f64_2x b)
{
    // NOTE(michiel): if mask return b else a
    f64_2x result;
#if MATS_USE_SSE4
    result.md = _mm_blendv_pd(a.md, b.md, mask.md);
#else
    result.md = _mm_or_pd(_mm_andnot_pd(mask.md, a.md),
                          _mm_and_pd(mask.md, b.md));
#endif
    return result;
}
