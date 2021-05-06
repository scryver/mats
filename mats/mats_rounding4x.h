//
// NOTE(michiel): Floor/Ceil/Round/Truncate 4x
//

internal f32_4x
floor32_4x(f32_4x x)
{
    f32_4x result;
    result.m = _mm_floor_ps(x.m);
    return result;
}

internal f32_4x
ceil32_4x(f32_4x x)
{
    f32_4x result;
    result.m = _mm_ceil_ps(x.m);
    return result;
}

internal f32_4x
round32_4x(f32_4x x)
{
    f32_4x result;
    result.m = _mm_round_ps(x.m, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    return result;
}

internal f32_4x
trunc32_4x(f32_4x x)
{
    f32_4x result;
    result.m = _mm_round_ps(x.m, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
    return result;
}
