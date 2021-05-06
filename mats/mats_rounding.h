//
// NOTE(michiel): Floor/Ceil/Round/Truncate
//

internal f32
floor32(f32 x)
{
    f32 result;
#if MATS_USE_SSE4
    WideMath m;
    m.m = _mm_floor_ss(_mm_set_ss(x), _mm_set_ss(x));
    result = m.e[0];
#else
    s32 i0 = (s32)u32f32(x).u;
    u32 ix = i0 & 0x7FFFFFFF;
    s32 j0 = (ix >> 23) - 0x7F;

    if (j0 < 23)
    {
        if (j0 < 0)
        {
            if ((gHugeF32 + x) > 0.0f)
            {
                if (i0 >= 0) {
                    i0 = 0;
                } else if (!FLT_UWORD_IS_ZERO(ix)) {
                    i0 = 0xBF800000;
                }
            }
        }
        else
        {
            u32 i = 0x007FFFFF >> j0;
            if ((i0 & i) == 0) {
                return x;
            }
            if ((gHugeF32 + x) > 0.0f) {
                if (i0 < 0) {
                    i0 += 0x00800000 >> j0;
                }
                i0 &= ~i;
            }
        }
        result = u32f32((u32)i0).f;
    }
    else
    {
        if (!FLT_UWORD_IS_FINITE(ix)) {
            result = x + x;
        } else {
            result = x;
        }
    }
#endif
    return result;
}

internal f32
ceil32(f32 x)
{
    f32 result;
#if MATS_USE_SSE4
    WideMath m;
    m.m = _mm_ceil_ss(_mm_set_ss(x), _mm_set_ss(x));
    result = m.e[0];
#else
    s32 i0 = (s32)u32f32(x).u;
    u32 ix = i0 & 0x7FFFFFFF;
    s32 j0 = (ix >> 23) - 0x7F;

    if (j0 < 23)
    {
        if (j0 < 0)
        {
            if ((gHugeF32 + x) > 0.0f)
            {
                if (i0 < 0) {
                    i0 = 0x80000000;
                } else if (!FLT_UWORD_IS_ZERO(ix)) {
                    i0 = 0x3F800000;
                }
            }
        }
        else
        {
            u32 i = 0x007FFFFF >> j0;
            if ((i0 & i) == 0) {
                return x;
            }
            if ((gHugeF32 + x) > 0.0f) {
                if (i0 > 0) {
                    i0 += 0x00800000 >> j0;
                }
                i0 &= ~i;
            }
        }
        result = u32f32((u32)i0).f;
    }
    else
    {
        if (!FLT_UWORD_IS_FINITE(ix)) {
            result = x + x;
        } else {
            result = x;
        }
    }
#endif
    return result;
}

internal f32
round32(f32 x)
{
    f32 result;
#if MATS_USE_SSE4
    // NOTE(michiel): Does round-to-even to break ties (0.5f exact), so it will differ from the 'standard' c/c++ roundf
    WideMath m;
    m.m = _mm_round_ss(_mm_set_ss(x), _mm_set_ss(x), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    result = m.e[0];
#else
    // TODO(michiel): Should we do here the same round-to-even instead of round-away-from-zero to break ties?
    u32 w = u32f32(x).u;

    s32 exponentLess127 = (s32)((w & 0x7F800000) >> 23) - 127;
    if (exponentLess127 < 23)
    {
        if (exponentLess127 < 0) {
            w &= 0x80000000;
            if (exponentLess127 == -1) {
                // NOTE(michiel): Result is +/- 1.0
                w |= (127 << 23);
            }
        }
        else
        {
            u32 exponentMask = 0x007FFFFF >> exponentLess127;
            if ((w & exponentMask) == 0) {
                return x;
            }
            w += 0x00400000 >> exponentLess127;
            w &= ~exponentMask;
        }
        result = u32f32(w).f;
    }
    else
    {
        if (exponentLess127 == 128) {
            result = x + x;
        } else {
            result = x;
        }
    }
#endif
    return result;
}

internal f32
trunc32(f32 x)
{
    f32 result;
#if MATS_USE_SSE4
    WideMath m;
    m.m = _mm_round_ss(_mm_set_ss(x), _mm_set_ss(x), _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
    result = m.e[0];
#else
#if 0
    // TODO(michiel): Or should we just do this?
    result = (f32)(s32)x;
#else
    u32 w = u32f32(x).u;
    s32 signBit = w & 0x80000000;

    s32 exponentLess127 = ((w & 0x7F800000) >> 23) - 127;
    if (exponentLess127 < 23)
    {
        if (exponentLess127 < 0)
        {
            // NOTE(michiel): -1 < x < 1, so result is +0 or -0
            result = u32f32((u32)signBit).f;
        }
        else
        {
            result = u32f32((u32)signBit | (w & ~(0x007fffff >> exponentLess127))).f;
        }
    }
    else
    {
        if (exponentLess127 == 128) {
            /* x is NaN or infinite. */
            result = x + x;
        } else {
            /* All bits in the fraction field are relevant. */
            result = x;
        }
    }
#endif
#endif
    return result;
}
