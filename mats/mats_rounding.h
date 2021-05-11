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

internal f32
modulus32(f32 num, f32 den)
{
    // NOTE(michiel): Returns 'num % den' where the returned value's sign is equal to the sign of the numerator,
    // and the magnitude less than the magnitude of the denominator.
    // This function computes the remainder from the division of numerator by denominator. Specifically, the return value
    // is numerator - n * denominator, where n is the quotient of numerator divided by denominator, rounded towards zero.
    // Thus, fmod (6.5, 2.3) returns 1.9, which is 6.5 minus 4.6.
    s32 hnum = (s32)u32f32(num).u;
    s32 hden = (s32)u32f32(den).u;

    u32 signNum = hnum & 0x80000000;
    hnum &= 0x7FFFFFFF; // |num|
    hden &= 0x7FFFFFFF; // |den|

    /* purge off exception values */
    if (FLT_UWORD_IS_ZERO(hden) ||
        !FLT_UWORD_IS_FINITE(hnum) ||
        FLT_UWORD_IS_NAN(hden)) {
        return (num * den) / (num * den);
    }

    if (hnum < hden) {
        return num;               // NOTE(michiel): |num| < |den|, return num
    }
    if (hnum == hden) {
        return u32f32(signNum).f; // NOTE(michiel): |num| == |den|, return 0 with the correct sign
    }

    /* Note: den cannot be zero if we reach here. */

    s32 i;
    /* determine ix = ilogb(x) */
    s32 expNum;
    if (FLT_UWORD_IS_SUBNORMAL(hnum))
    {
        for (expNum = -126, i = (hnum << 8);
             i > 0;
             i <<= 1)
        {
            --expNum;
        }
    }
    else
    {
        expNum = (s32)(hnum >> 23) - 127;
    }

    /* determine iy = ilogb(y) */
    s32 expDen;
    if (FLT_UWORD_IS_SUBNORMAL(hden))
    {
        for (expDen = -126, i = (hden << 8);
             i >= 0;
             i <<= 1)
        {
            --expDen;
        }
    }
    else
    {
        expDen = (s32)(hden >> 23) - 127;
    }

    /* set up {hnum,lx}, {hden,ly} and align den to num */
	if(expNum >= -126) {
	    hnum = 0x00800000 | (0x007FFFFF & hnum);
    } else {		/* subnormal x, shift x to normal */
	    s32 n = -126 - expNum;
	    hnum = hnum << n;
	}
    if (expDen >= -126) {
        hden = 0x00800000 | (0x007FFFFF & hden);
    } else {
        s32 n = -126 - expDen;
        hden = hden << n;
    }

    /* fixed point fmod */
	s32 n = expNum - expDen;
	while (n--)
    {
	    s32 hdiff = hnum - hden;
	    if (hdiff < 0) {
            hnum = hnum + hnum;
        } else {
	    	if (hdiff == 0) { 		/* return sign(x)*0 */
                return u32f32(signNum).f;
            }
            hnum = hdiff + hdiff;
	    }
	}

	s32 hdiff = hnum - hden;
	if (hdiff >= 0) {
        hnum = hdiff;
    }

    /* convert back to floating value and restore the sign */
	if (hnum == 0) {
        /* return sign(x)*0 */
        return u32f32(signNum).f;
    }

	while (hnum < 0x00800000) {		/* normalize x */
	    hnum = hnum + hnum;
        --expDen;
	}

    f32 result;
	if (expDen >= -126) {		/* normalize output */
        hnum = ((hnum - 0x00800000) | ((expDen + 127) << 23));
        result = u32f32(signNum | hnum).f;
	} else {		/* subnormal output */
	    /* If denormals are not supported, this code will generate a
	       zero representation.  */
	    s32 n = -126 - expDen;
	    hnum >>= n;
        result = u32f32(signNum | hnum).f;
	    //result *= 1.0f;		/* create necessary signal */
	}
	return result;		/* exact output */
}

internal f32
remainder32(f32 num, f32 den)
{
    // NOTE(michiel): Returns 'num % den' where the absolute value of the result is less than or equal to half the absolute value
    // of the denominator. returned value's sign is equal to the sign of the numerator,
    // This function is like modulus32 except that it rounds the internal quotient n to the nearest integer instead of towards
    // zero. For example, remainder (6.5, 2.3) returns -0.4, which is 6.5 minus 6.9.
    u32 hnum = u32f32(num).u;
    u32 hden = u32f32(den).u;

    u32 signNum = hnum & 0x80000000;
    hnum ^= signNum;    // |num|
    hden &= 0x7FFFFFFF; // |den|

    /* purge off exception values */
    if (FLT_UWORD_IS_ZERO(hden) ||
        !FLT_UWORD_IS_FINITE(hnum) ||
        FLT_UWORD_IS_NAN(hden))
    {
        return (num * den) / (num * den);
    }

    if (hden <= FLT_UWORD_HALF_MAX) {
        num = modulus32(num, den + den);
    }
    if ((hnum - hden) == 0) {
        return u32f32(signNum).f;
    }

    num = mats_fabs32(num);
    den = mats_fabs32(den);

    if (hden < 0x01000000)
    {
        if ((num + num) > den) {
            num -= den;
            if ((num + num) >= den) {
                num -= den;
            }
        }
    }
    else
    {
        f32 denHalf = 0.5f * den;
        if (num > denHalf) {
            num -= den;
            if (num >= denHalf) {
                num -= den;
            }
        }
    }

    hnum = u32f32(num).u;
    f32 result = u32f32(signNum ^ hnum).f;
    return result;
}
