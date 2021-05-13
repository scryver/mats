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

internal f32_4x
modulus32_4x(f32_4x num, f32_4x den)
{
    f32_4x result = num - (round32_4x(num / den) * den);
    return result;
}

internal f32
modulus32_nosafe(f32 num, f32 den)
{
    // NOTE(michiel): No 'nan', 'inf' support, and den may not be zero.

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

    if (hnum < hden) {
        return num;               // NOTE(michiel): |num| < |den|, return num
    }
    if (hnum == hden) {
        return u32f32(signNum).f; // NOTE(michiel): |num| == |den|, return 0 with the correct sign
    }

    /* determine ix = ilogb(x) */
    s32 expNum = (s32)(hnum >> 23) - 127;

    /* determine iy = ilogb(y) */
    s32 expDen = (s32)(hden >> 23) - 127;

    /* set up {hnum,lx}, {hden,ly} and align den to num */
    hnum = 0x00800000 | (0x007FFFFF & hnum);
    hden = 0x00800000 | (0x007FFFFF & hden);

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

    hnum = ((hnum - 0x00800000) | ((expDen + 127) << 23));
    f32 result = u32f32(signNum | hnum).f;
	return result;		/* exact output */
}
