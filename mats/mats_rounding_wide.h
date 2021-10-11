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

internal f32
modulus32_nosafe(f32 num, f32 den)
{
    // NOTE(michiel): No 'nan', 'inf' or subnormal support, and den may not be zero.

    // NOTE(michiel): Returns 'num % den' where the returned value's sign is equal to the sign of the numerator,
    // and the magnitude less than the magnitude of the denominator.
    // This function computes the remainder from the division of numerator by denominator. Specifically, the return value
    // is numerator - n * denominator, where n is the quotient of numerator divided by denominator, rounded towards zero.
    // Thus, fmod (6.5, 2.3) returns 1.9, which is 6.5 minus 4.6.
    s32 hnum = (s32)MATS_F32U(num).u;
    s32 hden = (s32)MATS_F32U(den).u;

    u32 signNum = hnum & MATS_F32_SIGN_MASK;
    hnum &= MATS_F32_ABS_MASK; // |num|
    hden &= MATS_F32_ABS_MASK; // |den|

    if (hnum < hden) {
        return num;               // NOTE(michiel): |num| < |den|, return num
    }
    if (hnum == hden) {
        return MATS_F32U(signNum).f; // NOTE(michiel): |num| == |den|, return 0 with the correct sign
    }

    /* determine ix = ilogb(x) */
    s32 expNum = (s32)(hnum >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;

    /* determine iy = ilogb(y) */
    s32 expDen = (s32)(hden >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;

    /* set up {hnum,lx}, {hden,ly} and align den to num */
    hnum = 0x00800000 | (MATS_F32_MANT_MASK & hnum);
    hden = 0x00800000 | (MATS_F32_MANT_MASK & hden);

    /* fixed point fmod */
	s32 n = expNum - expDen;
	while (n--)
    {
	    s32 hdiff = hnum - hden;
	    if (hdiff < 0) {
            hnum = hnum + hnum;
        } else {
	    	if (hdiff == 0) { 		/* return sign(x)*0 */
                return MATS_F32U(signNum).f;
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
        return MATS_F32U(signNum).f;
    }

	while (hnum < 0x00800000) {		/* normalize x */
	    hnum = hnum + hnum;
        --expDen;
	}

    hnum = ((hnum - 0x00800000) | ((expDen + MATS_F32_EXP_BIAS) << MATS_F32_EXP_SHIFT));
    f32 result = MATS_F32U(signNum | hnum).f;
	return result;		/* exact output */
}

internal f64_2x
floor64_2x(f64_2x x)
{
    f64_2x result;
    result.md = _mm_floor_pd(x.md);
    return result;
}

internal f64_2x
ceil64_2x(f64_2x x)
{
    f64_2x result;
    result.md = _mm_ceil_pd(x.md);
    return result;
}

internal f64_2x
round64_2x(f64_2x x)
{
    f64_2x result;
    result.md = _mm_round_pd(x.md, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    return result;
}

internal f64_2x
trunc64_2x(f64_2x x)
{
    f64_2x result;
    result.md = _mm_round_pd(x.md, _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
    return result;
}
