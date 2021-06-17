//
// NOTE(michiel): Floor/Ceil/Round/Truncate/Modulus/Remainder
//

//
// NOTE(michiel): 32-bit
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
    s32 i0 = MATS_S32_FROM_F32(x);
    u32 ix = i0 & MATS_F32_ABS_MASK;
    s32 j0 = (ix >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;

    if (j0 < MATS_F32_EXP_SHIFT)
    {
        if (j0 < 0)
        {
            if ((gHugeF32 + x) > 0.0f)
            {
                if (i0 >= 0) {
                    i0 = 0;
                } else if (!MATS_F32_UWORD_IS_ZERO(ix)) {
                    i0 = 0xBF800000;
                }
            }
        }
        else
        {
            u32 i = MATS_F32_MANT_MASK >> j0;
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
        result = MATS_F32_FROM_S32(i0);
    }
    else
    {
        if (!MATS_F32_UWORD_IS_FINITE(ix)) {
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
    s32 i0 = MATS_S32_FROM_F32(x);
    u32 ix = i0 & MATS_F32_ABS_MASK;
    s32 j0 = (ix >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;

    if (j0 < MATS_F32_EXP_SHIFT)
    {
        if (j0 < 0)
        {
            if ((gHugeF32 + x) > 0.0f)
            {
                if (i0 < 0) {
                    i0 = MATS_F32_SIGN_MASK;
                } else if (!MATS_F32_UWORD_IS_ZERO(ix)) {
                    i0 = 0x3F800000;
                }
            }
        }
        else
        {
            u32 i = MATS_F32_MANT_MASK >> j0;
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
        result = MATS_F32_FROM_S32(i0);
    }
    else
    {
        if (!MATS_F32_UWORD_IS_FINITE(ix)) {
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
    u32 w = MATS_U32_FROM_F32(x);

    s32 exponentLess127 = (s32)((w & MATS_F32_EXP_MASK) >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;
    if (exponentLess127 < MATS_F32_EXP_SHIFT)
    {
        if (exponentLess127 < 0) {
            w &= MATS_F32_SIGN_MASK;
            if (exponentLess127 == -1) {
                // NOTE(michiel): Result is +/- 1.0
                w |= (MATS_F32_EXP_BIAS << MATS_F32_EXP_SHIFT);
            }
        }
        else
        {
            u32 exponentMask = MATS_F32_MANT_MASK >> exponentLess127;
            if ((w & exponentMask) == 0) {
                return x;
            }
            w += 0x00400000 >> exponentLess127;
            w &= ~exponentMask;
        }
        result = MATS_F32_FROM_U32(w);
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
    u32 w = MATS_U32_FROM_F32(x);
    s32 signBit = w & MATS_F32_SIGN_MASK;

    s32 exponentLess127 = ((w & MATS_F32_EXP_MASK) >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;
    if (exponentLess127 < MATS_F32_EXP_SHIFT)
    {
        if (exponentLess127 < 0)
        {
            // NOTE(michiel): -1 < x < 1, so result is +0 or -0
            result = MATS_F32_FROM_S32(signBit);
        }
        else
        {
            result = MATS_F32_FROM_U32((u32)signBit | (w & ~(MATS_F32_MANT_MASK >> exponentLess127)));
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

    u32 signNum = hnum & MATS_F32_SIGN_MASK;
    hnum &= MATS_F32_ABS_MASK; // |num|
    hden &= MATS_F32_ABS_MASK; // |den|

    /* purge off exception values */
    if (MATS_F32_UWORD_IS_ZERO(hden) ||
        !MATS_F32_UWORD_IS_FINITE(hnum) ||
        MATS_F32_UWORD_IS_NAN(hden)) {
        return (num * den) / (num * den);
    }

    if (hnum < hden) {
        return num;               // NOTE(michiel): |num| < |den|, return num
    }
    if (hnum == hden) {
        return MATS_F32_FROM_U32(signNum); // NOTE(michiel): |num| == |den|, return 0 with the correct sign
    }

    /* Note: den cannot be zero if we reach here. */

    s32 i;
    /* determine ix = ilogb(x) */
    s32 expNum;
    if (MATS_F32_UWORD_IS_SUBNORMAL(hnum))
    {
        for (expNum = MATS_F32_EXP_MIN, i = (hnum << 8);
             i > 0;
             i <<= 1)
        {
            --expNum;
        }
    }
    else
    {
        expNum = (s32)(hnum >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;
    }

    /* determine iy = ilogb(y) */
    s32 expDen;
    if (MATS_F32_UWORD_IS_SUBNORMAL(hden))
    {
        for (expDen = MATS_F32_EXP_MIN, i = (hden << 8);
             i >= 0;
             i <<= 1)
        {
            --expDen;
        }
    }
    else
    {
        expDen = (s32)(hden >> MATS_F32_EXP_SHIFT) - MATS_F32_EXP_BIAS;
    }

    /* set up {hnum,lx}, {hden,ly} and align den to num */
	if(expNum >= MATS_F32_EXP_MIN) {
	    hnum = 0x00800000 | (MATS_F32_MANT_MASK & hnum);
    } else {		/* subnormal x, shift x to normal */
	    s32 n = MATS_F32_EXP_MIN - expNum;
	    hnum = hnum << n;
	}
    if (expDen >= MATS_F32_EXP_MIN) {
        hden = 0x00800000 | (MATS_F32_MANT_MASK & hden);
    } else {
        s32 n = MATS_F32_EXP_MIN - expDen;
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
	if (expDen >= MATS_F32_EXP_MIN) {		/* normalize output */
        hnum = ((hnum - 0x00800000) | ((expDen + MATS_F32_EXP_BIAS) << MATS_F32_EXP_SHIFT));
        result = u32f32(signNum | hnum).f;
	} else {		/* subnormal output */
	    /* If denormals are not supported, this code will generate a
	       zero representation.  */
	    s32 n = MATS_F32_EXP_MIN - expDen;
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

    u32 signNum = hnum & MATS_F32_SIGN_MASK;
    hnum ^= signNum;    // |num|
    hden &= MATS_F32_ABS_MASK; // |den|

    /* purge off exception values */
    if (MATS_F32_UWORD_IS_ZERO(hden) ||
        !MATS_F32_UWORD_IS_FINITE(hnum) ||
        MATS_F32_UWORD_IS_NAN(hden))
    {
        return (num * den) / (num * den);
    }

    if (hden <= MATS_F32_UWORD_HALF_MAX) {
        num = modulus32(num, den + den);
    }
    if ((hnum - hden) == 0) {
        return u32f32(signNum).f;
    }

    num = absolute32(num);
    den = absolute32(den);

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

//
// NOTE(michiel): 64-bit
//

internal f64
floor64(f64 x)
{
    f64 result;
#if MATS_USE_SSE4
    WideMath m;
    m.md = _mm_floor_sd(_mm_set_sd(x), _mm_set_sd(x));
    result = m.ed[0];
#else
    s64 i0 = MATS_S64_FROM_F64(x);
    u64 ix = i0 & MATS_F64_ABS_MASK;
    s32 j0 = (ix >> MATS_F64_EXP_SHIFT) - MATS_F64_EXP_BIAS;

    if (j0 < MATS_F64_EXP_SHIFT)
    {
        if (j0 < 0)
        {
            if (i0 >= 0) {
                i0 = 0;
            } else if (!MATS_F64_UWORD_IS_ZERO(ix)) {
                i0 = 0xBFF0000000000000LL;
            }
        }
        else
        {
            u64 i = MATS_F64_MANT_MASK >> j0;
            if ((i0 & i) == 0) {
                return x;
            }
            if (i0 < 0) {
                i0 += 0x0010000000000000LL >> j0;
            }
            i0 &= ~i;
        }
        result = MATS_F64_FROM_S64(i0);
    }
    else
    {
        if (!MATS_F64_UWORD_IS_FINITE(ix)) {
            result = x + x;
        } else {
            result = x;
        }
    }
#endif
    return result;
}

internal f64
ceil64(f64 x)
{
    f64 result;

#if MATS_USE_SSE4
    WideMath m;
    m.md = _mm_ceil_sd(_mm_set_sd(x), _mm_set_sd(x));
    result = m.ed[0];
#else
    s64 i0 = MATS_S64_FROM_F64(x);
    u64 ix = i0 & MATS_F64_ABS_MASK;
    s32 j0 = (ix >> MATS_F64_EXP_SHIFT) - MATS_F64_EXP_BIAS;

    if (j0 < MATS_F64_EXP_SHIFT)
    {
        if (j0 < 0)
        {
            if (i0 < 0) {
                i0 = MATS_F64_SIGN_MASK;
            } else if (!MATS_F64_UWORD_IS_ZERO(ix)) {
                i0 = 0x3FF0000000000000LL;
            }
        }
        else
        {
            u64 i = MATS_F64_MANT_MASK >> j0;
            if ((i0 & i) == 0) {
                return x;
            }
            if (i0 > 0) {
                i0 += 0x0010000000000000LL >> j0;
            }
            i0 &= ~i;
        }
        result = MATS_F64_FROM_S64(i0);
    }
    else
    {
        if (!MATS_F64_UWORD_IS_FINITE(ix)) {
            result = x + x;
        } else {
            result = x;
        }
    }
#endif
    return result;
}

internal f64
round64(f64 x)
{
    f64 result = 0.0;
#if MATS_USE_SSE4
    // NOTE(michiel): Does round-to-even to break ties (0.5f exact), so it will differ from the 'standard' c/c++ roundf
    WideMath m;
    m.md = _mm_round_sd(_mm_set_sd(x), _mm_set_sd(x), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
    result = m.ed[0];
#else
    u64 w = MATS_U64_FROM_F64(x);

    s32 j0 = (s32)((w & MATS_F64_EXP_MASK) >> MATS_F64_EXP_SHIFT) - MATS_F64_EXP_BIAS;
    if (j0 < MATS_F64_EXP_SHIFT)
    {
        if (j0 < 0) {
            w &= MATS_F32_SIGN_MASK;
            if (j0 == -1) {
                w |= (((u64)MATS_F64_EXP_BIAS) << MATS_F64_EXP_SHIFT);
            }
        }
        else
        {
            u64 exponentMask = MATS_F64_MANT_MASK >> j0;
            if ((w & exponentMask) == 0) {
                return x;
            }
            w += 0x0008000000000000ULL >> j0;
            w &= ~exponentMask;
        }
        result = MATS_F64_FROM_U64(w);
    }
    else
    {
        if (!MATS_F64_UWORD_IS_FINITE(w & MATS_F64_EXP_MASK)) {
            result = x + x;
        } else {
            result = x;
        }
    }
#endif
    return result;
}

internal f64
trunc64(f64 x)
{
    f64 result = 0.0;
#if MATS_USE_SSE4
    WideMath m;
    m.md = _mm_round_sd(_mm_set_sd(x), _mm_set_sd(x), _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC);
    result = m.ed[0];
#else
#if 0
    // TODO(michiel): Or should we just do this?
    result = (f64)(s64)x;
#else
    u64 w = MATS_U64_FROM_F64(x);
    s64 signBit = w & MATS_F64_SIGN_MASK;

    s32 j0 = ((w & MATS_F64_EXP_MASK) >> MATS_F64_EXP_SHIFT) - MATS_F64_EXP_BIAS;
    if (j0 < MATS_F64_EXP_SHIFT)
    {
        if (j0 < 0)
        {
            // NOTE(michiel): -1 < x < 1, so result is +0 or -0
            result = MATS_F64_FROM_S64(signBit);
        }
        else
        {
            result = MATS_F64_FROM_U64((u64)signBit | (w & ~(MATS_F64_MANT_MASK >> j0)));
        }
    }
    else
    {
        if (!MATS_F64_UWORD_IS_FINITE(w & MATS_F64_EXP_MASK)) {
            result = x + x;
        } else {
            result = x;
        }
    }
#endif
#endif
    return result;
}

internal f64
modulus64(f64 num, f64 den)
{
    // NOTE(michiel): Returns 'num % den' where the returned value's sign is equal to the sign of the numerator,
    // and the magnitude less than the magnitude of the denominator.
    // This function computes the remainder from the division of numerator by denominator. Specifically, the return value
    // is numerator - n * denominator, where n is the quotient of numerator divided by denominator, rounded towards zero.
    // Thus, fmod (6.5, 2.3) returns 1.9, which is 6.5 minus 4.6.

    s64 snum = MATS_S64_FROM_F64(num);
    s64 sden = MATS_S64_FROM_F64(den);

    s64 signNum = snum & MATS_F64_SIGN_MASK;

    u64 anum = snum & MATS_F64_ABS_MASK;
    u64 aden = sden & MATS_F64_ABS_MASK;

    // NOTE(michiel): Exception values
    if ((aden == 0) || (aden > 0x7FF0000000000000ULL) ||
        (anum >= 0x7FF0000000000000ULL))
    {
        return (num * den) / (num * den);
    }
    if (anum < aden) {
        return num;
    }
    if (anum == aden) {
        return MATS_F64_FROM_S64(signNum);
    }

    s32 expNum;
    if (anum < 0x0010000000000000ULL)
    {
        expNum = -MATS_F64_EXP_MIN;
        for (s64 i = anum << 11; i > 0; i <<= 1) {
            --expNum;
        }
    }
    else
    {
        expNum = (s32)(anum >> MATS_F64_EXP_SHIFT) - MATS_F64_EXP_BIAS;
    }

    s32 expDen;
    if (aden < 0x0010000000000000ULL)
    {
        expDen = -MATS_F64_EXP_MIN;
        for (s64 i = aden << 11; i > 0; i <<= 1) {
            --expDen;
        }
    }
    else
    {
        expDen = (s32)(aden >> MATS_F64_EXP_SHIFT) - MATS_F64_EXP_BIAS;
    }

    if (expNum >= MATS_F64_EXP_MIN)
    {
        anum = 0x0010000000000000ULL | (anum & MATS_F64_MANT_MASK);
    }
    else
    {
        s32 n = MATS_F64_EXP_MIN - expNum;
        anum = anum << n;
    }

    if (expDen >= MATS_F64_EXP_MIN)
    {
        aden = 0x0010000000000000ULL | (aden & MATS_F64_MANT_MASK);
    }
    else
    {
        s32 n = MATS_F64_EXP_MIN - expDen;
        aden = aden << n;
    }

    s32 n = expNum - expDen;
    while (n--)
    {
        s64 adiff = anum - aden;
        if (adiff < 0) {
            anum = anum + anum;
        } else {
            if (adiff == 0) {
                return MATS_F64_FROM_S64(signNum);
            }
            anum = adiff + adiff;
        }
    }

    s64 adiff = anum - aden;
    if (adiff >= 0) {
        anum = adiff;
    }

    if (anum == 0) {
        return MATS_F64_FROM_S64(signNum);
    }

    while (anum < 0x0010000000000000ULL) {
        anum = anum + anum;
        --expDen;
    }

    f64 result;
    if (expDen >= MATS_F64_EXP_MIN)
    {
        anum = ((anum - 0x0010000000000000ULL) | (((s64)expDen + MATS_F64_EXP_BIAS) << MATS_F64_EXP_SHIFT));
        result = MATS_F64_FROM_U64(signNum | anum);
    }
    else
    {
        s32 n = MATS_F64_EXP_MIN - expDen;
        anum >>= n;
        result = MATS_F64_FROM_U64(signNum | anum);
    }

	return result;
}

internal f64
remainder64(f64 num, f64 den)
{
    // NOTE(michiel): Returns 'num % den' where the absolute value of the result is less than or equal to half the absolute value
    // of the denominator. returned value's sign is equal to the sign of the numerator,
    // This function is like modulus64 except that it rounds the internal quotient n to the nearest integer instead of towards
    // zero. For example, remainder (6.5, 2.3) returns -0.4, which is 6.5 minus 6.9.
    s64 snum = MATS_S64_FROM_F64(num);
    s64 sden = MATS_S64_FROM_F64(den);

    s64 signNum = snum & MATS_F64_SIGN_MASK;

    u64 anum = snum & MATS_F64_ABS_MASK;
    u64 aden = sden & MATS_F64_ABS_MASK;

    if ((aden == 0) || (aden > 0x7FF0000000000000ULL) ||
        (anum >= 0x7FF0000000000000ULL))
    {
        return (num * den) / (num * den);
    }

    if (aden <= MATS_F64_UWORD_HALF_MAX) {
        num = modulus64(num, den + den);
    }
    if (anum == aden) {
        return MATS_F64_FROM_S64(signNum);
    }

    num = absolute64(num);
    den = absolute64(den);

    if (aden < 0x0020000000000000ULL)
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
        f64 denHalf = 0.5 * den;
        if (num > denHalf) {
            num -= den;
            if (num >= denHalf) {
                num -= den;
            }
        }
    }

    f64 result = MATS_F64_FROM_U64(signNum ^ MATS_U64_FROM_F64(num));
    return result;
}
