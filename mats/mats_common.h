
#define FLT_UWORD_IS_FINITE(x)     ((x) <  0x7f800000L)
#define FLT_UWORD_IS_NAN(x)        ((x) >  0x7f800000L)
#define FLT_UWORD_IS_INFINITE(x)   ((x) == 0x7f800000L)
#define FLT_UWORD_MAX              0x7fffffff
#define FLT_UWORD_LOG_MAX          0x42b17217
#define FLT_UWORD_HALF_MAX         (FLT_UWORD_MAX - (1L << 23))
#define FLT_LARGEST_EXP            (FLT_UWORD_MAX >> 23)
#define FLT_SMALLEST_EXP           -22

#define FLT_UWORD_IS_ZERO(x)       ((x) == 0)
#define FLT_UWORD_IS_SUBNORMAL(x)  ((x) <  0x00800000L)
#define FLT_UWORD_LOG_MIN          0x42cff1b5

internal f32
mats_xflow32(u32 sign, f32 base)
{
    // TODO(michiel): Could add errno
    f32 result = (sign ? -base : base) * base;
    return result;
}

internal f32
mats_underflow32(u32 sign)
{
    return mats_xflow32(sign, 0x1p-95f);
}

internal f32
mats_overflow32(u32 sign)
{
    return mats_xflow32(sign, 0x1p97f);
}

internal f32
mats_divzero32(u32 sign)
{
    f32 result = (sign ? -1.0f : 1.0f) / 0.0f;
    return result;
}

internal f32
mats_invalid32(f32 x)
{
    f32 result = (x - x) / (x - x);
    return result;
}

// TODO(michiel): Rename to absolute32
internal f32
mats_fabs32(f32 x)
{
    return u32f32(u32f32(x).u & 0x7FFFFFFF).f;
}

/* Top 12 bits of the float representation with the sign bit cleared.  */
internal u32
abstop12_(f32 x)
{
    return u32f32(x).u & 0x7ff00000;
}

internal b32
issignaling32(f32 x)
{
    u32 xu = u32f32(x).u;
    if (!IEEE_754_2008_SNAN)
        return (xu & 0x7fc00000) == 0x7fc00000;
    return 2 * (xu ^ 0x00400000) > 0xFF800000u;
}
