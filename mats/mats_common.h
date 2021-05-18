// TODO(michiel): Move out all constants

union mats_f32u
{
    f32 f;
    u32 u;
    s32 s;
};

#define MATS_U32_FROM_F32(x)       ((mats_f32u){.f = (x)}).u
#define MATS_S32_FROM_F32(x)       ((mats_f32u){.f = (x)}).s
#define MATS_F32_FROM_U32(x)       ((mats_f32u){.u = (x)}).f
#define MATS_F32_FROM_S32(x)       ((mats_f32u){.s = (x)}).f

#define FLT_UWORD_IS_FINITE(x)     ((x) <  MATS_F32_EXP_MASK)
#define FLT_UWORD_IS_NAN(x)        ((x) >  MATS_F32_EXP_MASK)
#define FLT_UWORD_IS_INFINITE(x)   ((x) == MATS_F32_EXP_MASK)
#define FLT_UWORD_LOG_MAX          0x42B17217
#define FLT_UWORD_LOG_2MAX         0x42B2D4fC
#define FLT_UWORD_HALF_MAX         (MATS_F32_ABS_MASK - (1L << MATS_F32_EXP_SHIFT))
#define FLT_LARGEST_EXP            (MATS_F32_ABS_MASK >> MATS_F32_EXP_SHIFT)
#define FLT_SMALLEST_EXP           -(MATS_F32_EXP_SHIFT - 1)

#define FLT_UWORD_IS_ZERO(x)       ((x) == 0)
#define FLT_UWORD_IS_SUBNORMAL(x)  ((x) <  0x00800000L)
#define FLT_UWORD_LOG_MIN          0x42CFF1B5

struct SinCos32
{
    f32 cos; // NOTE(michiel): cos first, because it is often used as x-axis
    f32 sin; // NOTE(michiel): sin after, because it is often used as y-axis
};

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

/* Top 12 bits of the float representation with the sign bit cleared.  */
internal u32
abstop12_(f32 x)
{
    return MATS_U32_FROM_F32(x) & 0x7FF00000;
}

internal b32
issignaling32(f32 x)
{
    u32 xu = MATS_U32_FROM_F32(x);
    if (!IEEE_754_2008_SNAN)
        return (xu & 0x7FC00000) == 0x7FC00000;
    return 2 * (xu ^ 0x00400000) > 0xFF800000u;
}

//
// NOTE(michiel): Public
//

internal f32
absolute32(f32 x)
{
    return MATS_F32_FROM_U32(MATS_U32_FROM_F32(x) & MATS_F32_ABS_MASK);
}
