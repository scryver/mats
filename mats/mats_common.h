// TODO(michiel): Move out all constants

union mats_f32u
{
    f32 f;
    u32 u;
    s32 s;
};

union mats_f64u
{
    f64 f;
    u64 u;
    s64 s;
    u32 u32s[2];
    s32 s32s[2];
};

#define MATS_U32_FROM_F32(x)       ((mats_f32u){.f = (x)}).u
#define MATS_S32_FROM_F32(x)       ((mats_f32u){.f = (x)}).s
#define MATS_F32_FROM_U32(x)       ((mats_f32u){.u = (x)}).f
#define MATS_F32_FROM_S32(x)       ((mats_f32u){.s = (x)}).f

#define MATS_U64_FROM_F64(x)       ((mats_f64u){.f = (x)}).u
#define MATS_S64_FROM_F64(x)       ((mats_f64u){.f = (x)}).s
#define MATS_F64_FROM_U64(x)       ((mats_f64u){.u = (x)}).f
#define MATS_F64_FROM_S64(x)       ((mats_f64u){.s = (x)}).f

#define MATS_F32_UWORD_IS_FINITE(x)     ((x) <  MATS_F32_EXP_MASK)
#define MATS_F32_UWORD_IS_NAN(x)        ((x) >  MATS_F32_EXP_MASK)
#define MATS_F32_UWORD_IS_INFINITE(x)   ((x) == MATS_F32_EXP_MASK)
#define MATS_F32_UWORD_LOG_MAX          0x42B17217
#define MATS_F32_UWORD_LOG_2MAX         0x42B2D4fC
#define MATS_F32_UWORD_HALF_MAX         (MATS_F32_ABS_MASK - (1L << MATS_F32_EXP_SHIFT))
#define MATS_F32_LARGEST_EXP            (MATS_F32_ABS_MASK >> MATS_F32_EXP_SHIFT)
#define MATS_F32_SMALLEST_EXP           -(MATS_F32_EXP_SHIFT - 1)

#define MATS_F32_UWORD_IS_ZERO(x)       ((x) == 0)
#define MATS_F32_UWORD_IS_SUBNORMAL(x)  ((x) <  0x00800000L)
#define MATS_F32_UWORD_LOG_MIN          0x42CFF1B5

#define MATS_F64_UWORD_IS_FINITE(x)     ((x) <  MATS_F64_EXP_MASK)
#define MATS_F64_UWORD_IS_NAN(x)        ((x) >  MATS_F64_EXP_MASK)
#define MATS_F64_UWORD_IS_INFINITE(x)   ((x) == MATS_F64_EXP_MASK)
//#define MATS_F64_UWORD_LOG_MAX          0x42B17217
//#define MATS_F64_UWORD_LOG_2MAX         0x42B2D4fC
#define MATS_F64_UWORD_HALF_MAX         (MATS_F64_ABS_MASK - (1L << MATS_F64_EXP_SHIFT))
#define MATS_F64_LARGEST_EXP            (MATS_F64_ABS_MASK >> MATS_F64_EXP_SHIFT)
#define MATS_F64_SMALLEST_EXP           -(MATS_F64_EXP_SHIFT - 1)

#define MATS_F64_UWORD_IS_ZERO(x)       ((x) == 0)
#define MATS_F64_UWORD_IS_SUBNORMAL(x)  ((x) <  0x0080000000000000LL)
//#define MATS_F64_UWORD_LOG_MIN          0x42CFF1B5

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

internal f64
mats_xflow64(u64 sign, f64 base)
{
    // TODO(michiel): Could add errno
    f64 result = (sign ? -base : base) * base;
    return result;
}

internal f64
mats_underflow64(u64 sign)
{
    return mats_xflow64(sign, 0x1p-767);
}

internal f64
mats_overflow64(u64 sign)
{
    return mats_xflow64(sign, 0x1p769);
}

internal f64
mats_divzero64(u64 sign)
{
    f64 result = (sign ? -1.0 : 1.0) / 0.0;
    return result;
}

internal f64
mats_invalid64(f64 x)
{
    f64 result = (x - x) / (x - x);
    return result;
}

internal b32
issignaling64(f64 x)
{
    u64 xu = MATS_U64_FROM_F64(x);
    if (!IEEE_754_2008_SNAN)
        return (xu & 0x7FF8000000000000) == 0x7FF8000000000000;
    return 2 * (xu ^ 0x0008000000000000) > 2 * 0x7FF8000000000000ULL;
}

internal u32
mats_top12(f64 x)
{
    return MATS_U64_FROM_F64(x) >> 52;
}

internal u32
mats_top16(f64 x)
{
    return MATS_U64_FROM_F64(x) >> 48;
}

//
// NOTE(michiel): Public
//

internal f32
absolute32(f32 x)
{
    return MATS_F32_FROM_U32(MATS_U32_FROM_F32(x) & MATS_F32_ABS_MASK);
}

internal f64
absolute64(f64 x)
{
    return MATS_F64_FROM_U64(MATS_U64_FROM_F64(x) & MATS_F64_ABS_MASK);
}
