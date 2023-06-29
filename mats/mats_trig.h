
func f64
reduce_fast_pi4(f64 x, int *np)
{
    /* Use scaled float to int conversion with explicit rounding.
       hpi_inv is prescaled by 2^24 so the quadrant ends up in bits 24..31.
       This avoids inaccuracies introduced by truncating negative values.  */
    f64 r = x * 0x1.45F306DC9C883p+23;
    s32 n = ((s32)r + 0x800000) >> 24;
    *np = n;
    return x - n * 0x1.921FB54442D18p0;
}

func f64
sinf_poly_q0(f64 x, f64 x2)
{
    // NOTE(michiel): x - s1*x^3 + s2*x^5 - s3*x^7
    f64 x3 = x * x2;                                // x^3
    f64 s1 = 0x1.1107605230bc4p-7 - x2 * 0x1.994eb3774cf24p-13; // s2 - s3*x^2
    f64 x7 = x3 * x2;                               // x^5
    f64 s = x - x3 * 0x1.555545995a603p-3;                // x - s1*x^3
    return s + x7 * s1;                             // x - s1*x^3 + x^5*(s2 - s3*x^2)
}

func f64
sinf_poly_q1(f64 x2)
{
    // NOTE(michiel): c0 - c1*x^2 + c2*x^4 - c3*x^6 + c4*x^8;
    f64 x4 = x2 * x2;                                 // x^4
    f64 c2 = -0x1.6c087e89a359dp-10 + x2 * 0x1.99343027bf8c3p-16; // -c3 + c4*x^2
    f64 c1 = 1.0f - x2 * 0x1.ffffffd0c621cp-2;              // c0 - c1*x^2
    f64 x6 = x4 * x2;                                 // x^6
    f64 c = c1 + x4 * 0x1.55553e1068f19p-5;                 // c0 - c1*x^2 + c2*x^4
    return c + x6 * c2;                               // c0 - c1*x^2 + c2*x^4 + x^6(-c3 + c4*x^2)
}

func f32
cos32(f32 y)
{
    MATS_ASSERT(absolute32(y) < 120.0f);
    f32 x = y;
    f32 result;
    if (abstop12_(y) < abstop12_(0x1p-12f))
    {
        result = 1.0f;
    }
    else
    {
        int n = 0;
        f64 xm = (abstop12_(y) < abstop12_(0.25f * F32_PI)) ? (f64)x : reduce_fast_pi4((f64)x, &n);
        f64 x2 = xm * xm;
        switch (n & 3)
        {
            default:
            case 0: { result =  (f32)sinf_poly_q1(x2); } break;
            case 1: { result = -(f32)sinf_poly_q0(xm, x2); } break;
            case 2: { result = -(f32)sinf_poly_q1(x2); } break;
            case 3: { result =  (f32)sinf_poly_q0(xm, x2); } break;
        }
    }

    return result;
}

func f32
sin32(f32 y)
{
    MATS_ASSERT(absolute32(y) < 120.0f);
    f32 x = y;

    f32 result;
    if (abstop12_(y) < abstop12_(0x1p-12f))
    {
        result = y;
    }
    else
    {
        int n = 0;
        f64 xm = (abstop12_(y) < abstop12_(0.25f * F32_PI)) ? (f64)x : reduce_fast_pi4((f64)x, &n);
        f64 x2 = xm * xm;
        switch (n & 3)
        {
            default:
            case 0: { result =  (f32)sinf_poly_q0(xm, x2); } break;
            case 1: { result =  (f32)sinf_poly_q1(x2); } break;
            case 2: { result = -(f32)sinf_poly_q0(xm, x2); } break;
            case 3: { result = -(f32)sinf_poly_q1(x2); } break;
        }
    }

    return result;
}

func SinCos32
sincos32(f32 y)
{
    MATS_ASSERT(absolute32(y) < 120.0f);
    f32 x = y;

    SinCos32 result;
    if (abstop12_(y) < abstop12_(0x1p-12f))
    {
        result.cos = 1.0f;
        result.sin = y;
    }
    else
    {
        int n = 0;
        f64 xm = (abstop12_(y) < abstop12_(0.25f * F32_PI)) ? x : reduce_fast_pi4(x, &n);
        f64 x2 = xm * xm;
        f64 sinP = sinf_poly_q0(xm, x2);
        f64 cosP = sinf_poly_q1(x2);
        switch (n & 3)
        {
            case 0: { result.cos =  (f32)cosP; result.sin =  (f32)sinP; } break;
            case 1: { result.cos = -(f32)sinP; result.sin =  (f32)cosP; } break;
            case 2: { result.cos = -(f32)cosP; result.sin = -(f32)sinP; } break;
            case 3: { result.cos =  (f32)sinP; result.sin = -(f32)cosP; } break;
        }
    }

    return result;
}

func f32
tan32_kernel(f32 x, s32 mod)
{
	s32 ix = (s32)MATS_F32U(x).u;
    u32 hx = ix & MATS_F32_ABS_MASK;
    if (hx < 0x31800000)
    {
        if ((s32)x == 0) {
            if ((hx | (mod + 1)) == 0) {
                return 1.0f / MATS_F32U(hx).f;
            } else {
                return (mod == 1) ? x : -1.0f / x;
            }
        }
    }
    f32 z;
    f32 w;
    if (hx >= 0x3F2CA140)
    {
        if (ix < 0) {
            x = -x;
        }
        z = gPiOver4F32_hi - x;
        w = gPiOver4F32_lo - (gPiOver4F32_hi - z - x);
        x = z + w;
    }
    z = x * x;
    w = z * z;
    f32 r = -1.8558637748e-05f * w;
    r = (r + 7.8179444245e-05f) * w;
    r = (r + 5.8804126456e-04f) * w;
    r = (r + 3.5920790397e-03f) * w;
    r = (r + 2.1869488060e-02f) * w;
    r = (r + 1.3333334029e-01f);

    f32 v =  2.5907305826e-05f * w;
    v = (v + 7.1407252108e-05f) * w;
    v = (v + 2.4646313977e-04f) * w;
    v = (v + 1.4562094584e-03f) * w;
    v = (v + 8.8632395491e-03f) * w;
    v = (v + 5.3968254477e-02f) * z;

    f32 s = z * x;
    r = z * (s * (r + v));
    r += 3.3333334327e-01f * s;
    w = x + r;

    if (hx >= 0x3F2CA140)
    {
        v = (f32)mod;
        return (f32)(1 - ((ix >> 30) & 2)) * (v - 2.0f * (x - (w * w / (w + v) - r)));
    }
    if (mod == 1) {
        return w;
    } else {
        z = w;
        u32 uz = MATS_F32U(z).u & 0xFFFFF000;
        z = MATS_F32U(uz).f;
        v = r - (z - x);
        f32 t = -1.0f / w;
        f32 a = t;
        t = MATS_F32U(MATS_F32U(t).u & 0xFFFFF000).f;
        s = 1.0f + t * z;
        return t + a * (s + t * v);
    }
}

func f32
tan32(f32 y)
{
    MATS_ASSERT(absolute32(y) < 120.0f);
    f32 result;
    u32 uy = MATS_F32U(y).u & MATS_F32_ABS_MASK;
    if (uy < 0x3F490FDA)
    {
        result = tan32_kernel(y, 1);
    }
    else
    {
        int n = 0;
        f32 x = (f32)reduce_fast_pi4(y, &n);
        result = tan32_kernel(x, 1 - ((n & 1) << 1));
    }
    return result;
}

// NOTE(michiel): The smaller poly were found in the musl c library version of the sun93 code
#define MATS_ASINCOS_USE_SMALL_POLY 1
#define MATS_ATAN_USE_SMALL_POLY    1

#if MATS_ASINCOS_USE_SMALL_POLY
#define pS0   1.6666586697e-01f
#define pS1  -4.2743422091e-02f
#define pS2  -8.6563630030e-03f
#define qS1  -7.0662963390e-01f
#define poly(z) \
f32 p = pS2 * z; \
p = (p + pS1) * z; \
p = (p + pS0) * z; \
f32 q = qS1 * z; \
q = (q + 1.0f); \
f32 r = p / q

#else
#define pS0   1.6666667163e-01f
#define pS1  -3.2556581497e-01f
#define pS2   2.0121252537e-01f
#define pS3  -4.0055535734e-02f
#define pS4   7.9153501429e-04f
#define pS5   3.4793309169e-05f
#define qS1  -2.4033949375e+00f
#define qS2   2.0209457874e+00f
#define qS3  -6.8828397989e-01f
#define qS4   7.7038154006e-02f
#define poly(z) \
f32 p = pS5 * z; \
f32 q = qS4 * z; \
p = (p + pS4) * z; \
q = (q + qS3) * z; \
p = (p + pS3) * z; \
q = (q + qS2) * z; \
p = (p + pS2) * z; \
p = (p + pS1) * z; \
q = (q + qS1) * z; \
p = (p + pS0) * z; \
q = (q + 1.0f); \
f32 r = p / q

#endif

func f32
acos32(f32 x)
{
#define pi        3.1415925026e+00f
#define pio2_hi   1.5707962513e+00f
#define pio2_lo   7.5497894159e-08f

	s32 ix = (s32)MATS_F32U(x).u;
    u32 hx = ix & MATS_F32_ABS_MASK;

    if (hx == 0x3F800000) {
        if (ix > 0) {
            return 0.0f;                 // NOTE(michiel): x == 1.0f
        } else {
            return pi + 2.0f * pio2_lo;  // NOTE(michiel): x == -1.0f
        }
    } else if (hx > 0x3F800000) {
        return (x - x) / (x - x);        // NOTE(michiel): |x| > 1.0f
    }

    if (hx < 0x3F000000)
    {
        if (hx <= 0x23000000) {
            return pio2_hi + pio2_lo;    // NOTE(michiel): |x| <= 2^-57
        }
        f32 x2 = x * x;
        poly(x2);
        return pio2_hi - (x - (pio2_lo - x * r)); // NOTE(michiel): |x| < 0.5f
    }
    else if (ix < 0)
    {
        f32 z = (1.0f + x) * 0.5f;
        poly(z);
        f32 s = sqrt32(z);
        f32 w = r * s - pio2_lo;
        return pi - 2.0f * (s + w); // NOTE(michiel): x <= -0.5f
    }
    else
    {
        f32 z = (1.0f - x) * 0.5f;
        f32 s = sqrt32(z);
        f32 df = MATS_F32U(MATS_F32U(s).u & 0xFFFFF000).f;
        f32 c = (z - df * df) / (s + df);
        poly(z);
        f32 w = r * s + c;
        return 2.0f * (df + w); // NOTE(michiel): x >= 0.5f
    }

#undef pio2_hi
#undef pio2_lo
#undef pi
}

func f32
asin32(f32 x)
{
#define pio2_hi 1.57079637050628662109375f
#define pio2_lo -4.37113900018624283e-8f
#define pio4_hi 0.785398185253143310546875f

    s32 ix = (s32)MATS_F32U(x).u;
    u32 hx = ix & MATS_F32_ABS_MASK;

    if (hx == 0x3F800000) {
        return x * pio2_hi + x * pio2_lo;
    } else if (hx > 0x3F800000) {
        return (x - x) / (x - x);
    } else if (hx < 0x3F000000) {
        if (hx < 0x32000000) {
            if (gHugeF32 + x > 1.0f) {
                return x;
            }
        } else {
            f32 x2 = x * x;
            poly(x2);
            return x + x * r;
        }
    }

    f32 w = 1.0f - MATS_F32U(hx).f;
    f32 t = w * 0.5f;
    poly(t);
    f32 s = sqrt32(t);
    if (hx >= 0x3F79999A)
    {
        t = pio2_hi - (2.0f * (s + s * r) - pio2_lo);
    }
    else
    {
        w = s;
        u32 iw = MATS_F32U(w).u;
        w = MATS_F32U(iw & 0xFFFFF000).f;
        f32 c = (t - w * w) / (s + w);
        p = 2.0f * s * r - (pio2_lo - 2.0f * c);
        q = pio4_hi - 2.0f * w;
        t = pio4_hi - (p - q);
    }

    return (ix < 0) ? -t : t;
#undef pio4_hi
#undef pio2_hi
#undef pio2_lo
}

#undef poly
#undef pS0
#undef pS1
#undef pS2
#undef pS3
#undef pS4
#undef pS5
#undef qS1
#undef qS2
#undef qS3
#undef qS4

func f32
atan32(f32 x)
{
    s32 ix = (s32)MATS_F32U(x).u;
    u32 hx = ix & MATS_F32_ABS_MASK;

    if (hx >= 0x50800000)
    {
        // NOTE(michiel): |x| >= 2^34
        if (MATS_F32_UWORD_IS_NAN(hx)) {
            return x + x;
        } else if (ix > 0) {
            return 1.5707962513e+00f + 7.5497894159e-08f;
        } else {
            return -1.5707962513e+00f - 7.5497894159e-08f;
        }
    }
    else
    {
        s32 id = 0;
        if (hx < 0x3EE00000) { // NOTE(michiel): |x| < 0.4375
            if (hx < 0x31000000) { // NOTE(michiel): |x| < 2^-29
                if ((gHugeF32 + x) > 1.0f) {
                    return x; // NOTE(michiel): Raise inexact
                }
            }
            id = -1;
        }
        else
        {
            x = MATS_F32U(hx).f;
            if (hx < 0x3F980000) { // NOTE(michiel): |x| < 1.1875
                if (hx < 0x3F300000) { // NOTE(michiel): 7/16 <= |x| < 11/16
                    id = 0;
                    x = (2.0f * x - 1.0f) / (2.0f + x);
                } else {
                    id = 1;
                    x = (x - 1.0f) / (x + 1.0f);
                }
            } else {
                if (hx < 0x401C0000) {
                    id = 2;
                    x = (x - 1.5f) / (1.0f + 1.5f * x);
                } else {
                    id = 3;
                    x = -1.0f / x;
                }
            }
        }

        f32 x2 = x * x;
#if MATS_ATAN_USE_SMALL_POLY
        f32 x4 = x2 * x2;
        f32 s1 =   6.1687607318e-02f * x4;
        f32 s2 = - 1.0648017377e-01f * x4;
        s1 = (s1 + 1.4253635705e-01f) * x4;
        s2 = (s2 - 1.9999158382e-01f) * x4;
        s1 = (s1 + 3.3333328366e-01f) * x2;
#else
#if MATS_USE_SSE
        WideMath x2w; x2w.m = _mm_set1_ps(x2);
        WideMath x4w; x4w.m = _mm_mul_ps(x2w.m, x2w.m);

        WideMath coef0; coef0.m = _mm_setr_ps(3.3333334327e-01f, 0.0f, 0.0f, 0.0f);
        WideMath coef1; coef1.m = _mm_setr_ps(1.4285714924e-01f, -2.0000000298e-01f, 0, 0);
        WideMath coef2; coef2.m = _mm_setr_ps(9.0908870101e-02f, -1.1111110449e-01f, 0, 0);
        WideMath coef3; coef3.m = _mm_setr_ps(6.6610731184e-02f, -7.6918758452e-02f, 0, 0);
        WideMath coef4; coef4.m = _mm_setr_ps(4.9768779427e-02f, -5.8335702866e-02f, 0, 0);
        WideMath coef5; coef5.m = _mm_setr_ps(1.6285819933e-02f, -3.6531571299e-02f, 0, 0);
        WideMath s1s2; s1s2.m = _mm_mul_ps(coef5.m, x4w.m);
        s1s2.m = _mm_mul_ps(_mm_add_ps(s1s2.m, coef4.m), x4w.m);
        s1s2.m = _mm_mul_ps(_mm_add_ps(s1s2.m, coef3.m), x4w.m);
        s1s2.m = _mm_mul_ps(_mm_add_ps(s1s2.m, coef2.m), x4w.m);
        s1s2.m = _mm_mul_ps(_mm_add_ps(s1s2.m, coef1.m), x4w.m);
        s1s2.m = _mm_add_ps(s1s2.m, coef0.m);

        f32 s1 = s1s2.e[0] * x2;
        f32 s2 = s1s2.e[1];
#else
        f32 x4 = x2 * x2;
        f32 s1 =   1.6285819933e-02f * x4;
        f32 s2 = - 3.6531571299e-02f * x4;
        s1 = (s1 + 4.9768779427e-02f) * x4;
        s2 = (s2 - 5.8335702866e-02f) * x4;
        s1 = (s1 + 6.6610731184e-02f) * x4;
        s2 = (s2 - 7.6918758452e-02f) * x4;
        s1 = (s1 + 9.0908870101e-02f) * x4;
        s2 = (s2 - 1.1111110449e-01f) * x4;
        s1 = (s1 + 1.4285714924e-01f) * x4;
        s2 = (s2 - 2.0000000298e-01f) * x4;
        s1 = (s1 + 3.3333334327e-01f) * x2;
#endif
#endif

        f32 result = x * (s1 + s2);
        if (id < 0) {
            result = x - result;
        } else {
            result = gAtanHiF32[id] - ((result - gAtanLoF32[id]) - x);
            result = (ix < 0) ? -result : result;
        }

        return result;
    }
}

func f32
atan2_32(f32 y, f32 x)
{
    s32 ix = MATS_S32_FROM_F32(x);
    u32 hx = ix & MATS_F32_ABS_MASK;
    s32 iy = MATS_S32_FROM_F32(y);
    u32 hy = iy & MATS_F32_ABS_MASK;

    if (MATS_F32_UWORD_IS_NAN(hx) || MATS_F32_UWORD_IS_NAN(hy)) {
        return x + y;
    }
    if (ix == 0x3F800000) {
        return atan32(y); // NOTE(michiel): x == 1.0f
    }
    s32 m = ((iy >> 31) & 1) | ((ix >> 30) & 2); // NOTE(michiel): 2 * sign(x) + sign(y)
    if (MATS_F32_UWORD_IS_ZERO(hy))
    {
        switch (m)
        {
            case 0:
            case 1: { return y; } break;                  // NOTE(michiel): atan(+/-0,+anything) = +/-0
            case 2: { return  gPiF32 + gTinyF32; } break; // NOTE(michiel): atan(+0,  -anything) =  pi
            case 3: { return -gPiF32 - gTinyF32; } break; // NOTE(michiel): atan(-0,  -anything) = -pi
        }
    }

    if (MATS_F32_UWORD_IS_ZERO(hx))
    {
        return (iy < 0) ? -gPiOver2F32 - gTinyF32 : gPiOver2F32 + gTinyF32;
    }

    if (MATS_F32_UWORD_IS_INFINITE(hx))
    {
        if (MATS_F32_UWORD_IS_INFINITE(hy))
        {
            switch (m)
            {
                case 0: { return  gPiOver4F32 + gTinyF32; } break;        // NOTE(michiel): atan(+inf, +inf) =   pi / 4
                case 1: { return -gPiOver4F32 - gTinyF32; } break;        // NOTE(michiel): atan(-inf, +inf) =  -pi / 4
                case 2: { return  3.0f * gPiOver4F32 + gTinyF32; } break; // NOTE(michiel): atan(+inf, -inf) =  3pi / 4
                case 3: { return -3.0f * gPiOver4F32 - gTinyF32; } break; // NOTE(michiel): atan(-inf, -inf) = -3pi / 4
            }
        }
        else
        {
            switch (m)
            {
                case 0: { return  0.0f; } break;              // NOTE(michiel): atan(+anything,+inf) =  0
                case 1: { return -0.0f; } break;              // NOTE(michiel): atan(-anything,+inf) = -0
                case 2: { return  gPiF32 + gTinyF32; } break; // NOTE(michiel): atan(+anything,-inf) =  pi
                case 3: { return -gPiF32 - gTinyF32; } break; // NOTE(michiel): atan(-anything,-inf) = -pi
            }
        }
    }

    if (MATS_F32_UWORD_IS_INFINITE(hy)) {
        return (iy < 0) ? -gPiOver2F32 - gTinyF32 : gPiOver2F32 + gTinyF32;
    }

    f32 result;
    s32 k = ((s32)hy - (s32)hx) >> MATS_F32_EXP_SHIFT;
    if (k > 60) {
        result = gPiOver2F32 + 0.5f * gPiF32_lo;
    } else if ((ix < 0) && (k < -60)) {
        result = 0.0f;
    } else {
        result = atan32(absolute32(y / x));
    }
    switch (m)
    {
        case 0: {} break;
        case 1: { result = MATS_F32U(MATS_F32U(result).u ^ MATS_F32_SIGN_MASK).f; } break;
        case 2: { result = gPiF32 - (result - gPiF32_lo); } break;
        case 3: { result = (result - gPiF32_lo) - gPiF32; } break;
    }
    return result;
}

//
// NOTE(michiel): Hyperbolic functions
//

func f32
cosh32(f32 x)
{
    s32 ix = (s32)MATS_F32U(x).u;
    ix &= MATS_F32_ABS_MASK;

    if (!MATS_F32_UWORD_IS_FINITE(ix)) {
        return x * x;
    }

    if (ix < 0x3EB17218)
    {
        // NOTE(michiel): |x| in [0, ln(2)/2], return 1 + expm1(|x|)^2 / (2*exp(|x|))
        f32 t = expm1_32(absolute32(x));
        f32 w = 1.0f + t;
        if (ix < 0x24000000) {
            // NOTE(michiel): |x| < 2^-55
            return w;
        } else {
            return 1.0f + (t * t) / (w + w);
        }
    }
    else if (ix < 0x41B00000)
    {
        // NOTE(michiel): |x| in [ln(2)/2, 22], return (exp(|x|) + 1/exp(|x|))/2
        f32 t = exp32(absolute32(x));
        return 0.5f * t + 0.5f / t;
    }
    else if (ix <= MATS_F32_UWORD_LOG_MAX)
    {
        // NOTE(michiel): |x| in [22, log(FLT_MAX)], return exp(|x|)/2
        return 0.5f * exp32(absolute32(x));
    }
    else if (ix <= MATS_F32_UWORD_LOG_2MAX)
    {
        // NOTE(michiel): |x| in[log(FLT_MAX), overflowThreshold], return exp(|x|/2)^2/2
        f32 w = exp32(0.5f*absolute32(x));
        f32 t = 0.5f * w;
        return t * w;
    }
    else
    {
        return mats_overflow32(0);
    }
}

func f32
sinh32(f32 x)
{
    s32 jx = (s32)MATS_F32U(x).u;
    s32 ix = jx & MATS_F32_ABS_MASK;

    if (!MATS_F32_UWORD_IS_FINITE(ix)) {
        return x + x;
    }

    f32 h = (jx < 0) ? -0.5f : 0.5f;

    if (ix < 0x41B00000)
    {
        // NOTE(michiel): |x| in [0, 22], return sign(x) * 0.5 * (E + E/(E + 1))
        if (ix < 0x31800000) {
            // NOTE(michiel): |x| < 2^-28
            return x;
        }
        f32 t = expm1_32(absolute32(x));
        if (ix < 0x3F800000) {
            return h * (2.0f * t - t * t / (t + 1.0f));
        } else {
            return h * (t + t / (t + 1.0f));
        }
    }
    else if (ix <= MATS_F32_UWORD_LOG_MAX)
    {
        // NOTE(michiel): |x| in [22, log(FLT_MAX)], return exp(|x|)/2
        return h * exp32(absolute32(x));
    }
    else if (ix <= MATS_F32_UWORD_LOG_2MAX)
    {
        // NOTE(michiel): |x| in[log(FLT_MAX), overflowThreshold], return exp(|x|/2)^2/2
        f32 w = exp32(0.5f*absolute32(x));
        f32 t = h * w;
        return t * w;
    }
    else
    {
        return x * 1.0e37f;
    }
}

func SinCos32
sinhcosh32(f32 x)
{
    SinCos32 result;
    if (absolute32(x) <= 0.5f)
    {
        result.cos = cosh32(x);
        result.sin = sinh32(x);
    }
    else
    {
        f32 e = exp32(x);
        f32 ei = 0.5f / e;
        e  = 0.5f * e;
        result.cos = e + ei;
        result.sin = e - ei;
    }
    return result;
}

func f32
tanh32(f32 x)
{
    s32 jx = (s32)MATS_F32U(x).u;
    s32 ix = jx & MATS_F32_ABS_MASK;

    if (!MATS_F32_UWORD_IS_FINITE(ix))
    {
        if (jx >= 0) {
            return 1.0f / x + 1.0f; // NOTE(michiel): tanh(+/-inf) = +/-1
        } else {
            return 1.0f / x - 1.0f; // NOTE(michiel): tanh(NaN) = NaN
        }
    }

    f32 result;
    if (ix < 0x41B00000)
    {
        // NOTE(michiel): |x| < 22
        if (ix < 0x24000000) {
            // NOTE(michiel): |x| < 2^-55
            return x * (1.0f + x); // NOTE(michiel): tanh(small) = small
        }

        if (ix >= 0x3F800000)
        {
            // NOTE(michiel): |x| >= 1
            f32 t = expm1_32(2.0f * absolute32(x));
            result = 1.0f - 2.0f / (t + 2.0f);
        }
        else
        {
            f32 t = expm1_32(-2.0f * absolute32(x));
            result = -t / (t + 2.0f);
        }
    }
    else
    {
        // NOTE(michiel): |x| >= 22, return +/-1
        result = 1.0f - gTinyF32;
    }
    return jx < 0 ? -result : result;
}

func f32
acosh32(f32 x)
{
    s32 jx = (s32)MATS_F32U(x).u;

    if (jx < 0x3F800000) {
        // NOTE(michiel): x < 1
        return (x - x) / (x - x);
    } else if (jx >= 0x4D800000) {
        // NOTE(michiel): x >= 2^28
        if (!MATS_F32_UWORD_IS_FINITE(jx)) {
            // NOTE(michiel): Inf or NaN
            return x + x;
        } else {
            return log32(x) + gLn2F32;
        }
    } else if (jx == 0x3F800000) {
        // NOTE(michiel): x == 1
        return 0.0f;
    } else if (jx > 0x40000000) {
        // NOTE(michiel): 2 < x < 2^28
        f32 t = x * x;
        return log32(2.0f * x - 1.0f / (x + sqrt32(t - 1.0f)));
    } else {
        // NOTE(michiel): 1 < x <= 2
        f32 t = x - 1.0f;
        return log1p32(t + sqrt32(2.0f * t + t * t));
    }
}

func f32
asinh32(f32 x)
{
    s32 jx = (s32)MATS_F32U(x).u;
    s32 ix = jx & MATS_F32_ABS_MASK;

    f32 result;
    if (!MATS_F32_UWORD_IS_FINITE(ix)) {
        // NOTE(michiel): Inf or NaN
        return x + x;
    } else if (ix < 0x31800000) {
        // NOTE(michiel): |x| < 2^-28
        return x;
    } else if (ix >= 0x4D800000) {
        // NOTE(michiel): x > 2^28
        result = log32(absolute32(x)) + gLn2F32;
    } else if (ix > 0x40000000) {
        // NOTE(michiel): 2 < x < 2^28
        f32 t = absolute32(x);
        result = log32(2.0f * t + 1.0f / (sqrt32(x * x + 1.0f) + t));
    } else {
        // NOTE(michiel): 2^-28 <= |x| <= 2
        f32 t = x * x;
        result = log1p32(absolute32(x) + t / (1.0f + sqrt32(1.0f + t)));
    }

    return jx > 0 ? result : -result;
}

func f32
atanh32(f32 x)
{
    s32 jx = (s32)MATS_F32U(x).u;
    s32 ix = jx & MATS_F32_ABS_MASK;

    if (ix > 0x3F800000) {
        // NOTE(michiel): |x| > 1
        return (x - x) / (x - x);
    } else if (ix == 0x3F800000) {
        // NOTE(michiel): |x| == 1
        return x / 0.0f;
    } else if (ix < 0x31800000) {
        // NOTE(michiel): |x| < 2^-28
        return x;
    } else {
        x = MATS_F32U((u32)ix).f;
        f32 result;
        if (ix < 0x3F000000) {
            // NOTE(michiel): |x| < 0.5
            f32 t = x + x;
            result = 0.5f * log1p32(t + t * x / (1.0f - x));
        } else {
            result = 0.5f * log1p32((x + x) / (1.0f - x));
        }
        return jx >= 0 ? result : -result;
    }
}

//
// NOTE(michiel): 64 bit
//

global const s32 gRemPiInitJK[] = {
    2, 3, 4, 6,
};

func s32
kernel_rem_pi_over_2(f64 *x, f64 *y, s32 e0, s32 nx, s32 prec, const s32 *twoOverPi)
{
    s32 jk = gRemPiInitJK[prec];
    s32 jp = jk;

    s32 jx = nx - 1;
    s32 jv = (e0 - 3) / 24;
    if (jv < 0) {
        jv = 0;
    }
    s32 q0 = e0 - 24 * (jv + 1);

    s32 j = jv - jx;
    s32 m = jx + jk;
    f64 f[20];
    for (s32 i = 0; i < m; ++i, ++j)
    {
        f[i] = (j < 0) ? 0.0 : (f64)twoOverPi[j];
    }

    f64 q[20];
    for (s32 i = 0; i < jk; ++i)
    {
        f64 fw = 0.0;
        for (j = 0; j <= jx; ++j) {
            fw += x[j] * f[jx + i - j];
        }
        q[i] = fw;
    }

    s32 jz = jk;

    s32 n,iq[20],i,k,ih;
    f64 z,fq[20];

    recompute:
    z = q[jz];
    i = 0;
    for (j = jz; j > 0; --j, ++i)
    {
        f64 fw = (f64)((s32)(g2powMin24F64 * z));
        iq[i] = (s32)(z - g2pow24F64 * fw);
        z = q[j - 1] + fw;
    }

    z = scalbn64(z, q0);
    z -= 8.0 * floor64(z * 0.125);
    n = (s32)z;
    z -= (f64)n;

    ih = 0;
    if (q0 > 0)
    {
        i = (iq[jz - 1] >> (24 - q0));
        n += i;
        iq[jz - 1] -= i << (24 - q0);
        ih = iq[jz - 1] >> (23 - q0);
    }
    else if (q0 == 0)
    {
        ih = iq[jz - 1] >> 23;
    }
    else if (z >= 0.5)
    {
        ih = 2;
    }

    if (ih > 0)
    {
	    ++n;
        s32 carry = 0;
	    for (i = 0; i < jz; ++i)
        {
            j = iq[i];
            if (carry == 0)
            {
                if (j != 0) {
                    carry = 1;
                    iq[i] = 0x01000000 - j;
                }
            }
            else
            {
                iq[i] = 0x00FFFFFF - j;
            }
	    }
	    if (q0 > 0)
        {
            /* rare case: chance is 1 in 12 */
	        switch(q0)
            {
                case 1: { iq[jz-1] &= 0x007FFFFF; } break;
                case 2: { iq[jz-1] &= 0x003FFFFF; } break;
                default: {} break;
	        }
	    }
	    if (ih == 2)
        {
            z = 1.0 - z;
            if (carry != 0) {
                z -= scalbn64(1.0, q0);
            }
	    }
	}

    /* check if recomputation is needed */
	if (z == 0.0)
    {
	    j = 0;
	    for (i = jz - 1; i >= jk; --i) {
            j |= iq[i];
        }
	    if (j == 0)
        {
            /* need recomputation */
            for (k = 1; iq[jk - k] == 0; ++k) {}   /* k = no. of terms needed */

            for (i = jz + 1; i <= jz + k; ++i)
            {
                /* add q[jz+1] to q[jz+k] */
                f[jx + i] = (f64)twoOverPi[jv + i];
                f64 fw = 0.0;
                for (j = 0; j <= jx; ++j) {
                    fw += x[j]*f[jx + i - j];
                }
                q[i] = fw;
            }
            jz += k;
            goto recompute;
	    }
	}

    /* chop off zero terms */
	if (z == 0.0)
    {
	    --jz;
        q0 -= 24;
	    while (iq[jz] == 0) {
            jz--;
            q0-=24;
        }
	}
    else
    {
        /* break z into 24-bit if necessary */
	    z = scalbn64(z, -q0);
	    if (z >= g2pow24F64)
        {
            f64 fw = (f64)((s32)(g2powMin24F64 * z));
            iq[jz] = (s32)(z - g2pow24F64 * fw);
            ++jz;
            q0 += 24;
            iq[jz] = (s32)fw;
	    }
        else
        {
            iq[jz] = (s32)z;
        }
	}

    /* convert integer "bit" chunk to floating-point value */
	f64 fw = scalbn64(1.0, q0);
	for (i = jz; i >= 0; --i) {
	    q[i] = fw * (f64)iq[i];
        fw *= g2powMin24F64;
	}

    /* compute PIo2[0,...,jp]*q[jz,...,0] */
	for (i = jz; i >= 0; --i) {
        fw = 0.0;
	    for (k = 0; (k <= jp) && (k <= (jz - i)); ++k) {
            fw += gPiOver2F64_Table[k] * q[i + k];
        }
	    fq[jz - i] = fw;
	}

    /* compress fq[] into y[] */
	switch(prec)
    {
        default:
	    case 0: {
            fw = 0.0;
            for (i = jz; i >= 0; --i) {
                fw += fq[i];
            }
            y[0] = (ih == 0) ? fw : -fw;
        } break;

	    case 1:
	    case 2: {
            fw = 0.0;
            for (i = jz; i >= 0; --i) {
                fw += fq[i];
            }
            y[0] = (ih == 0) ? fw : -fw;
            fw = fq[0] - fw;

            for (i = 1; i <= jz; ++i) {
                fw += fq[i];
            }
            y[1] = (ih == 0) ? fw: -fw;
		} break;

	    case 3: {
            /* painful */
            for (i = jz; i > 0; --i)
            {
                fw = fq[i - 1] + fq[i];
                fq[i] += fq[i - 1] - fw;
                fq[i - 1] = fw;
            }
            for (i = jz; i > 1; --i)
            {
                fw = fq[i - 1] + fq[i];
                fq[i] += fq[i - 1] - fw;
                fq[i - 1] = fw;
            }
            fw = 0.0;
            for (i = jz; i >= 2; --i) {
                fw += fq[i];
            }

            if (ih == 0) {
                y[0] =  fq[0];
                y[1] =  fq[1];
                y[2] =  fw;
            } else {
                y[0] = -fq[0];
                y[1] = -fq[1];
                y[2] = -fw;
            }
        } break;
    }

    return n & 7;
}

func f64
kernel_cos64(f64 x, f64 y)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;
    if(top < 0x3E400000)
    {
        /* if x < 2**27 */
        if(((int)x)==0) {
            return 1.0;         /* generate inexact */
        }
    }

    f64 z = x * x;
#define C1   4.16666666666666019037e-02 /* 0x3FA55555, 0x5555554C */
#define C2  -1.38888888888741095749e-03 /* 0xBF56C16C, 0x16C15177 */
#define C3   2.48015872894767294178e-05 /* 0x3EFA01A0, 0x19CB1590 */
#define C4  -2.75573143513906633035e-07 /* 0xBE927E4F, 0x809C52AD */
#define C5   2.08757232129817482790e-09 /* 0x3E21EE9E, 0xBDB4B1C4 */
#define C6  -1.13596475577881948265e-11 /* 0xBDA8FAE9, 0xBE8838D4 */

    f64 z2 = z * z;
    f64 p1 = z2 * (C2 + z2 * (C4 + z2 * C6));
    f64 p2 = z * (C1 + z2 * (C3 + z2 * C5));
    f64 r = p1 + p2;

#undef C6
#undef C5
#undef C4
#undef C3
#undef C2
#undef C1

    f64 result;
    if (top < 0x3FD33333)
    {
        result = 1.0 - (0.5 * z - (z * r - x * y));
    }
    else
    {
        f64 qx;
        if (top > 0x3FE90000)
        {
            // NOTE(michiel): x > 0.78125
            qx = 0.28125;
        }
        else
        {
            // NOTE(michiel): x / 4
            qx = MATS_F64_FROM_U64((u64)(top - 0x00200000) << 32);
        }
        f64 hz = 0.5 * z - qx;
        f64 a = 1.0 - qx;
        result = a - (hz - (z * r - x * y));
    }
    return result;
}

func f64
kernel_sin64(f64 x, f64 y, s32 iy)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;
    if(top < 0x3E400000)
    {
        /* if x < 2**27 */
        if(((int)x)==0) {
            return x;         /* generate inexact */
        }
    }

    f64 z = x * x;
#define S1  -1.66666666666666324348e-01 /* 0xBFC55555, 0x55555549 */
#define S2   8.33333333332248946124e-03 /* 0x3F811111, 0x1110F8A6 */
#define S3  -1.98412698298579493134e-04 /* 0xBF2A01A0, 0x19C161D5 */
#define S4   2.75573137070700676789e-06 /* 0x3EC71DE3, 0x57B1FE7D */
#define S5  -2.50507602534068634195e-08 /* 0xBE5AE5E6, 0x8A2B9CEB */
#define S6   1.58969099521155010221e-10 /* 0x3DE5D93A, 0x5ACFD57C */

    f64 z2 = z * z;
    f64 p1 = S2 + z2 * (S4 + z2 * S6);
    f64 p2 = z * (S3 + z2 * S5);
    f64 r = p1 + p2;

#undef S6
#undef S5
#undef S4
#undef S3
#undef S2
    f64 v = x * z;
    f64 result;
    if (iy == 0)
    {
        result = x + v * (S1 + z * r);
    }
    else
    {
        result = x - ((z * (0.5 * y - v * r) - y) - v * S1);
    }
#undef S1

    return result;
}

func f64
kernel_tan64(f64 x, f64 y, s32 iy)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;

    if (top < 0x3E300000)
    {
        if ((s32)x == 0)
        {
            u32 low = ix & 0xFFFFFFFF;
            if (((top | low) | (iy + 1)) == 0)
            {
                return 1.0 / absolute64(x);
            }
            else
            {
                if (iy == 1)
                {
                    return x;
                }
                else
                {
                    f64 w = x + y;
                    f64 z = MATS_F64_FROM_U64(MATS_U64_FROM_F64(w) & 0xFFFFFFFF00000000ULL);
                    f64 v = y - (z - x);
                    f64 a = -1.0 / w;
                    f64 t = MATS_F64_FROM_U64(MATS_U64_FROM_F64(a) & 0xFFFFFFFF00000000ULL);
                    f64 s = 1.0 + t * z;
                    return t + a * (s + t * v);
                }
            }
        }
    }

    if (top >= 0x3FE59428)
    {
        if (ix < 0) {
            x = -x;
            y = -y;
        }
        f64 z = gPiOver4F64 - x;
        f64 w = gPiOver4F64_lo - y;
        x = z + w;
        y = 0.0;
    }

    f64 z = x * x;
    f64 w = z * z;

#define T1   3.33333333333334091986e-01 /* 0x3FD5555555555563 */
#define T2   1.33333333333201242699e-01 /* 0x3FC111111110FE7A */
#define T3   5.39682539762260521377e-02 /* 0x3FABA1BA1BB341FE */
#define T4   2.18694882948595424599e-02 /* 0x3F9664F48406D637 */
#define T5   8.86323982359930005737e-03 /* 0x3F8226E3E96E8493 */
#define T6   3.59207910759131235356e-03 /* 0x3F6D6D22C9560328 */
#define T7   1.45620945432529025516e-03 /* 0x3F57DBC8FEE08315 */
#define T8   5.88041240820264096874e-04 /* 0x3F4344D8F2F26501 */
#define T9   2.46463134818469906812e-04 /* 0x3F3026F71A8D1068 */
#define T10  7.81794442939557092300e-05 /* 0x3F147E88A03792A6 */
#define T11  7.14072491382608190305e-05 /* 0x3F12B80F32F0A7E9 */
#define T12 -1.85586374855275456654e-05 /* 0xBEF375CBDB605373 */
#define T13  2.59073051863633712884e-05 /* 0x3EFB2A7074BF7AD4 */

    f64 r = T2 + w * (T4 + w * (T6 + w * (T8 + w * (T10 + w * T12))));
    f64 v = z * (T3 + w * (T5 + w * (T7 + w * (T9 + w * (T11 + w * T13)))));
    f64 s = z * x;
    r = y + z * (s * (r + v) + y);
    r += T1 * s;
    w = x + r;

#undef T13
#undef T12
#undef T11
#undef T10
#undef T9
#undef T8
#undef T7
#undef T6
#undef T5
#undef T4
#undef T3
#undef T2
#undef T1

    if (top >= 0x3FE59428)
    {
        v = (f64)iy;
        return (f64)(1 - ((ix >> 62) & 0x2)) * (v - 2.0 * (x - (w * w / (w + v) - r)));
    }
    if (iy == 1)
    {
        return w;
    }
    else
    {
        z = MATS_F64_FROM_U64(MATS_U64_FROM_F64(w) & 0xFFFFFFFF00000000ULL);
        v = r - (z - x);
        f64 a = -1.0 / w;
        f64 t = MATS_F64_FROM_U64(MATS_U64_FROM_F64(a) & 0xFFFFFFFF00000000ULL);
        s = 1.0 + t * z;
        return t + a * (s + t * v);
    }
}

func s32
ieee754_rem_pi_over_2(f64 x, f64 *y)
{
    s64 sx = MATS_S64_FROM_F64(x);
    s32 hx = sx >> 32;
    s32 ix = hx & 0x7FFFFFFF;

    if (ix <= 0x3FE921FB)
    {
        // NOTE(michiel): No reduction needed
        y[0] = x;
        y[1] = 0;
        return 0;
    }
    if (ix < 0x4002D97C)
    {
        if (hx > 0)
        {
            f64 z = x - gPiOver2F64_1;
            if (ix != 0x3FF921FB)
            {
                y[0] = z - gPiOver2F64_1t;
                y[1] = (z - y[0]) - gPiOver2F64_1t;
            }
            else
            {
                z -= gPiOver2F64_2;
                y[0] = z - gPiOver2F64_2t;
                y[1] = (z - y[0]) - gPiOver2F64_2t;
            }
            return 1;
        }
        else
        {
            f64 z = x + gPiOver2F64_1;
            if (ix != 0x3FF921FB)
            {
                y[0] = z + gPiOver2F64_1t;
                y[1] = (z - y[0]) + gPiOver2F64_1t;
            }
            else
            {
                z += gPiOver2F64_2;
                y[0] = z + gPiOver2F64_2t;
                y[1] = (z - y[0]) + gPiOver2F64_2t;
            }
            return -1;
        }
    }
    if (ix <= 0x413921FB)
    {
        f64 t = absolute64(x);
        s32 n = (s32)(t * g2OverPiF64 + 0.5);
        f64 fn = (f64)n;
        f64 r = t - fn * gPiOver2F64_1;
        f64 w = fn * gPiOver2F64_1t;

        if ((n < 32) && (ix != gNPiOver2F64_hw[n - 1]))
        {
            y[0] = r - w;
        }
        else
        {
            s32 j = ix >> 20;
            y[0] = r - w;
            u32 high = MATS_U64_FROM_F64(y[0]) >> 32;
            s32 i = j - ((high >> 20) & 0x7FF);
            if (i > 16)
            {
                t = r;
                w = fn * gPiOver2F64_2;
                r = t - w;
                w = fn * gPiOver2F64_2t - ((t - r) - w);
                y[0] = r - w;
                high = MATS_U64_FROM_F64(y[0]) >> 32;
                i = j - ((high >> 20) & 0x7FF);
                if (i > 49)
                {
                    t = r;
                    w = fn * gPiOver2F64_3;
                    r = t - w;
                    w = fn * gPiOver2F64_3t - ((t - r) - w);
                    y[0] = r - w;
                }
            }
        }
        y[1] = (r - y[0]) - w;

        if (hx < 0)
        {
            y[0] = -y[0];
            y[1] = -y[1];
            return -n;
        }
        else
        {
            return n;
        }
    }

    if (ix >= 0x7FF00000)
    {
        y[0] = y[1] = x - x;
        return 0;
    }

    u32 lowZ = sx & 0xFFFFFFFF;
    s32 e0 = (s32)((ix >> 20) - 1046); // NOTE(michiel): ilogb(z) - 23
    s32 highZ = ix - (e0 << 20);
    f64 z = MATS_F64_FROM_U64(((u64)highZ << 32) | (u64)lowZ);

    f64 tx[3];
    tx[0] = (f64)((s32)z);
    z = (z - tx[0]) * g2pow24F64;
    tx[1] = (f64)((s32)z);
    z = (z - tx[1]) * g2pow24F64;
    tx[2] = z;

    s32 nx = 3;
    while (tx[nx - 1] == 0) {
        --nx;
    }

    s32 n = kernel_rem_pi_over_2(tx, y, e0, nx, 2, g2OverPiF64_Table);

    if (hx < 0)
    {
        y[0] = -y[0];
        y[1] = -y[1];
        return -n;
    }
    else
    {
        return n;
    }
}

func f64
cos64(f64 x)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;

    f64 result;
    if (top <= 0x3FE921FB)
    {
        result = kernel_cos64(x, 0.0);
    }
    else if (top >= 0x7FF00000)
    {
        result = x - x;
    }
    else
    {
        f64 y[2];
        s32 n = ieee754_rem_pi_over_2(x, y);
        switch (n & 0x3)
        {
            case 0 : { result =  kernel_cos64(y[0], y[1]); } break;
            case 1 : { result = -kernel_sin64(y[0], y[1], 1); } break;
            case 2 : { result = -kernel_cos64(y[0], y[1]); } break;
            default: { result =  kernel_sin64(y[0], y[1], 1); } break;
        }
    }

    return result;
}

func f64
sin64(f64 x)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;

    f64 result;
    if (top <= 0x3FE921FB)
    {
        result = kernel_sin64(x, 0.0, 0);
    }
    else if (top >= 0x7FF00000)
    {
        result = x - x;
    }
    else
    {
        f64 y[2];
        s32 n = ieee754_rem_pi_over_2(x, y);
        switch (n & 0x3)
        {
            case 0 : { result =  kernel_sin64(y[0], y[1], 1); } break;
            case 1 : { result =  kernel_cos64(y[0], y[1]); } break;
            case 2 : { result = -kernel_sin64(y[0], y[1], 1); } break;
            default: { result = -kernel_cos64(y[0], y[1]); } break;
        }
    }

    return result;
}

func SinCos64
sincos64(f64 x)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;

    SinCos64 result;
    if (top <= 0x3FE921FB)
    {
        result.cos = kernel_cos64(x, 0.0);
        result.sin = kernel_sin64(x, 0.0, 0);
    }
    else if (top >= 0x7FF00000)
    {
        result.cos = result.sin = x - x;
    }
    else
    {
        f64 y[2];
        s32 n = ieee754_rem_pi_over_2(x, y);
        f64 c = kernel_cos64(y[0], y[1]);
        f64 s = kernel_sin64(y[0], y[1], 1);
        switch (n & 0x3)
        {
            case 0 : {} break;
            case 1 : {
                f64 t = c;
                c = -s;
                s = t;
            } break;
            case 2 : {
                c = -c;
                s = -s;
            } break;
            default: {
                f64 t = c;
                c = s;
                s = -t;
            } break;
        }
        result.cos = c;
        result.sin = s;
    }

    return result;
}

func f64
tan64(f64 x)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;

    f64 result;
    if (top <= 0x3FE921FB)
    {
        result = kernel_tan64(x, 0.0, 1);
    }
    else if (top >= 0x7FF00000)
    {
        result = x - x;
    }
    else
    {
        f64 y[2];
        s32 n = ieee754_rem_pi_over_2(x, y);
        n = 1 - ((n & 1) << 1); // NOTE(michiel): +/- 1
        result = kernel_tan64(y[0], y[1], n);
    }

    return result;
}

func f64
acos64(f64 x)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;

    if (top >= 0x3FF00000)
    {
        u32 low = ix & 0xFFFFFFFF;
        if (((top - 0x3FF00000) | low) == 0)
        {
            if (ix > 0) {
                return 0.0;
            } else {
                return gPiF64 + 2.0 * gPiOver2F64_lo;
            }
        }
        return (x - x) / (x - x);
    }

#define pS0  1.66666666666666657415e-01    /* 0x3FC5555555555555 */
#define pS1 -3.25565818622400915405e-01    /* 0xBFD4D61203EB6F7D */
#define pS2  2.01212532134862925881e-01    /* 0x3FC9C1550E884455 */
#define pS3 -4.00555345006794114027e-02    /* 0xBFA48228B5688F3B */
#define pS4  7.91534994289814532176e-04    /* 0x3F49EFE07501B288 */
#define pS5  3.47933107596021167570e-05    /* 0x3F023DE10DFDF709 */
#define qS1 -2.40339491173441421878e+00    /* 0xC0033A271C8A2D4B */
#define qS2  2.02094576023350569471e+00    /* 0x40002AE59C598AC8 */
#define qS3 -6.88283971605453293030e-01    /* 0xBFE6066C1B8D0159 */
#define qS4  7.70381505559019352791e-02    /* 0x3FB3B8C5B12E9282 */

    if (top < 0x3FE00000)
    {
        // NOTE(michiel): |x| < 0.5
        if (top <= 0x3C600000) {
            return gPiOver2F64 + gPiOver2F64_lo;
        }
        f64 z = x * x;
        f64 p = z * (pS0 + z * (pS1 + z * (pS2 + z * (pS3 + z * (pS4 + z * pS5)))));
        f64 q = 1.0 + z * (qS1 + z * (qS2 + z * (qS3 + z * qS4)));
        f64 r = p / q;
        return gPiOver2F64 - (x - (gPiOver2F64_lo - x * r));
    }
    else if (ix < 0)
    {
        f64 z = (1.0 + x) * 0.5;
        f64 p = z * (pS0 + z * (pS1 + z * (pS2 + z * (pS3 + z * (pS4 + z * pS5)))));
        f64 q = 1.0 + z * (qS1 + z * (qS2 + z * (qS3 + z * qS4)));
        f64 s = sqrt64(z);
        f64 r = p / q;
        f64 w = r * s - gPiOver2F64_lo;
        return gPiF64 - 2.0 * (s + w);
    }
    else
    {
        f64 z = (1.0 - x) * 0.5;
        f64 s = sqrt64(z);
        f64 df = MATS_F64_FROM_U64(MATS_U64_FROM_F64(s) & 0xFFFFFFFF00000000ULL);
        f64 c = (z - df * df) / (s + df);
        f64 p = z * (pS0 + z * (pS1 + z * (pS2 + z * (pS3 + z * (pS4 + z * pS5)))));
        f64 q = 1.0 + z * (qS1 + z * (qS2 + z * (qS3 + z * qS4)));
        f64 r = p / q;
        f64 w = r * s + c;
        return 2.0 * (df + w);
    }
#undef qS4
#undef qS3
#undef qS2
#undef qS1
#undef pS5
#undef pS4
#undef pS3
#undef pS2
#undef pS1
#undef pS0
}

func f64
asin64(f64 x)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;

#define pS0  1.66666666666666657415e-01    /* 0x3FC5555555555555 */
#define pS1 -3.25565818622400915405e-01    /* 0xBFD4D61203EB6F7D */
#define pS2  2.01212532134862925881e-01    /* 0x3FC9C1550E884455 */
#define pS3 -4.00555345006794114027e-02    /* 0xBFA48228B5688F3B */
#define pS4  7.91534994289814532176e-04    /* 0x3F49EFE07501B288 */
#define pS5  3.47933107596021167570e-05    /* 0x3F023DE10DFDF709 */
#define qS1 -2.40339491173441421878e+00    /* 0xC0033A271C8A2D4B */
#define qS2  2.02094576023350569471e+00    /* 0x40002AE59C598AC8 */
#define qS3 -6.88283971605453293030e-01    /* 0xBFE6066C1B8D0159 */
#define qS4  7.70381505559019352791e-02    /* 0x3FB3B8C5B12E9282 */
    if (top >= 0x3FF00000)
    {
        u32 low = ix & 0xFFFFFFFF;
        if (((top - 0x3FF00000) | low) == 0)
        {
            return x * gPiOver2F64 + x * gPiOver2F64_lo;
        }
        return (x - x) / (x - x);
    }
    else if (top < 0x3FE00000)
    {
        if (top < 0x3E400000)
        {
            return x;    /* generate inexact */
        }
        else
        {
            f64 t = x * x;
            f64 p = t * (pS0 + t * (pS1 + t * (pS2 + t * (pS3 + t * (pS4 + t * pS5)))));
            f64 q = 1.0 + t * (qS1 + t * (qS2 + t * (qS3 + t * qS4)));
            f64 w = p / q;
            return x + x * w;
        }
    }

    f64 w = 1.0 - absolute64(x);
    f64 t = 0.5 * w;
    f64 p = t * (pS0 + t * (pS1 + t * (pS2 + t * (pS3 + t * (pS4 + t * pS5)))));
    f64 q = 1.0 + t * (qS1 + t * (qS2 + t * (qS3 + t * qS4)));
    f64 s = sqrt64(t);

    if (top >= 0x3FEF3333)
    {
        w = p / q;
        t = gPiOver2F64 - (2.0 * (s + s * w) - gPiOver2F64_lo);
    }
    else
    {
        w = MATS_F64_FROM_U64(MATS_U64_FROM_F64(s) & 0xFFFFFFFF00000000ULL);
        f64 c = (t - w * w) / (s + w);
        f64 r = p / q;
        p = 2.0 * s * r - (gPiOver2F64_lo - 2.0 * c);
        q = gPiOver4F64 - 2.0 * w;
        t = gPiOver4F64 - (p - q);
    }
    return (ix < 0) ? -t : t;
#undef qS4
#undef qS3
#undef qS2
#undef qS1
#undef pS5
#undef pS4
#undef pS3
#undef pS2
#undef pS1
#undef pS0
}

func f64
atan64(f64 x)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;

    s32 id = 0;
    if (top >= 0x44100000)
    {
        u32 low = ix & 0xFFFFFFFF;
        if ((top > 0x7FF00000) ||
            ((top == 0x7FF00000) && low))
        {
            return x + x;
        }
        if (top > 0) {
            return gAtanHiF64[3] + gAtanLoF64[3];
        } else {
            return -gAtanHiF64[3] - gAtanLoF64[3];
        }
    }
    if (top < 0x3FDC0000)
    {
        if (top < 0x3e200000) {
            return x;   /* generate inexact */
        }
        id = -1;
    }
    else
    {
        x = absolute64(x);
        if (top < 0x3FF30000)
        {
            if (top < 0x3FE60000) {
                id = 0;
                x = (2.0 * x - 1.0) / (2.0 + x);
            } else {
                id = 1;
                x = (x - 1.0) / (x + 1.0);
            }
        }
        else
        {
            if (top < 0x40038000) {
                id = 2;
                x = (x - 1.5) / (1.0 + 1.5 * x);
            } else {
                id = 3;
                x = -1.0 / x;
            }
        }
    }

    f64 z = x * x;
    f64 w = z * z;

#define T1  3.33333333333329318027e-01    /* 0x3FD555555555550D */
#define T2 -1.99999999998764832476e-01    /* 0xBFC999999998EBC4 */
#define T3  1.42857142725034663711e-01    /* 0x3FC24924920083FF */
#define T4 -1.11111104054623557880e-01    /* 0xBFBC71C6FE231671 */
#define T5  9.09088713343650656196e-02    /* 0x3FB745CDC54C206E */
#define T6 -7.69187620504482999495e-02    /* 0xBFB3B0F2AF749A6D */
#define T7  6.66107313738753120669e-02    /* 0x3FB10D66A0D03D51 */
#define T8 -5.83357013379057348645e-02    /* 0xBFADDE2D52DEFD9A */
#define T9  4.97687799461593236017e-02    /* 0x3FA97B4B24760DEB */
#define T10 -3.65315727442169155270e-02    /* 0xBFA2B4442C6A6C2F */
#define T11  1.62858201153657823623e-02    /* 0x3F90AD3AE322DA11 */

    f64 s1 = z * (T1 + w * (T3 + w * (T5 + w * (T7 + w * (T9 + w * T11)))));
    f64 s2 = w * (T2 + w * (T4 + w * (T6 + w * (T8 + w * T10))));

#undef T11
#undef T10
#undef T9
#undef T8
#undef T7
#undef T6
#undef T5
#undef T4
#undef T3
#undef T2
#undef T1

    if (id < 0) {
        return x - x * (s1 + s2);
    } else {
        z = gAtanHiF64[id] - ((x * (s1 + s2) - gAtanLoF64[id]) - x);
        return ix < 0 ? -z : z;
    }
}

func f64
atan2_64(f64 y, f64 x)
{
    s64 sx = MATS_S64_FROM_F64(x);
    s64 sy = MATS_S64_FROM_F64(y);

    u32 ix = (sx >> 32) & 0x7FFFFFFF;
    u32 iy = (sy >> 32) & 0x7FFFFFFF;
    u32 lx = (sx & 0xFFFFFFFF);
    u32 ly = (sy & 0xFFFFFFFF);

    if (((ix | ((lx | -lx) >> 31)) > 0x7FF00000) ||
        ((iy | ((ly | -ly) >> 31)) > 0x7FF00000))
    {
        return x + y;
    }
    if (sx == 0x3FF0000000000000LL) {
        // NOTE(michiel): x == 1.0
        return atan64(y);
    }

    s32 m = ((sx >> 62) & 0x2) | ((sy >> 63) & 0x1);

    // NOTE(michiel): y == 0
    if ((iy | ly) == 0)
    {
        switch (m)
        {
            case 0:
            case 1: { return y; } break;
            case 2: { return  gPiF64 + gTinyF64; } break;
            case 3: { return -gPiF64 - gTinyF64; } break;
        }
    }

    // NOTE(michiel): x == 0
    if ((ix | lx) == 0)
    {
        return (sy < 0) ? -gPiOver2F64 - gTinyF64 : gPiOver2F64 + gTinyF64;
    }

    // NOTE(michiel): x == inf
    if (ix == 0x7FF00000)
    {
        if (iy == 0x7FF00000)
        {
            switch (m)
            {
                case 0: { return  gPiOver4F64 + gTinyF64; } break;
                case 1: { return -gPiOver4F64 - gTinyF64; } break;
                case 2: { return  3.0 * gPiOver4F64 + gTinyF64; } break;
                case 3: { return -3.0 * gPiOver4F64 - gTinyF64; } break;
            }
        }
        else
        {
            switch (m)
            {
                case 0: { return  0.0; } break;
                case 1: { return -0.0; } break;
                case 2: { return  gPiF64 + gTinyF64; } break;
                case 3: { return -gPiF64 - gTinyF64; } break;
            }
        }
    }

    // NOTE(michiel): y == inf
    if (iy == 0x7FF00000)
    {
        return (sy < 0) ? -gPiOver2F64 - gTinyF64 : gPiOver2F64 + gTinyF64;
    }

    s32 k = (s32)(iy - ix) >> 20;
    f64 z;
    if (k > 60) {
        z = gPiOver2F64 + 0.5 * gPiF64_lo;
    } else if ((sx < 0) && (k < -60)) {
        z = 0.0;
    } else {
        z = atan64(absolute64(y / x));
    }
    switch (m)
    {
        default:
        case 0: { return z; } break;
        case 1: { return MATS_F64_FROM_U64(MATS_U64_FROM_F64(z) ^ 0x8000000000000000LL); } break;
        case 2: { return gPiF64 - (z - gPiF64_lo); } break;
        case 3: { return (z - gPiF64_lo) - gPiF64; } break;
    }
}

func f64
cosh64(f64 x)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;

    if (top >= 0x7FF00000) {
        // NOTE(michiel): Inf or NaN
        return x * x;
    }
    else if (top < 0x3FD62E43)
    {
        f64 t = expm1_64(absolute64(x));
        f64 w = 1.0 + t;
        if (top < 0x3C800000) {
            return w;
        } else {
            return 1.0 + (t * t) / (w + w);
        }
    }
    else if (top < 0x40360000)
    {
        f64 t = exp64(absolute64(x));
        return 0.5 * t + 0.5 / t;
    }
    else if (top < 0x40862E42)
    {
        return 0.5 * exp64(absolute64(x));
    }
    else
    {
        u32 low = ix & 0xFFFFFFFF;
        if ((top < 0x408633CE) ||
            ((top == 0x408633CE) && (low <= 0x8FB9F87DU)))
        {
            f64 w = exp64(0.5 * absolute64(x));
            f64 t = 0.5 * w;
            return t * w;
        }
        else
        {
            return mats_overflow64(0);
        }
    }
}

func f64
sinh64(f64 x)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;

    f64 h = (ix < 0) ? -0.5 : 0.5;

    if (top >= 0x7FF00000)
    {
        // NOTE(michiel): Inf or NaN
        return x * x;
    }
    else if (top < 0x40360000)
    {
        if (top < 0x3E300000) {
            return x;      /* generate inexact */
        }
        f64 t = expm1_64(absolute64(x));
        if (top < 0x3FF00000) {
            return h * (2.0 * t - t * t / (t + 1.0));
        } else {
            return h * (t + t / (t + 1.0));
        }
    }
    else if (top < 0x40862E42)
    {
        return h * exp64(absolute64(x));
    }
    else
    {
        u32 low = ix & 0xFFFFFFFF;
        if ((top < 0x408633CE) ||
            ((top == 0x408633CE) && (low <= 0x8FB9F87DU)))
        {
            f64 w = exp64(0.5 * absolute64(x));
            f64 t = h * w;
            return t * w;
        }
        else
        {
            return x * 1.0e307;
        }
    }
}

func SinCos64
sinhcosh64(f64 x)
{
    SinCos64 result;
    if (absolute64(x) <= 0.5)
    {
        result.cos = cosh64(x);
        result.sin = sinh64(x);
    }
    else
    {
        f64 e = exp64(x);
        f64 ei = 0.5 / e;
        e  = 0.5 * e;
        result.cos = e + ei;
        result.sin = e - ei;
    }
    return result;
}

func f64
tanh64(f64 x)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;

    f64 result;
    if (top >= 0x7FF00000)
    {
        if (ix < 0) {
            return 1.0 / x - 1.0;
        } else {
            return 1.0 / x + 1.0;
        }
    }
    else if (top < 0x40360000)
    {
        if (top < 0x3C800000) {
            return x * (1.0 + x);
        } else if (top >= 0x3FF00000) {
            f64 t = expm1_64(2.0 * absolute64(x));
            result = 1.0 - 2.0 / (t + 2.0);
        } else {
            f64 t = expm1_64(-2.0 * absolute64(x));
            result = -t / (t + 2.0);
        }
    }
    else
    {
        result = 1.0 - gTinyF64;
    }
    return (ix < 0) ? -result : result;
}

func f64
acosh64(f64 x)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32);

    if (top < 0x3FF00000)
    {
        // NOTE(michiel): x < 1.0
        return (x - x) / (x - x);
    }
    else if (top >= 0x41B00000)
    {
        if (top >= 0x7FF00000) {
            return x + x;
        } else {
            return log64(x) + gLn2F64;
        }
    }
    else if (ix == 0x3FF0000000000000LL)
    {
        return 0.0;
    }
    else if (top > 0x40000000)
    {
        f64 t = x * x;
        return log64(2.0 * x - 1.0 / (x + sqrt64(t - 1.0)));
    }
    else
    {
        f64 t = x - 1.0;
        return log1p64(t + sqrt64(2.0 * t + t * t));
    }
}

func f64
asinh64(f64 x)
{
    s64 ix = MATS_S64_FROM_F64(x);
    s32 top = (ix >> 32) & 0x7FFFFFFF;

    f64 result;
    if (top >= 0x7FF00000)
    {
        return x + x;
    }
    else if (top < 0x3E300000)
    {
        return x;  /* generate inexact */
    }
    else if (top > 0x41B00000)
    {
        result = log64(absolute64(x)) + gLn2F64;
    }
    else if (top > 0x40000000)
    {
        f64 t = absolute64(x);
        result = log64(2.0 * t + 1.0 / (sqrt64(x * x + 1.0) + t));
    }
    else
    {
        f64 t = x * x;
        result = log1p64(absolute64(x) + t / (1.0 + sqrt64(1.0 + t)));
    }

    return (ix < 0) ? -result : result;
}

func f64
atanh64(f64 x)
{
    s64 sx = MATS_S64_FROM_F64(x);

    u32 ix = (sx >> 32) & 0x7FFFFFFF;
    u32 lx = (sx & 0xFFFFFFFF);

    if ((ix | ((lx | -lx) >> 31)) > 0x7FF00000)
    {
        // NOTE(michiel): |x| > 1.0
        return (x - x) / (x - x);
    }
    else if (ix == 0x3FF00000)
    {
        return x / 0.0;
    }
    else if (ix < 0x3E300000)
    {
        return x; /* generate inexact */
    }
    else
    {
        x = absolute64(x);
        f64 result;
        if (ix < 0x3FE00000) {
            result = x + x;
            result = 0.5 * log1p64(result + result * x / (1.0 - x));
        } else {
            result = 0.5 * log1p64((x + x) / (1.0 - x));
        }
        return (sx < 0) ? -result : result;
    }
}
