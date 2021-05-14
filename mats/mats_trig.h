
internal f32
reduce_fast_pi4(f32 x, int *np)
{
    /* Use scaled float to int conversion with explicit rounding.
       hpi_inv is prescaled by 2^24 so the quadrant ends up in bits 24..31.
       This avoids inaccuracies introduced by truncating negative values.  */
    f32 r = x * 0x1.45F306p+23f;
    s32 n = ((s32)r + 0x800000) >> 24;
    *np = n;
    return x - n * 0x1.921FB4p0f;
}

internal f32
sinf_poly_q0(f32 x, f32 x2)
{
    // NOTE(michiel): x - s1*x^3 + s2*x^5 - s3*x^7
    f32 x3 = x * x2;                                // x^3
    f32 s1 = 0x1.110760p-7f - x2 * 0x1.994eb2p-13f; // s2 - s3*x^2
    f32 x7 = x3 * x2;                               // x^5
    f32 s = x - x3 * 0x1.555544p-3f;                // x - s1*x^3
    return s + x7 * s1;                             // x - s1*x^3 + x^5*(s2 - s3*x^2)
}

internal f32
sinf_poly_q1(f32 x2)
{
    // NOTE(michiel): c0 - c1*x^2 + c2*x^4 - c3*x^6 + c4*x^8;
    f32 x4 = x2 * x2;                                 // x^4
    f32 c2 = -0x1.6c087ep-10f + x2 * 0x1.993430p-16f; // -c3 + c4*x^2
    f32 c1 = 1.0f - x2 * 0x1.fffffep-2f;              // c0 - c1*x^2
    f32 x6 = x4 * x2;                                 // x^6
    f32 c = c1 + x4 * 0x1.55553ep-5f;                 // c0 - c1*x^2 + c2*x^4
    return c + x6 * c2;                               // c0 - c1*x^2 + c2*x^4 + x^6(-c3 + c4*x^2)
}

internal f64
reduce_fast_pi4_prec(f64 x, int *np)
{
    /* Use scaled float to int conversion with explicit rounding.
       hpi_inv is prescaled by 2^24 so the quadrant ends up in bits 24..31.
       This avoids inaccuracies introduced by truncating negative values.  */
    f64 r = x * 0x1.45F306DC9C883p+23;
    s32 n = ((s32)r + 0x800000) >> 24;
    *np = n;
    return x - n * 0x1.921FB54442D18p0;
}

internal f64
sinf_poly_q0_prec(f64 x, f64 x2)
{
    // NOTE(michiel): x - s1*x^3 + s2*x^5 - s3*x^7
    f64 x3 = x * x2;                                // x^3
    f64 s1 = 0x1.1107605230bc4p-7 - x2 * 0x1.994eb3774cf24p-13; // s2 - s3*x^2
    f64 x7 = x3 * x2;                               // x^5
    f64 s = x - x3 * 0x1.555545995a603p-3;                // x - s1*x^3
    return s + x7 * s1;                             // x - s1*x^3 + x^5*(s2 - s3*x^2)
}

internal f64
sinf_poly_q1_prec(f64 x2)
{
    // NOTE(michiel): c0 - c1*x^2 + c2*x^4 - c3*x^6 + c4*x^8;
    f64 x4 = x2 * x2;                                 // x^4
    f64 c2 = -0x1.6c087e89a359dp-10 + x2 * 0x1.99343027bf8c3p-16; // -c3 + c4*x^2
    f64 c1 = 1.0f - x2 * 0x1.ffffffd0c621cp-2;              // c0 - c1*x^2
    f64 x6 = x4 * x2;                                 // x^6
    f64 c = c1 + x4 * 0x1.55553e1068f19p-5;                 // c0 - c1*x^2 + c2*x^4
    return c + x6 * c2;                               // c0 - c1*x^2 + c2*x^4 + x^6(-c3 + c4*x^2)
}

internal f32
cos32_fast(f32 y)
{
    // NOTE(michiel): Sloppy
    i_expect(absolute32(y) < 120.0f);
    f32 x = y;
    f32 result;
    if (abstop12_(y) < abstop12_(0x1p-12f))
    {
        result = 1.0f;
    }
    else
    {
        int n = 0;
        f32 xm = (abstop12_(y) < abstop12_(0.25f * F32_PI)) ? x : reduce_fast_pi4(x, &n);
        f32 x2 = xm * xm;
        switch (n & 3)
        {
            default:
            case 0: { result =  sinf_poly_q1(x2); } break;
            case 1: { result = -sinf_poly_q0(xm, x2); } break;
            case 2: { result = -sinf_poly_q1(x2); } break;
            case 3: { result =  sinf_poly_q0(xm, x2); } break;
        }
    }

    return result;
}

internal f32
sin32_fast(f32 y)
{
    // NOTE(michiel): Sloppy
    i_expect(absolute32(y) < 120.0f);
    f32 x = y;

    f32 result;
    if (abstop12_(y) < abstop12_(0x1p-12f))
    {
        result = y;
    }
    else
    {
        int n = 0;
        f32 xm = (abstop12_(y) < abstop12_(0.25f * F32_PI)) ? x : reduce_fast_pi4(x, &n);
        f32 x2 = xm * xm;
        switch (n & 3)
        {
            default:
            case 0: { result =  sinf_poly_q0(xm, x2); } break;
            case 1: { result =  sinf_poly_q1(x2); } break;
            case 2: { result = -sinf_poly_q0(xm, x2); } break;
            case 3: { result = -sinf_poly_q1(x2); } break;
        }
    }

    return result;
}

internal f32
cos32(f32 y)
{
    i_expect(absolute32(y) < 120.0f);
    f32 x = y;
    f32 result;
    if (abstop12_(y) < abstop12_(0x1p-12f))
    {
        result = 1.0f;
    }
    else
    {
        int n = 0;
        f64 xm = (abstop12_(y) < abstop12_(0.25f * F32_PI)) ? (f64)x : reduce_fast_pi4_prec((f64)x, &n);
        f64 x2 = xm * xm;
        switch (n & 3)
        {
            default:
            case 0: { result =  sinf_poly_q1_prec(x2); } break;
            case 1: { result = -sinf_poly_q0_prec(xm, x2); } break;
            case 2: { result = -sinf_poly_q1_prec(x2); } break;
            case 3: { result =  sinf_poly_q0_prec(xm, x2); } break;
        }
    }

    return result;
}

internal f32
sin32(f32 y)
{
    i_expect(absolute32(y) < 120.0f);
    f32 x = y;

    f32 result;
    if (abstop12_(y) < abstop12_(0x1p-12f))
    {
        result = y;
    }
    else
    {
        int n = 0;
        f64 xm = (abstop12_(y) < abstop12_(0.25f * F32_PI)) ? (f64)x : reduce_fast_pi4_prec((f64)x, &n);
        f64 x2 = xm * xm;
        switch (n & 3)
        {
            default:
            case 0: { result =  sinf_poly_q0_prec(xm, x2); } break;
            case 1: { result =  sinf_poly_q1_prec(x2); } break;
            case 2: { result = -sinf_poly_q0_prec(xm, x2); } break;
            case 3: { result = -sinf_poly_q1_prec(x2); } break;
        }
    }

    return result;
}

internal v2
sincos32_fast(f32 y)
{
    // NOTE(michiel): x = cos, y = sin
    i_expect(absolute32(y) < 120.0f);
    f32 x = y;

    v2 result;
    if (abstop12_(y) < abstop12_(0x1p-12f))
    {
        result.x = 1.0f;
        result.y = y;
    }
    else
    {
        int n = 0;
        f32 xm = (abstop12_(y) < abstop12_(0.25f * F32_PI)) ? x : reduce_fast_pi4(x, &n);
        f32 x2 = xm * xm;
        f32 sinP = sinf_poly_q0(xm, x2);
        f32 cosP = sinf_poly_q1(x2);
        switch (n & 3)
        {
            case 0: { result.x =  cosP; result.y =  sinP; } break;
            case 1: { result.x = -sinP; result.y =  cosP; } break;
            case 2: { result.x = -cosP; result.y = -sinP; } break;
            case 3: { result.x =  sinP; result.y = -cosP; } break;
        }
    }

    return result;
}

internal v2
sincos32(f32 y)
{
    // NOTE(michiel): x = cos, y = sin
    i_expect(absolute32(y) < 120.0f);
    f32 x = y;

    v2 result;
    if (abstop12_(y) < abstop12_(0x1p-12f))
    {
        result.x = 1.0f;
        result.y = y;
    }
    else
    {
        int n = 0;
        f64 xm = (abstop12_(y) < abstop12_(0.25f * F32_PI)) ? x : reduce_fast_pi4_prec(x, &n);
        f64 x2 = xm * xm;
        f64 sinP = sinf_poly_q0_prec(xm, x2);
        f64 cosP = sinf_poly_q1_prec(x2);
        switch (n & 3)
        {
            case 0: { result.x =  cosP; result.y =  sinP; } break;
            case 1: { result.x = -sinP; result.y =  cosP; } break;
            case 2: { result.x = -cosP; result.y = -sinP; } break;
            case 3: { result.x =  sinP; result.y = -cosP; } break;
        }
    }

    return result;
}

internal f32
tan32_kernel(f32 x, s32 mod)
{
	s32 ix = (s32)u32f32(x).u;
    u32 hx = ix & 0x7FFFFFFF;
    if (hx < 0x31800000)
    {
        if ((s32)x == 0) {
            if ((hx | (mod + 1)) == 0) {
                return 1.0f / u32f32(hx).f;
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
        u32 uz = u32f32(z).u & 0xFFFFF000;
        z = u32f32(uz).f;
        v = r - (z - x);
        f32 t = -1.0f / w;
        f32 a = t;
        t = u32f32(u32f32(t).u & 0xFFFFF000).f;
        s = 1.0f + t * z;
        return t + a * (s + t * v);
    }
}

internal f32
tan32(f32 y)
{
    i_expect(absolute32(y) < 120.0f);
    f32 result;
    u32 uy = u32f32(y).u & 0x7FFFFFFF;
    if (uy < 0x3F490FDA)
    {
        result = tan32_kernel(y, 1);
    }
    else
    {
        int n = 0;
        f32 x = reduce_fast_pi4_prec(y, &n);
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

internal f32
acos32(f32 x)
{
#define pi        3.1415925026e+00f
#define pio2_hi   1.5707962513e+00f
#define pio2_lo   7.5497894159e-08f

	s32 ix = (s32)u32f32(x).u;
    u32 hx = ix & 0x7FFFFFFF;

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
        f32 df = u32f32(u32f32(s).u & 0xFFFFF000).f;
        f32 c = (z - df * df) / (s + df);
        poly(z);
        f32 w = r * s + c;
        return 2.0f * (df + w); // NOTE(michiel): x >= 0.5f
    }

#undef pio2_hi
#undef pio2_lo
#undef pi
}

internal f32
asin32(f32 x)
{
#define pio2_hi 1.57079637050628662109375f
#define pio2_lo -4.37113900018624283e-8f
#define pio4_hi 0.785398185253143310546875f

    s32 ix = (s32)u32f32(x).u;
    u32 hx = ix & 0x7FFFFFFF;

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

    f32 w = 1.0f - u32f32(hx).f;
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
        u32 iw = u32f32(w).u;
        w = u32f32(iw & 0xFFFFF000).f;
        f32 c = (t - w * w) / (s + w);
        p = 2.0f * s * r - (pio2_lo - 2.0f * c);
        q = pio4_hi - 2.0f * w;
        t = pio4_hi - (p - q);
    }

    if (ix > 0) {
        return t;
    } else {
        return -t;
    }

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

global const f32 gAtanHiF32[] = {
    4.6364760399e-01, /* atan(0.5)hi 0x3eed6338 */
    7.8539812565e-01, /* atan(1.0)hi 0x3f490fda */
    9.8279368877e-01, /* atan(1.5)hi 0x3f7b985e */
    1.5707962513e+00, /* atan(inf)hi 0x3fc90fda */
};

global const f32 gAtanLoF32[] = {
    5.0121582440e-09, /* atan(0.5)lo 0x31ac3769 */
    3.7748947079e-08, /* atan(1.0)lo 0x33222168 */
    3.4473217170e-08, /* atan(1.5)lo 0x33140fb4 */
    7.5497894159e-08, /* atan(inf)lo 0x33a22168 */
};

internal f32
atan32(f32 x)
{
    s32 ix = (s32)u32f32(x).u;
    u32 hx = ix & 0x7FFFFFFF;

    if (hx >= 0x50800000)
    {
        // NOTE(michiel): |x| >= 2^34
        if (FLT_UWORD_IS_NAN(hx)) {
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
            x = u32f32(hx).f;
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

internal f32
atan2_32(f32 y, f32 x)
{
    s32 ix = (s32)u32f32(x).u;
    u32 hx = ix & 0x7FFFFFFF;
    s32 iy = (s32)u32f32(y).u;
    u32 hy = iy & 0x7FFFFFFF;
    if (FLT_UWORD_IS_NAN(hx) || FLT_UWORD_IS_NAN(hy)) {
        return x + y;
    }
    if (ix == 0x3F800000) {
        return atan32(y); // NOTE(michiel): x == 1.0f
    }
    s32 m = ((iy >> 31) & 1) | ((ix >> 30) & 2); // NOTE(michiel): 2 * sign(x) + sign(y)
    if (FLT_UWORD_IS_ZERO(hy))
    {
        switch (m)
        {
            case 0:
            case 1: { return y; } break;                  // NOTE(michiel): atan(+/-0,+anything) = +/-0
            case 2: { return  gPiF32 + gTinyF32; } break; // NOTE(michiel): atan(+0,  -anything) =  pi
            case 3: { return -gPiF32 - gTinyF32; } break; // NOTE(michiel): atan(-0,  -anything) = -pi
        }
    }

    if (FLT_UWORD_IS_ZERO(hx))
    {
        return (iy < 0) ? -gPiOver2F32 - gTinyF32 : gPiOver2F32 + gTinyF32;
    }

    if (FLT_UWORD_IS_INFINITE(hx))
    {
        if (FLT_UWORD_IS_INFINITE(hy))
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

    if (FLT_UWORD_IS_INFINITE(hy)) {
        return (iy < 0) ? -gPiOver2F32 - gTinyF32 : gPiOver2F32 + gTinyF32;
    }

    f32 result;
    s32 k = ((s32)hy - (s32)hx) >> 23;
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
        case 1: { result = u32f32(u32f32(result).u ^ 0x80000000).f; } break;
        case 2: { result = gPiF32 - (result - gPiF32_lo); } break;
        case 3: { result = (result - gPiF32_lo) - gPiF32; } break;
    }
    return result;
}

