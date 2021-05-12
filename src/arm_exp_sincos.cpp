
internal f32
reduce_fast_exp(f32 x, int *np)
{
    f32 hpi_inv = 0x1.45F306p+23f;
    /* Use scaled float to int conversion with explicit rounding.
       hpi_inv is prescaled by 2^24 so the quadrant ends up in bits 24..31.
       This avoids inaccuracies introduced by truncating negative values.  */
    f32 r = x * hpi_inv;
    s32 n = ((s32)r + 0x800000) >> 24;
    *np = n;
    //return x - n * 0x1.921FB54442D18p0;
    return x - n * 0x1.921FB4p0f;
}

internal f32
sinf_exp_poly_q0(f32 x, f32 x2)
{
    // NOTE(michiel): x - s1*x^3 + s2*x^5 - s3*x^7
    f32 x3 = x * x2;                                // x^3
    f32 s1 = 0x1.110760p-7f - x2 * 0x1.994eb2p-13f; // s2 - s3*x^2
    f32 x7 = x3 * x2;                               // x^5
    f32 s = x - x3 * 0x1.555544p-3f;                // x - s1*x^3
    return s + x7 * s1;                             // x - s1*x^3 + x^5*(s2 - s3*x^2)
}

internal f32
sinf_exp_poly_q1(f32 x2)
{
    // NOTE(michiel): c0 - c1*x^2 + c2*x^4 - c3*x^6 + c4*x^8;
    f32 x4 = x2 * x2;                                 // x^4
    f32 c2 = -0x1.6c087ep-10f + x2 * 0x1.993430p-16f; // -c3 + c4*x^2
    f32 c1 = 1.0f - x2 * 0x1.fffffep-2f;              // c0 - c1*x^2
    f32 x6 = x4 * x2;                                 // x^6
    f32 c = c1 + x4 * 0x1.55553ep-5f;                 // c0 - c1*x^2 + c2*x^4
    return c + x6 * c2;                               // c0 - c1*x^2 + c2*x^4 + x^6(-c3 + c4*x^2)
}

global f64 gSinCosSigns[4] = { 1.0, -1.0, -1.0, 1.0 };

internal f32
arm_exp_cosf(f32 y)
{
    i_expect(absolute(y) < 120.0f);
    f32 x = y;

    f32 result;
    if (abstop12(y) < abstop12(0x1p-12f))
    {
        result = 1.0f;
    }
    else
    {
        int n = 0;
        f32 xm = (abstop12(y) < abstop12(gPiOver4)) ? x : reduce_fast_exp(x, &n);
        f32 x2 = xm * xm;
        switch (n & 3)
        {
            default:
            case 0: { result =  sinf_exp_poly_q1(x2); } break;
            case 1: { result = -sinf_exp_poly_q0(xm, x2); } break;
            case 2: { result = -sinf_exp_poly_q1(x2); } break;
            case 3: { result =  sinf_exp_poly_q0(xm, x2); } break;
        }
    }

#if 0
    f32 result;
    if (abstop12(y) < abstop12(gPiOver4))
    {
        if (abstop12(y) < abstop12(0x1p-12f)) {
            result = 1.0f;
        } else {
            f32 x2 = x * x;
            result = sinf_exp_poly_q1(x2);
        }
    }
    else
    {
        int n = 0;
        x = reduce_fast_exp(x, &n);
        f32 x2 = x * x;
        switch (n & 3)
        {
            case 0: { result =  sinf_exp_poly_q1(x2); } break;
            case 1: { result = -sinf_exp_poly_q0(x, x2); } break;
            case 2: { result = -sinf_exp_poly_q1(x2); } break;
            case 3: { result =  sinf_exp_poly_q0(x, x2); } break;
        }
    }
#endif

#if 0
    else if (abstop12 (y) < abstop12 (120.0f))
    {
        int n = 0;
        x = reduce_fast_exp(x, &n);
        f32 x2 = x * x;
        switch (n & 3)
        {
            case 0: { result =  sinf_exp_poly_q1(x2); } break;
            case 1: { result =  sinf_exp_poly_q0(-x, x2); } break;
            case 2: { result = -sinf_exp_poly_q1(x2); } break;
            case 3: { result =  sinf_exp_poly_q0( x, x2); } break;
        }
    }
    else
    {
        // NOTE(michiel): Don't support greater values, we want speed so modulus it yourself
        result = y - y;
    }
#endif

    return result;
}

internal f32
arm_exp_sinf(f32 y)
{
    i_expect(absolute(y) < 120.0f);
    f32 x = y;

    f32 result;
    if (abstop12(y) < abstop12(0x1p-12f))
    {
        result = y;
    }
    else
    {
        int n = 0;
        f32 xm = (abstop12(y) < abstop12(gPiOver4)) ? x : reduce_fast_exp(x, &n);
        f32 x2 = xm * xm;
        switch (n & 3)
        {
            default:
            case 0: { result =  sinf_exp_poly_q0(xm, x2); } break;
            case 1: { result =  sinf_exp_poly_q1(x2); } break;
            case 2: { result = -sinf_exp_poly_q0(xm, x2); } break;
            case 3: { result = -sinf_exp_poly_q1(x2); } break;
        }
    }

#if 0
    f32 result;
    if (abstop12(y) < abstop12(gPiOver4))
    {
        f32 x2 = x * x;
        if (abstop12(y) < abstop12(0x1p-12f))
        {
            /* Force underflow for tiny y.  */
            if (abstop12(y) < abstop12(0x1p-126f)) {
                force_eval_float(x2);
            }
            result = y;
        }
        else
        {
            result = sinf_exp_poly_q0(x, x2);
        }
    }
    else
    {
        int n = 0;
        x = reduce_fast_exp(x, &n);
        f32 x2 = x * x;
        switch (n & 3)
        {
            case 0: { result =  sinf_exp_poly_q0(x, x2); } break;
            case 1: { result =  sinf_exp_poly_q1(x2); } break;
            case 2: { result = -sinf_exp_poly_q0(x, x2); } break;
            case 3: { result = -sinf_exp_poly_q1(x2); } break;
        }
    }
#endif

#if 0
    else if (abstop12(y) < abstop12(120.0f))
    {
        int n = 0;
        x = reduce_fast_exp(x, &n);
        f32 x2 = x * x;
        switch (n & 3)
        {
            case 0: { result =  sinf_exp_poly_q0( x, x2); } break;
            case 1: { result =  sinf_exp_poly_q1(x2); } break;
            case 2: { result =  sinf_exp_poly_q0(-x, x2); } break;
            case 3: { result = -sinf_exp_poly_q1(x2); } break;
        }
    }
    else
    {
        // NOTE(michiel): Don't support greater values, we want speed so modulus it yourself
        result = y - y;
    }
#endif

    return result;
}
