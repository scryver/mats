
//
// NOTE(michiel): Helpers
//
global const f64 gPiReduce1 = 3.140625;
global const f64 gPiReduce2 = 9.67502593994140625E-4;
global const f64 gPiReduce3 = 1.509957990978376432E-7;

internal f32
reduce_pi32(f32 x)
{
    f32 t = x / F32_PI;
    if (t >= 0.0f) {
        t += 0.5f;
    } else {
        t -= 0.5f;
    }

    t = (f32)(s32)t;  /* the multiple */
    f32 result = ((x - t * gPiReduce1) - t * gPiReduce2) - t * gPiReduce3;
    return result;
}

//
//
//

struct c32
{
    f32 real;
    f32 imag;
};

global const c32 gImaginary = {0, 1};

internal c32
complex32(f32 real, f32 imag)
{
    c32 result;
    result.real = real;
    result.imag = imag;
    return result;
}

internal c32
operator -(c32 a)
{
    c32 result;
    result.real = -a.real;
    result.imag = -a.imag;
    return result;
}

internal c32 &
operator +=(c32 &a, c32 b)
{
    a.real += b.real;
    a.imag += b.imag;
    return a;
}

internal c32
operator +(c32 a, c32 b)
{
    c32 result = a;
    result += b;
    return result;
}

internal c32 &
operator -=(c32 &a, c32 b)
{
    a.real -= b.real;
    a.imag -= b.imag;
    return a;
}

internal c32
operator -(c32 a, c32 b)
{
    c32 result = a;
    result -= b;
    return result;
}

internal c32
operator *(c32 a, c32 b)
{
    c32 result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

internal c32 &
operator *=(c32 &a, c32 b)
{
    a = a * b;
    return a;
}

internal c32
operator /(c32 a, c32 b)
{
    f32 divisor = (b.real * b.real + b.imag * b.imag);
    c32 result;
    result.real = (a.real * b.real + a.imag * b.imag) / divisor;
    result.imag = (a.imag * b.real - a.real * b.imag) / divisor;
    return result;
}

internal c32 &
operator /=(c32 &a, c32 b)
{
    a = a / b;
    return a;
}

internal c32 &
operator +=(c32 &a, f32 b)
{
    a.real += b;
    return a;
}

internal c32
operator +(c32 a, f32 b)
{
    c32 result = a;
    result.real += b;
    return result;
}

internal c32
operator +(f32 a, c32 b)
{
    c32 result = b;
    result.real += a;
    return result;
}

internal c32 &
operator -=(c32 &a, f32 b)
{
    a.real -= b;
    return a;
}

internal c32
operator -(c32 a, f32 b)
{
    c32 result = a;
    result.real -= b;
    return result;
}

internal c32
operator -(f32 a, c32 b)
{
    c32 result = b;
    result.real -= a;
    return result;
}

internal c32 &
operator *=(c32 &a, f32 b)
{
    a.real *= b;
    a.imag *= b;
    return a;
}

internal c32
operator *(c32 a, f32 b)
{
    c32 result = a;
    result *= b;
    return result;
}

internal c32
operator *(f32 a, c32 b)
{
    c32 result = b;
    result *= a;
    return result;
}

internal c32 &
operator /=(c32 &a, f32 b)
{
    a.real /= b;
    a.imag /= b;
    return a;
}

internal c32
operator /(c32 a, f32 b)
{
    c32 result = a;
    result /= b;
    return result;
}

internal c32
operator /(f32 a, c32 b)
{
    c32 result = b;
    result /= a;
    return result;
}

internal c32
conjugate32(c32 c)
{
    c32 result = c;
    result.imag = -result.imag;
    return result;
}

internal f32
absolute32(c32 c)
{
    return hypot32(c.real, c.imag);
}

internal f32
argument32(c32 c)
{
    return atan2_32(c.imag, c.real);
}

//
// NOTE(michiel): Complex elementary functions
//

internal c32
sqrt32(c32 c)
{
    f32 x = c.real;
    f32 y = c.imag;

    if (y == 0.0f)
    {
        if (x < 0.0f) {
            return complex32(0, sqrt32(-x));
        } else if (x == 0.0f) {
            return complex32(0, y);
        } else {
            return complex32(sqrt32(x), y);
        }
    }
    else if (x == 0.0f)
    {
        f32 r = absolute32(y);
        r = sqrt32(0.5f * r);
        return complex32(r, (y > 0.0f) ? r : -r);
    }

    f32 scale = 1.0f;
    if ((absolute32(x) > 4.0f) || (absolute32(y) > 4.0f))
    {
        // NOTE(michiel): Rescale to avoid internal overflow
        x *= 0.25f;
        y *= 0.25f;
        scale = 2.0f;
    }
    else
    {
        x *= 6.7108864e7f; // NOTE(michiel): 2^26
        y *= 6.7108864e7f;
        scale = 1.220703125e-4f; // NOTE(michiel): 2^-13
    }

    f32 r = absolute32(complex32(x, y));
    f32 t;
    if (x > 0.0f)
    {
        t = sqrt32(0.5f * r + 0.5f * x);
        r = scale * absolute32((0.5f * y) / t);
        t *= scale;
    }
    else
    {
        r = sqrt32(0.5f * r - 0.5f * x);
        t = scale * absolute32((0.5f * y) / r);
        r *= scale;
    }

    c32 result = complex32(t, (y < 0.0f) ? -r : r);
    return result;
}

internal c32
exp32(c32 c)
{
    SinCos32 cs = sincos32(c.imag);
    f32 r = exp32(c.real);
    c32 result = complex32(r * cs.cos, r * cs.sin);
    return result;
}

internal c32
log32(c32 c)
{
    f32 rr = absolute32(c);
    f32 p  = log32(rr);
    rr = argument32(c);
    c32 result = complex32(p, rr);
    return result;
}

internal c32
log10_32(c32 c)
{
    f32 rr = absolute32(c);
    f32 p  = log10_32(rr);
    rr = argument32(c) * gInvLn10F32;
    c32 result = complex32(p, rr);
    return result;
}

internal c32
pow32(c32 x, c32 y)
{

    f32 real = y.real;
    f32 imag = y.imag;

    f32 absX = absolute32(x);
    if (absX == 0.0f) {
        return complex32(0, 0);
    }

    f32 argX = argument32(x);
    f32 r = pow32(absX, real);
    f32 theta = real * argX;

    if (imag != 0.0f)
    {
        r = r * exp32(-imag * argX);
        theta = theta + imag * log32(absX);
    }

    SinCos32 th = sincos32(theta);
    c32 result = complex32(r * th.cos, r * th.sin);
    return result;
}

//
// NOTE(michiel): Complex trigonometric functions
//

internal c32
cos32(c32 c)
{
    SinCos32 cs = sincos32(c.real);
    SinCos32 csh = sinhcosh32(c.imag);
    c32 result = complex32(cs.cos * csh.cos, -(cs.sin * csh.sin));
    return result;
}

internal c32
sin32(c32 c)
{
    SinCos32 cs = sincos32(c.real);
    SinCos32 csh = sinhcosh32(c.imag);
    c32 result = complex32(cs.sin * csh.cos, cs.cos * csh.sin);
    return result;
}

internal f32
tan_serie32(c32 c)
{
    // NOTE(michiel): Taylor-series expansion for cosh(2y) - cos(2x)
    f32 x = absolute32(2.0f * c.real);
    f32 y = absolute32(2.0f * c.imag);

    x = reduce_pi32(x);

    x = x * x;
    y = y * y;

    f32 x2 = 1.0f;
    f32 y2 = 1.0f;
    f32 f  = 1.0f;
    f32 rn = 0.0f;
    f32 d  = 0.0f;

    // NOTE(michiel): While error is greater than machine precision
    f32 t;
    do
    {
        rn += 1.0f;
        f  *= rn;
        rn += 1.0f;
        f  *= rn;
        x2 *= x;
        y2 *= y;
        t   = y2 + x2;
        t  /= f;
        d  += f;

        rn += 1.0f;
        f  *= rn;
        rn += 1.0f;
        f  *= rn;
        x2 *= x;
        y2 *= y;
        t   = y2 - x2;
        t  /= f;
        d  += t;
    } while (absolute32(t / d) > 3.0e-8f);
    return d;
}

internal c32
tan32(c32 c)
{
    SinCos32 cs  = sincos32(2.0f * c.real);
    SinCos32 csh = sinhcosh32(2.0f * c.imag);
    f32 d = cs.cos + csh.cos;

    if (absolute32(d) < 0.25f) {
        d = tan_serie32(c);
    }

    if (d == 0.0f) {
        return complex32(gHugeF32, gHugeF32);
    }

    c32 result = complex32(cs.sin / d, csh.sin / d);
    return result;
}

internal c32
asin32(c32 c)
{
    f32 x = c.real;
    f32 y = c.imag;

    c32 ct = c * gImaginary;

    c32 zz = complex32((x - y) * (x + y), 2.0f * x * y);
    zz.real = 1.0f - zz.real;
    zz.imag = -zz.imag;

    c32 z2 = sqrt32(zz);

    zz = ct + z2;
    zz = log32(zz);

    c32 result = zz * (-gImaginary);
    return result;
}

internal c32
acos32(c32 c)
{
    c32 result = asin32(c);
    result = (F32_PI_OVER_2 - result.real) - complex32(0.0f, result.imag);
    return result;
}

internal c32
atan32(c32 c)
{
    c32 result;
    f32 x = c.real;
    f32 y = c.imag;

    if ((x == 0.0f) && (y > 1.0f)) {
        return complex32(gHugeF32, gHugeF32);
    }

    f32 x2 = x * x;
    f32 a = 1.0f - x2 - (y * y);
    if (a == 0.0f) {
        return complex32(gHugeF32, gHugeF32);
    }

    f32 t = 0.5f * atan2_32(2.0f * x, a);
    f32 w = reduce_pi32(t);

    t = y - 1.0f;
    a = x2 + (t * t);
    if (a == 0.0f) {
        return complex32(gHugeF32, gHugeF32);
    }

    t = y + 1.0f;
    a = (x2 + (t * t)) / a;
    result = complex32(w, 0.25f * log(a));
    return result;
}

//
// NOTE(michiel): Complex hyperbolic functions
//

internal c32
cosh32(c32 c)
{
    SinCos32 csh = sinhcosh32(c.real);
    SinCos32 cs = sincos32(c.imag);
    c32 result = complex32(cs.cos * csh.cos, cs.sin * csh.sin);
    return result;
}

internal c32
sinh32(c32 c)
{
    SinCos32 csh = sinhcosh32(c.real);
    SinCos32 cs = sincos32(c.imag);
    c32 result = complex32(csh.sin * cs.cos, csh.cos * cs.sin);
    return result;
}

internal c32
tanh32(c32 c)
{
    SinCos32 csh = sinhcosh32(2.0f * c.real);
    SinCos32 cs  = sincos32(2.0f * c.imag);
    f32 d = cs.cos + csh.cos;
    c32 result = complex32(csh.sin / d, cs.sin / d);
    return result;
}

internal c32
acosh32(c32 c)
{
    c32 result = log32(c + sqrt32(c + 1.0f) * sqrt32(c - 1.0f));
    return result;
}

internal c32
asinh32(c32 c)
{
    c32 result = complex32(0, -1.0f) * asin32(c * complex32(0, 1.0f));
    return result;
}

internal c32
atanh32(c32 c)
{
    c32 result = complex32(0, -1.0f) * atan32(c * complex32(0, 1));
    return result;
}
