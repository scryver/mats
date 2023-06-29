
struct c32
{
    f32 real;
    f32 imag;
};

global const c32 gImaginary32 = {0, 1};

struct c64
{
    f64 real;
    f64 imag;
};

global const c64 gImaginary64 = {0, 1};

//
// NOTE(michiel): Helpers
//
global const f64 gPiReduce1F32 = 3.140625;
global const f64 gPiReduce2F32 = 9.67502593994140625e-4;
global const f64 gPiReduce3F32 = 1.509957990978376432e-7;

global const f64 gPiReduce1F64 = 3.14159265160560607910e0;
global const f64 gPiReduce2F64 = 1.98418714791870343106e-9;
global const f64 gPiReduce3F64 = 1.14423774522196636802e-17;

func f32
reduce_pi32(f32 x)
{
    f32 t = x / gPiF32;
    if (t >= 0.0f) {
        t += 0.5f;
    } else {
        t -= 0.5f;
    }

    t = (f32)(s32)t;  /* the multiple */
    f32 result = (f32)(((x - t * gPiReduce1F32) - t * gPiReduce2F32) - t * gPiReduce3F32);
    return result;
}

func f64
reduce_pi64(f64 x)
{
    f64 t = x / gPiF32;
    if (t >= 0.0f) {
        t += 0.5f;
    } else {
        t -= 0.5f;
    }

    t = (f64)(s64)t;  /* the multiple */
    f64 result = ((x - t * gPiReduce1F64) - t * gPiReduce2F64) - t * gPiReduce3F64;
    return result;
}

//
//
//

func c32
complex32(f32 real, f32 imag)
{
    c32 result;
    result.real = real;
    result.imag = imag;
    return result;
}

func c32
complex32_from_mag_phase(f32 magnitude, f32 phase)
{
    return complex32(magnitude * cos32(phase), magnitude * sin32(phase));
}

func c32
operator -(c32 a)
{
    c32 result;
    result.real = -a.real;
    result.imag = -a.imag;
    return result;
}

func c32 &
operator +=(c32 &a, c32 b)
{
    a.real += b.real;
    a.imag += b.imag;
    return a;
}

func c32
operator +(c32 a, c32 b)
{
    c32 result = a;
    result += b;
    return result;
}

func c32 &
operator -=(c32 &a, c32 b)
{
    a.real -= b.real;
    a.imag -= b.imag;
    return a;
}

func c32
operator -(c32 a, c32 b)
{
    c32 result = a;
    result -= b;
    return result;
}

func c32
operator *(c32 a, c32 b)
{
    c32 result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

func c32 &
operator *=(c32 &a, c32 b)
{
    a = a * b;
    return a;
}

func c32
operator /(c32 a, c32 b)
{
    f32 divisor = (b.real * b.real + b.imag * b.imag);
    c32 result;
    result.real = (a.real * b.real + a.imag * b.imag) / divisor;
    result.imag = (a.imag * b.real - a.real * b.imag) / divisor;
    return result;
}

func c32 &
operator /=(c32 &a, c32 b)
{
    a = a / b;
    return a;
}

func c32 &
operator +=(c32 &a, f32 b)
{
    a.real += b;
    return a;
}

func c32
operator +(c32 a, f32 b)
{
    c32 result = a;
    result.real += b;
    return result;
}

func c32
operator +(f32 a, c32 b)
{
    c32 result = b;
    result.real += a;
    return result;
}

func c32 &
operator -=(c32 &a, f32 b)
{
    a.real -= b;
    return a;
}

func c32
operator -(c32 a, f32 b)
{
    c32 result = a;
    result.real -= b;
    return result;
}

func c32
operator -(f32 a, c32 b)
{
    c32 result = b;
    result.real -= a;
    return result;
}

func c32 &
operator *=(c32 &a, f32 b)
{
    a.real *= b;
    a.imag *= b;
    return a;
}

func c32
operator *(c32 a, f32 b)
{
    c32 result = a;
    result *= b;
    return result;
}

func c32
operator *(f32 a, c32 b)
{
    c32 result = b;
    result *= a;
    return result;
}

func c32 &
operator /=(c32 &a, f32 b)
{
    a.real /= b;
    a.imag /= b;
    return a;
}

func c32
operator /(c32 a, f32 b)
{
    c32 result = a;
    result /= b;
    return result;
}

func c32
operator /(f32 a, c32 b)
{
    c32 result = complex32(a, 0.0);
    result /= b;
    return result;
}

func c32
conjugate32(c32 c)
{
    c32 result = c;
    result.imag = -result.imag;
    return result;
}

func f32
absolute32(c32 c)
{
    return hypot32(c.real, c.imag);
}

func f32
argument32(c32 c)
{
    return atan2_32(c.imag, c.real);
}

//
// NOTE(michiel): Complex elementary functions
//

func c32
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

func c32
exp32(c32 c)
{
    SinCos32 cs = sincos32(c.imag);
    f32 r = exp32(c.real);
    c32 result = complex32(r * cs.cos, r * cs.sin);
    return result;
}

func c32
log32(c32 c)
{
    f32 rr = absolute32(c);
    f32 p  = log32(rr);
    rr = argument32(c);
    c32 result = complex32(p, rr);
    return result;
}

func c32
log10_32(c32 c)
{
    f32 rr = absolute32(c);
    f32 p  = log10_32(rr);
    rr = argument32(c) * gInvLn10F32;
    c32 result = complex32(p, rr);
    return result;
}

func c32
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

func c32
cos32(c32 c)
{
    SinCos32 cs = sincos32(c.real);
    SinCos32 csh = sinhcosh32(c.imag);
    c32 result = complex32(cs.cos * csh.cos, -(cs.sin * csh.sin));
    return result;
}

func c32
sin32(c32 c)
{
    SinCos32 cs = sincos32(c.real);
    SinCos32 csh = sinhcosh32(c.imag);
    c32 result = complex32(cs.sin * csh.cos, cs.cos * csh.sin);
    return result;
}

func f32
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

func c32
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

func c32
asin32(c32 c)
{
    f32 x = c.real;
    f32 y = c.imag;

    c32 ct = c * gImaginary32;

    c32 zz = complex32((x - y) * (x + y), 2.0f * x * y);
    zz.real = 1.0f - zz.real;
    zz.imag = -zz.imag;

    c32 z2 = sqrt32(zz);

    zz = ct + z2;
    zz = log32(zz);

    c32 result = zz * (-gImaginary32);
    return result;
}

func c32
acos32(c32 c)
{
    c32 result = asin32(c);
    result = (gPiOver2F32 - result.real) - complex32(0.0f, result.imag);
    return result;
}

func c32
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
    result = complex32(w, 0.25f * log32(a));
    return result;
}

//
// NOTE(michiel): Complex hyperbolic functions
//

func c32
cosh32(c32 c)
{
    SinCos32 csh = sinhcosh32(c.real);
    SinCos32 cs = sincos32(c.imag);
    c32 result = complex32(cs.cos * csh.cos, cs.sin * csh.sin);
    return result;
}

func c32
sinh32(c32 c)
{
    SinCos32 csh = sinhcosh32(c.real);
    SinCos32 cs = sincos32(c.imag);
    c32 result = complex32(csh.sin * cs.cos, csh.cos * cs.sin);
    return result;
}

func c32
tanh32(c32 c)
{
    SinCos32 csh = sinhcosh32(2.0f * c.real);
    SinCos32 cs  = sincos32(2.0f * c.imag);
    f32 d = cs.cos + csh.cos;
    c32 result = complex32(csh.sin / d, cs.sin / d);
    return result;
}

func c32
acosh32(c32 c)
{
    c32 result = log32(c + sqrt32(c + 1.0f) * sqrt32(c - 1.0f));
    return result;
}

func c32
asinh32(c32 c)
{
    c32 result = complex32(0, -1.0f) * asin32(c * complex32(0, 1.0f));
    return result;
}

func c32
atanh32(c32 c)
{
    c32 result = complex32(0, -1.0f) * atan32(c * complex32(0, 1));
    return result;
}

//
//
// 64 bit
//
//

func c64
complex64(f64 real, f64 imag)
{
    c64 result;
    result.real = real;
    result.imag = imag;
    return result;
}

func c64
complex64_from_mag_phase(f64 magnitude, f64 phase)
{
    return complex64(magnitude * cos64(phase), magnitude * sin64(phase));
}

func c64
operator -(c64 a)
{
    c64 result;
    result.real = -a.real;
    result.imag = -a.imag;
    return result;
}

func c64 &
operator +=(c64 &a, c64 b)
{
    a.real += b.real;
    a.imag += b.imag;
    return a;
}

func c64
operator +(c64 a, c64 b)
{
    c64 result = a;
    result += b;
    return result;
}

func c64 &
operator -=(c64 &a, c64 b)
{
    a.real -= b.real;
    a.imag -= b.imag;
    return a;
}

func c64
operator -(c64 a, c64 b)
{
    c64 result = a;
    result -= b;
    return result;
}

func c64
operator *(c64 a, c64 b)
{
    c64 result;
    result.real = a.real * b.real - a.imag * b.imag;
    result.imag = a.real * b.imag + a.imag * b.real;
    return result;
}

func c64 &
operator *=(c64 &a, c64 b)
{
    a = a * b;
    return a;
}

func c64
operator /(c64 a, c64 b)
{
    f64 divisor = (b.real * b.real + b.imag * b.imag);
    c64 result;
    result.real = (a.real * b.real + a.imag * b.imag) / divisor;
    result.imag = (a.imag * b.real - a.real * b.imag) / divisor;
    return result;
}

func c64 &
operator /=(c64 &a, c64 b)
{
    a = a / b;
    return a;
}

func c64 &
operator +=(c64 &a, f64 b)
{
    a.real += b;
    return a;
}

func c64
operator +(c64 a, f64 b)
{
    c64 result = a;
    result.real += b;
    return result;
}

func c64
operator +(f64 a, c64 b)
{
    c64 result = b;
    result.real += a;
    return result;
}

func c64 &
operator -=(c64 &a, f64 b)
{
    a.real -= b;
    return a;
}

func c64
operator -(c64 a, f64 b)
{
    c64 result = a;
    result.real -= b;
    return result;
}

func c64
operator -(f64 a, c64 b)
{
    c64 result = b;
    result.real -= a;
    return result;
}

func c64 &
operator *=(c64 &a, f64 b)
{
    a.real *= b;
    a.imag *= b;
    return a;
}

func c64
operator *(c64 a, f64 b)
{
    c64 result = a;
    result *= b;
    return result;
}

func c64
operator *(f64 a, c64 b)
{
    c64 result = b;
    result *= a;
    return result;
}

func c64 &
operator /=(c64 &a, f64 b)
{
    a.real /= b;
    a.imag /= b;
    return a;
}

func c64
operator /(c64 a, f64 b)
{
    c64 result = a;
    result /= b;
    return result;
}

func c64
operator /(f64 a, c64 b)
{
    c64 result = complex64(a, 0.0);
    result /= b;
    return result;
}

func c64
conjugate64(c64 c)
{
    c64 result = c;
    result.imag = -result.imag;
    return result;
}

func f64
absolute64(c64 c)
{
    return hypot64(c.real, c.imag);
}

func f64
argument64(c64 c)
{
    return atan2_64(c.imag, c.real);
}

//
// NOTE(michiel): Complex elementary functions
//

func c64
sqrt64(c64 c)
{
    f64 x = c.real;
    f64 y = c.imag;

    if (y == 0.0)
    {
        if (x < 0.0) {
            return complex64(0, sqrt64(-x));
        } else if (x == 0.0) {
            return complex64(0, y);
        } else {
            return complex64(sqrt64(x), y);
        }
    }
    else if (x == 0.0)
    {
        f64 r = absolute64(y);
        r = sqrt64(0.5 * r);
        return complex64(r, (y > 0.0) ? r : -r);
    }

    f64 scale = 1.0;
    if ((absolute64(x) > 4.0) || (absolute64(y) > 4.0))
    {
        // NOTE(michiel): Rescale to avoid internal overflow
        x *= 0.25;
        y *= 0.25;
        scale = 2.0;
    }
    else
    {
        x *= 1.8014398509481984e16; // NOTE(michiel): 2^54
        y *= 1.8014398509481984e16;
        scale = 7.450580596923828125e-9; // NOTE(michiel): 2^-27
    }

    f64 r = absolute64(complex64(x, y));
    f64 t;
    if (x > 0.0)
    {
        t = sqrt64(0.5 * r + 0.5 * x);
        r = scale * absolute64((0.5 * y) / t);
        t *= scale;
    }
    else
    {
        r = sqrt64(0.5 * r - 0.5 * x);
        t = scale * absolute64((0.5 * y) / r);
        r *= scale;
    }

    c64 result = complex64(t, (y < 0.0) ? -r : r);
    return result;
}

func c64
exp64(c64 c)
{
    SinCos64 cs = sincos64(c.imag);
    f64 r = exp64(c.real);
    c64 result = complex64(r * cs.cos, r * cs.sin);
    return result;
}

func c64
log64(c64 c)
{
    f64 rr = absolute64(c);
    f64 p  = log64(rr);
    rr = argument64(c);
    c64 result = complex64(p, rr);
    return result;
}

func c64
log10_64(c64 c)
{
    f64 rr = absolute64(c);
    f64 p  = log10_64(rr);
    rr = argument64(c) * gInvLn10F64;
    c64 result = complex64(p, rr);
    return result;
}

func c64
pow64(c64 x, c64 y)
{

    f64 real = y.real;
    f64 imag = y.imag;

    f64 absX = absolute64(x);
    if (absX == 0.0) {
        return complex64(0, 0);
    }

    f64 argX = argument64(x);
    f64 r = pow64(absX, real);
    f64 theta = real * argX;

    if (imag != 0.0)
    {
        r = r * exp64(-imag * argX);
        theta = theta + imag * log64(absX);
    }

    SinCos64 th = sincos64(theta);
    c64 result = complex64(r * th.cos, r * th.sin);
    return result;
}

//
// NOTE(michiel): Complex trigonometric functions
//

func c64
cos64(c64 c)
{
    SinCos64 cs = sincos64(c.real);
    SinCos64 csh = sinhcosh64(c.imag);
    c64 result = complex64(cs.cos * csh.cos, -(cs.sin * csh.sin));
    return result;
}

func c64
sin64(c64 c)
{
    SinCos64 cs = sincos64(c.real);
    SinCos64 csh = sinhcosh64(c.imag);
    c64 result = complex64(cs.sin * csh.cos, cs.cos * csh.sin);
    return result;
}

func f64
tan_serie64(c64 c)
{
    // NOTE(michiel): Taylor-series expansion for cosh(2y) - cos(2x)
    f64 x = absolute64(2.0 * c.real);
    f64 y = absolute64(2.0 * c.imag);

    x = reduce_pi64(x);

    x = x * x;
    y = y * y;

    f64 x2 = 1.0;
    f64 y2 = 1.0;
    f64 f  = 1.0;
    f64 rn = 0.0;
    f64 d  = 0.0;

    // NOTE(michiel): While error is greater than machine precision
    f64 t;
    do
    {
        rn += 1.0;
        f  *= rn;
        rn += 1.0;
        f  *= rn;
        x2 *= x;
        y2 *= y;
        t   = y2 + x2;
        t  /= f;
        d  += f;

        rn += 1.0;
        f  *= rn;
        rn += 1.0;
        f  *= rn;
        x2 *= x;
        y2 *= y;
        t   = y2 - x2;
        t  /= f;
        d  += t;
    } while (absolute64(t / d) > 1.1e-16);
    return d;
}

func c64
tan64(c64 c)
{
    SinCos64 cs  = sincos64(2.0 * c.real);
    SinCos64 csh = sinhcosh64(2.0 * c.imag);
    f64 d = cs.cos + csh.cos;

    if (absolute64(d) < 0.25) {
        d = tan_serie64(c);
    }

    if (d == 0.0) {
        return complex64(gHugeF64, gHugeF64);
    }

    c64 result = complex64(cs.sin / d, csh.sin / d);
    return result;
}

func c64
asin64(c64 c)
{
    f64 x = c.real;
    f64 y = c.imag;

    c64 ct = c * gImaginary64;

    c64 zz = complex64((x - y) * (x + y), 2.0 * x * y);
    zz.real = 1.0 - zz.real;
    zz.imag = -zz.imag;

    c64 z2 = sqrt64(zz);

    zz = ct + z2;
    zz = log64(zz);

    c64 result = zz * (-gImaginary64);
    return result;
}

func c64
acos64(c64 c)
{
    c64 result = asin64(c);
    result = (gPiOver2F64 - result.real) - complex64(0.0, result.imag);
    return result;
}

func c64
atan64(c64 c)
{
    c64 result;
    f64 x = c.real;
    f64 y = c.imag;

    if ((x == 0.0) && (y > 1.0)) {
        return complex64(gHugeF64, gHugeF64);
    }

    f64 x2 = x * x;
    f64 a = 1.0 - x2 - (y * y);
    if (a == 0.0) {
        return complex64(gHugeF64, gHugeF64);
    }

    f64 t = 0.5 * atan2_64(2.0 * x, a);
    f64 w = reduce_pi64(t);

    t = y - 1.0;
    a = x2 + (t * t);
    if (a == 0.0) {
        return complex64(gHugeF64, gHugeF64);
    }

    t = y + 1.0;
    a = (x2 + (t * t)) / a;
    result = complex64(w, 0.25 * log64(a));
    return result;
}

//
// NOTE(michiel): Complex hyperbolic functions
//

func c64
cosh64(c64 c)
{
    SinCos64 csh = sinhcosh64(c.real);
    SinCos64 cs = sincos64(c.imag);
    c64 result = complex64(cs.cos * csh.cos, cs.sin * csh.sin);
    return result;
}

func c64
sinh64(c64 c)
{
    SinCos64 csh = sinhcosh64(c.real);
    SinCos64 cs = sincos64(c.imag);
    c64 result = complex64(csh.sin * cs.cos, csh.cos * cs.sin);
    return result;
}

func c64
tanh64(c64 c)
{
    SinCos64 csh = sinhcosh64(2.0 * c.real);
    SinCos64 cs  = sincos64(2.0 * c.imag);
    f64 d = cs.cos + csh.cos;
    c64 result = complex64(csh.sin / d, cs.sin / d);
    return result;
}

func c64
acosh64(c64 c)
{
    c64 result = log64(c + sqrt64(c + 1.0) * sqrt64(c - 1.0));
    return result;
}

func c64
asinh64(c64 c)
{
    c64 result = complex64(0, -1.0) * asin64(c * complex64(0, 1.0));
    return result;
}

func c64
atanh64(c64 c)
{
    c64 result = complex64(0, -1.0) * atan64(c * complex64(0, 1));
    return result;
}
