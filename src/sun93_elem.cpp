/*
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 */

/* ef_sqrtf.c -- float version of e_sqrt.c.
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 */

internal f32
sun93_sqrtf(f32 x)
{
	s32 ix = (s32)u32f32(x).u;
    u32 hx = ix & 0x7FFFFFFF;

    /* take care of Inf and NaN */
	if(!FLT_UWORD_IS_FINITE(hx)) {
	    return x * x + x;		// sqrt(NaN)=NaN, sqrt(+inf)=+inf, sqrt(-inf)=sNaN
    }
    /* take care of zero and -ves */
	if (FLT_UWORD_IS_ZERO(hx)) {
        return x; // sqrt(+-0) = +-0
    }
	if (ix < 0) {
        return (x - x) / (x - x);		/* sqrt(-ve) = sNaN */
    }

    /* normalize x */
	s32 m = (ix >> 23);
	if (FLT_UWORD_IS_SUBNORMAL(hx))
    {   /* subnormal x */
        s32 i;
	    for (i = 0; (ix & 0x00800000) == 0; ++i) {
            ix <<= 1;
        }
	    m -= i - 1;
	}
	m -= 127;	/* unbias exponent */
	ix = (ix & 0x007FFFFF) | 0x00800000;
	if (m & 1) {	/* odd m, double x to make it even */
	    ix += ix;
    }
	m >>= 1;	/* m = [m/2] */

    /* generate sqrt(x) bit by bit */
	ix += ix;
	s32 q = 0;        /* q = sqrt(x) */
    s32 s = 0;
	u32 r = 0x01000000;  /* r = moving bit from right to left */

	while (r != 0)
    {
	    s32 t = s + r;
	    if (t <= ix)
        {
            s    = t + r;
            ix  -= t;
            q   += r;
	    }
	    ix += ix;
	    r >>= 1;
	}

    /* use floating add to find out rounding direction */
	if (ix != 0)
    {
	    f32 z = 1.0f - gTinyF32; /* trigger inexact flag */
	    if (z >= 1.0f)
        {
	        z = 1.0f + gTinyF32;
            if (z > 1.0f) {
                q += 2;
            } else {
                q += (q & 1);
            }
	    }
	}
	ix = (q >> 1) + 0x3f000000L;
	ix += (m << 23);

    f32 result = u32f32((u32)ix).f;
	return result;
}

global const f32 gAtanHi[] = {
    4.6364760399e-01, /* atan(0.5)hi 0x3eed6338 */
    7.8539812565e-01, /* atan(1.0)hi 0x3f490fda */
    9.8279368877e-01, /* atan(1.5)hi 0x3f7b985e */
    1.5707962513e+00, /* atan(inf)hi 0x3fc90fda */
};

global const f32 gAtanLo[] = {
    5.0121582440e-09, /* atan(0.5)lo 0x31ac3769 */
    3.7748947079e-08, /* atan(1.0)lo 0x33222168 */
    3.4473217170e-08, /* atan(1.5)lo 0x33140fb4 */
    7.5497894159e-08, /* atan(inf)lo 0x33a22168 */
};

global const f32 gAtanF32[] = {
    3.3333334327e-01, /* 0x3eaaaaaa */
    -2.0000000298e-01, /* 0xbe4ccccd */
    1.4285714924e-01, /* 0x3e124925 */
    -1.1111110449e-01, /* 0xbde38e38 */
    9.0908870101e-02, /* 0x3dba2e6e */
    -7.6918758452e-02, /* 0xbd9d8795 */
    6.6610731184e-02, /* 0x3d886b35 */
    -5.8335702866e-02, /* 0xbd6ef16b */
    4.9768779427e-02, /* 0x3d4bda59 */
    -3.6531571299e-02, /* 0xbd15a221 */
    1.6285819933e-02, /* 0x3c8569d7 */
};

internal f32
sun93_atanf(f32 x)
{
	s32 ix = (s32)u32f32(x).u;
    u32 hx = ix & 0x7FFFFFFF;

    if (ix >= 0x50800000)
    {
        // NOTE(michiel): |x| >= 2^34
        if (FLT_UWORD_IS_NAN(ix)) {
            return x + x;
        } else if (hx > 0) {
            return gAtanHi[3] + gAtanLo[3];
        } else {
            return -gAtanHi[3] - gAtanLo[3];
        }
    }
    else
    {
        s32 id = 0;
        if (ix < 0x3EE00000) { // NOTE(michiel): |x| < 0.4375
            if (ix < 0x31000000) { // NOTE(michiel): |x| < 2^-29
                if ((gHugeF32 + x) > 1.0f) {
                    return x; // NOTE(michiel): Raise inexact
                }
            }
            id = -1;
        }
        else
        {
            x = sun93_fabsf(x);
            if (ix < 0x3F980000) { // NOTE(michiel): |x| < 1.1875
                if (ix < 0x3F300000) { // NOTE(michiel): 7/16 <= |x| < 11/16
                    id = 0;
                    x = (2.0f * x - 1.0f) / (2.0f * x);
                } else {
                    id = 1;
                    x = (x - 1.0f) / (x + 1.0f);
                }
            } else {
                if (ix < 0x401C0000) {
                    id = 2;
                    x = (x - 1.5f) / (1.0f + 1.5f * x);
                } else {
                    id = 3;
                    x = -1.0f / x;
                }
            }
        }

        f32 x2 = x * x;
        f32 x4 = x2 * x2;

        f32 s1 = gAtanF32[4] * x4 + gAtanF32[2];
        f32 s2 = gAtanF32[3] * x4 + gAtanF32[1];
        s1 = s1 * x4 + gAtanF32[0];
        s2 = s2 * x4;
        s1 = s1 * x2;

        f32 result;
        if (id < 0) {
            result = x - x * (s1 + s2);
        } else {
            result = gAtanHi[id] - ((x * (s1 + s2) - gAtanLo[id]) - x);
            result = (hx < 0) ? -result : result;
        }

        return result;
    }
}

#if 0
static const float
pi_o_4  = 7.8539818525e-01, /* 0x3f490fdb */
pi_o_2  = 1.5707963705e+00, /* 0x3fc90fdb */
pi      = 3.1415927410e+00,  /* 0x40490fdb */
pi_lo   = -8.7422776573e-08; /* 0xb3bbbd2e */

internal f32
sun93_atan2f(f32 y, f32 x)
{
	f32 z;
    s32 k,m,hx,hy,ix,iy;

	GET_FLOAT_WORD(hx,x);
	ix = hx&0x7fffffff;
	GET_FLOAT_WORD(hy,y);
	iy = hy&0x7fffffff;
	if(FLT_UWORD_IS_NAN(ix)||
	   FLT_UWORD_IS_NAN(iy))	/* x or y is NaN */
        return x+y;
	if(hx==0x3f800000) return atanf(y);   /* x=1.0 */
	m = ((hy>>31)&1)|((hx>>30)&2);	/* 2*sign(x)+sign(y) */

    /* when y = 0 */
	if(FLT_UWORD_IS_ZERO(iy)) {
	    switch(m) {
            case 0:
            case 1: return y; 	/* atan(+-0,+anything)=+-0 */
            case 2: return  pi+tiny;/* atan(+0,-anything) = pi */
            case 3: return -pi-tiny;/* atan(-0,-anything) =-pi */
	    }
	}
    /* when x = 0 */
	if(FLT_UWORD_IS_ZERO(ix)) return (hy<0)?  -pi_o_2-tiny: pi_o_2+tiny;

    /* when x is INF */
	if(FLT_UWORD_IS_INFINITE(ix)) {
	    if(FLT_UWORD_IS_INFINITE(iy)) {
            switch(m) {
                case 0: return  pi_o_4+tiny;/* atan(+INF,+INF) */
                case 1: return -pi_o_4-tiny;/* atan(-INF,+INF) */
                case 2: return  (float)3.0*pi_o_4+tiny;/*atan(+INF,-INF)*/
                case 3: return (float)-3.0*pi_o_4-tiny;/*atan(-INF,-INF)*/
            }
	    } else {
            switch(m) {
                case 0: return  zero  ;	/* atan(+...,+INF) */
                case 1: return -zero  ;	/* atan(-...,+INF) */
                case 2: return  pi+tiny  ;	/* atan(+...,-INF) */
                case 3: return -pi-tiny  ;	/* atan(-...,-INF) */
            }
	    }
	}
    /* when y is INF */
	if(FLT_UWORD_IS_INFINITE(iy)) return (hy<0)? -pi_o_2-tiny: pi_o_2+tiny;

    /* compute y/x */
	k = (iy-ix)>>23;
	if(k > 60) z=pi_o_2+(float)0.5*pi_lo; 	/* |y/x| >  2**60 */
	else if(hx<0&&k<-60) z=0.0; 	/* |y|/x < -2**60 */
	else z=atanf(fabsf(y/x));	/* safe to do y/x */
	switch (m) {
	    case 0: return       z  ;	/* atan(+,+) */
	    case 1: {
            __uint32_t zh;
            GET_FLOAT_WORD(zh,z);
            SET_FLOAT_WORD(z,zh ^ 0x80000000);
        }
        return       z  ;	/* atan(-,+) */
	    case 2: return  pi-(z-pi_lo);/* atan(+,-) */
	    default: /* case 3 */
        return  (z-pi_lo)-pi;/* atan(-,-) */
	}
}
#endif
