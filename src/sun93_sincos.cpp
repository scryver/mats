/* kf_cos.c -- float version of k_cos.c
 * Conversion to float by Ian Lance Taylor, Cygnus Support, ian@cygnus.com.
 */

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

global const f32 gCos1F32       =  4.1666667908e-02f; /* 0x3d2aaaab */
global const f32 gCos2F32       = -1.3888889225e-03f; /* 0xbab60b61 */
global const f32 gCos3F32       =  2.4801587642e-05f; /* 0x37d00d01 */
global const f32 gCos4F32       = -2.7557314297e-07f; /* 0xb493f27c */
global const f32 gCos5F32       =  2.0875723372e-09f; /* 0x310f74f6 */
global const f32 gCos6F32       = -1.1359647598e-11f; /* 0xad47d74e */

global const f32 gSin1F32       = -1.6666667163e-01f; /* 0xbe2aaaab */
global const f32 gSin2F32       =  8.3333337680e-03f; /* 0x3c088889 */
global const f32 gSin3F32       = -1.9841270114e-04f; /* 0xb9500d01 */
global const f32 gSin4F32       =  2.7557314297e-06f; /* 0x3638ef1b */
global const f32 gSin5F32       = -2.5050759689e-08f; /* 0xb2d72f34 */
global const f32 gSin6F32       =  1.5896910177e-10f; /* 0x2f2ec9d3 */

global const s32 gInitJK[] = {4,7,9}; /* initial value for jk */

global const f32 gPiOver2[] = {
    1.5703125000e+00f, /* 0x3fc90000 */
    4.5776367188e-04f, /* 0x39f00000 */
    2.5987625122e-05f, /* 0x37da0000 */
    7.5437128544e-08f, /* 0x33a20000 */
    6.0026650317e-11f, /* 0x2e840000 */
    7.3896444519e-13f, /* 0x2b500000 */
    5.3845816694e-15f, /* 0x27c20000 */
    5.6378512969e-18f, /* 0x22d00000 */
    8.3009228831e-20f, /* 0x1fc40000 */
    3.2756352257e-22f, /* 0x1bc60000 */
    6.3331015649e-25f, /* 0x17440000 */
};

/*
 * Table of constants for 2/pi, 396 Hex digits (476 decimal) of 2/pi
 */
global const s32 g2OverPiArr[] = {
    0xA2, 0xF9, 0x83, 0x6E, 0x4E, 0x44, 0x15, 0x29, 0xFC,
    0x27, 0x57, 0xD1, 0xF5, 0x34, 0xDD, 0xC0, 0xDB, 0x62,
    0x95, 0x99, 0x3C, 0x43, 0x90, 0x41, 0xFE, 0x51, 0x63,
    0xAB, 0xDE, 0xBB, 0xC5, 0x61, 0xB7, 0x24, 0x6E, 0x3A,
    0x42, 0x4D, 0xD2, 0xE0, 0x06, 0x49, 0x2E, 0xEA, 0x09,
    0xD1, 0x92, 0x1C, 0xFE, 0x1D, 0xEB, 0x1C, 0xB1, 0x29,
    0xA7, 0x3E, 0xE8, 0x82, 0x35, 0xF5, 0x2E, 0xBB, 0x44,
    0x84, 0xE9, 0x9C, 0x70, 0x26, 0xB4, 0x5F, 0x7E, 0x41,
    0x39, 0x91, 0xD6, 0x39, 0x83, 0x53, 0x39, 0xF4, 0x9C,
    0x84, 0x5F, 0x8B, 0xBD, 0xF9, 0x28, 0x3B, 0x1F, 0xF8,
    0x97, 0xFF, 0xDE, 0x05, 0x98, 0x0F, 0xEF, 0x2F, 0x11,
    0x8B, 0x5A, 0x0A, 0x6D, 0x1F, 0x6D, 0x36, 0x7E, 0xCF,
    0x27, 0xCB, 0x09, 0xB7, 0x4F, 0x46, 0x3F, 0x66, 0x9E,
    0x5F, 0xEA, 0x2D, 0x75, 0x27, 0xBA, 0xC7, 0xEB, 0xE5,
    0xF1, 0x7B, 0x3D, 0x07, 0x39, 0xF7, 0x8A, 0x52, 0x92,
    0xEA, 0x6B, 0xFB, 0x5F, 0xB1, 0x1F, 0x8D, 0x5D, 0x08,
    0x56, 0x03, 0x30, 0x46, 0xFC, 0x7B, 0x6B, 0xAB, 0xF0,
    0xCF, 0xBC, 0x20, 0x9A, 0xF4, 0x36, 0x1D, 0xA9, 0xE3,
    0x91, 0x61, 0x5E, 0xE6, 0x1B, 0x08, 0x65, 0x99, 0x85,
    0x5F, 0x14, 0xA0, 0x68, 0x40, 0x8D, 0xFF, 0xD8, 0x80,
    0x4D, 0x73, 0x27, 0x31, 0x06, 0x06, 0x15, 0x56, 0xCA,
    0x73, 0xA8, 0xC9, 0x60, 0xE2, 0x7B, 0xC0, 0x8C, 0x6B,
};

/* This array is like the one in e_rem_pio2.c, but the numbers are
   single precision and the last 8 bits are forced to 0.  */
global const s32 gNPiOver2_hw[] = {
    0x3fc90f00, 0x40490f00, 0x4096cb00, 0x40c90f00, 0x40fb5300, 0x4116cb00,
    0x412fed00, 0x41490f00, 0x41623100, 0x417b5300, 0x418a3a00, 0x4196cb00,
    0x41a35c00, 0x41afed00, 0x41bc7e00, 0x41c90f00, 0x41d5a000, 0x41e23100,
    0x41eec200, 0x41fb5300, 0x4203f200, 0x420a3a00, 0x42108300, 0x4216cb00,
    0x421d1400, 0x42235c00, 0x4229a500, 0x422fed00, 0x42363600, 0x423c7e00,
    0x4242c700, 0x42490f00
};

internal f32
sun93_fabsf(f32 x)
{
	u32 ix = u32f32(x).u;
    f32 result = u32f32(ix & 0x7FFFFFFF).f;
    return result;
}

internal f32
sun93_copysignf(f32 x, f32 y)
{
	u32 ix = u32f32(x).u;
    u32 iy = u32f32(y).u;
    f32 result = u32f32((ix & 0x7FFFFFFF) | (iy & 0x80000000)).f;
    return result;
}

internal f32 sun93_floorf(f32 x)
{
	s32 i0 = (s32)u32f32(x).u;
	u32 ix = (i0 & 0x7fffffff);
	s32 j0 = (ix >> 23) - 0x7f;
	if (j0 < 23)
    {
	    if (j0 < 0)
        {   /* raise inexact if x != 0 */
            if (gHugeF32 + x > 0.0f)
            {
                /* return 0*sign(x) if |x|<1 */
                if (i0 >= 0) {
                    i0 = 0;
                } else if (!FLT_UWORD_IS_ZERO(ix)) {
                    i0 = 0xbf800000;
                }
            }
	    }
        else
        {
            u32 i = (0x007fffff) >> j0;
            if ((i0 & i) == 0) {
                return x; /* x is integral */
            }
            if (gHugeF32 + x > 0.0f)
            {   /* raise inexact flag */
                if (i0 < 0) {
                    i0 += (0x00800000) >> j0;
                }
                i0 &= (~i);
            }
	    }
	}
    else
    {
	    if (!FLT_UWORD_IS_FINITE(ix)) {
            return x+x;	/* inf or NaN */
        } else {
            return x;	  /* x is integral */
        }
	}
    x = u32f32((u32)i0).f;
	return x;
}

internal f32
sun93_scalbnf(f32 x, int n)
{
    s32 ix = (s32)u32f32(x).u;
    u32 hx = ix & 0x7FFFFFFF;
    s32 k  = hx >> 23;		/* extract exponent */
	if (FLT_UWORD_IS_ZERO(hx)) {
	    return x;
    }
    if (!FLT_UWORD_IS_FINITE(hx)) {
	    return x+x;		/* NaN or Inf */
    }
    if (FLT_UWORD_IS_SUBNORMAL(hx))
    {
	    x *= g2pow25F32;
        ix = (s32)u32f32(x).u;
	    k = ((ix & 0x7f800000) >> 23) - 25;
        if (n < UNDERFLOW_INT) {
            return gTinyF32 * x; 	/*underflow*/
        }
    }
    k = k + n;
    if (k > FLT_LARGEST_EXP) {
        return gHugeF32 * sun93_copysignf(gHugeF32, x); /* overflow  */
    }
    if (k > 0) 				/* normal result */
    {
        x = u32f32((ix & 0x807fffff) | (k << 23)).f;
        return x;
    }
    if (k < FLT_SMALLEST_EXP) {
        if (n > OVERFLOW_INT) { 	/* in case integer overflow in n+k */
            return gHugeF32 * sun93_copysignf(gHugeF32, x);	/*overflow*/
        } else {
            return gTinyF32 * sun93_copysignf(gTinyF32, x);	/*underflow*/
        }
    }
    k += 25;				/* subnormal result */
    x = u32f32((ix & 0x807fffff) | (k << 23)).f;
    return x * g2powMin25F32;
}

internal s32
__kernel_rem_pio2f(f32 *x, f32 *y, s32 e0, s32 nx, s32 prec, const s32 *ipio2)
{
    s32 n,iq[20],k,ih;
    f32 z,fw,f[20],fq[20],q[20];

    /* initialize jk*/
    s32 jk = gInitJK[prec];
    s32 jp = jk;

    /* determine jx,jv,q0, note that 3>q0 */
    s32 jx = nx - 1;
    s32 jv = (e0 - 3) / 8;
    if (jv < 0) {
        jv = 0;
    }
    s32 q0 = e0 - 8 * (jv + 1);

    /* set up f[0] to f[jx+jk] where f[jx+jk] = ipio2[jv+jk] */
    s32 j = jv - jx;
    s32 m = jx + jk;
    for(s32 i = 0; i <= m; i++, j++) {
        f[i] = (j < 0) ? gZeroF32 : (f32)ipio2[j];
    }

    /* compute q[0],q[1],...q[jk] */
    for (s32 i = 0; i <= jk; i++) {
        for (j = 0, fw=0.0f; j <= jx; j++) {
            fw += x[j] * f[jx+i-j];
        }
        q[i] = fw;
    }

    s32 jz = jk;
    recompute:
    /* distill q[] into iq[] reversingly */
    s32 loopI = 0;
    for (j = jz, z = q[jz]; j > 0; loopI++, j--)
    {
        fw    =  (f32)((s32)(g2powMin8F32 * z));
        iq[loopI] =  (s32)(z - g2pow8F32 * fw);
        z     =  q[j - 1] + fw;
    }

    /* compute n */
    z  = sun93_scalbnf(z, (s32)q0);	/* actual value of z */
    z -= (f32)8.0 * sun93_floorf(z * (f32)0.125f);	/* trim off integer >= 8 */
    n  = (s32)z;
    z -= (f32)n;
    ih = 0;
    if (q0 > 0) {	/* need iq[jz-1] to determine n */
        s32 i  = (iq[jz-1] >> (8 - q0));
        n += i;
        iq[jz-1] -= i << (8 - q0);
        ih = iq[jz-1] >> (7 - q0);
    } else if (q0 == 0) {
        ih = iq[jz-1] >> 7;
    } else if(z >= 0.5f) {
        ih=2;
    }

    if (ih > 0)
    {	/* q > 0.5 */
        n += 1;
        s32 carry = 0;
        for (s32 i = 0; i < jz; i++)
        {   /* compute 1-q */
            j = iq[i];
            if (carry == 0)
            {
                if (j != 0) {
                    carry = 1;
                    iq[i] = 0x100 - j;
                }
            } else {
                iq[i] = 0xff - j;
            }
        }
        if (q0 > 0) {		/* rare case: chance is 1 in 12 */
            switch (q0) {
                case 1: { iq[jz-1] &= 0x7f; } break;
                case 2: { iq[jz-1] &= 0x3f; } break;
            }
        }
        if (ih == 2) {
            z = gOneF32 - z;
            if (carry != 0) {
                z -= sun93_scalbnf(gOneF32,(s32)q0);
            }
        }
    }

    /* check if recomputation is needed */
    if (z == gZeroF32)
    {
        s32 jj = 0;
        for (s32 i = jz - 1; i >= jk; i--) {
            jj |= iq[i];
        }
        if (jj == 0)
        {   /* need recomputation */
            for (k = 1; iq[jk-k] == 0; k++);   /* k = no. of terms needed */

            for (s32 i= jz + 1; i <= jz + k; i++)
            {   /* add q[jz+1] to q[jz+k] */
                f[jx+i] = (f32)ipio2[jv+i];
                for (s32 ij = 0, fw=0.0f; ij <= jx; ij++) {
                    fw += x[ij] * f[jx+i-ij];
                }
                q[i] = fw;
            }
            jz += k;
            goto recompute;
        }
    }

    /* chop off zero terms */
    if (z == 0.0f)
    {
        jz -= 1;
        q0 -= 8;
        while (iq[jz] == 0) {
            jz--;
            q0-=8;
        }
    }
    else
    { /* break z into 8-bit if necessary */
        z = sun93_scalbnf(z, -(s32)q0);
        if (z >= g2pow8F32)
        {
            fw = (f32)((s32)(g2powMin8F32 * z));
            iq[jz] = (s32)(z-g2pow8F32 * fw);
            jz += 1;
            q0 += 8;
            iq[jz] = (s32)fw;
        }
        else
        {
            iq[jz] = (s32)z;
        }
    }

    /* convert integer "bit" chunk to floating-point value */
    fw = sun93_scalbnf(gOneF32, (s32)q0);
    for (s32 i = jz; i >= 0; i--) {
        q[i] = fw * (f32)iq[i];
        fw *= g2powMin8F32;
    }

    /* compute gPiOver2[0,...,jp]*q[jz,...,0] */
    for (s32 i = jz; i >= 0; i--) {
        for (fw=0.0f, k = 0; k <= jp && k <= jz - i; k++) {
            fw += gPiOver2[k] * q[i + k];
        }
        fq[jz - i] = fw;
    }

    /* compress fq[] into y[] */
    switch(prec)
    {
        case 0: {
            fw = 0.0;
            for (s32 i = jz; i >= 0; i--) {
                fw += fq[i];
            }
            y[0] = (ih == 0) ? fw : -fw;
        } break;
        case 1:
        case 2: {
            fw = 0.0;
            for (s32 i = jz; i >= 0; i--) {
                fw += fq[i];
            }
            y[0] = (ih == 0) ? fw: -fw;
            fw = fq[0] - fw;
            for (s32 i = 1; i <= jz; i++) {
                fw += fq[i];
            }
            y[1] = (ih == 0) ? fw: -fw;
        } break;
        case 3: {	/* painful */
            for (s32 i = jz; i > 0; i--) {
                fw        = fq[i - 1] + fq[i];
                fq[i]    += fq[i - 1] - fw;
                fq[i - 1] = fw;
            }
            for (s32 i = jz; i > 1; i--) {
                fw        = fq[i - 1] + fq[i];
                fq[i]    += fq[i - 1] - fw;
                fq[i - 1] = fw;
            }
            fw = 0.0f;
            for (s32 i = jz; i >= 2; i--) {
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

internal f32
__kernel_cosf(f32 x, f32 y)
{
	u32 ix = u32f32(x).u;
	ix &= 0x7fffffff;			/* ix = |x|'s high word*/
    /* if x < 2**27 */
	if (ix < 0x32000000)
    {
	    if((s32)x == 0) {
            return gOneF32;		/* generate inexact */
        }
	}
	f32 z  = x * x;
	f32 r  = z * (gCos1F32 + z * (gCos2F32 + z * (gCos3F32 + z * (gCos4F32 + z * (gCos5F32 + z * gCos6F32)))));
    /* if |x| < 0.3 */
	if (ix < 0x3e99999a)
    {
	    return gOneF32 - (0.5f * z - (z*r - x*y));
    }
	else
    {
        f32 qx;
        /* x > 0.78125 */
	    if (ix > 0x3f480000)
        {
            qx = 0.28125f;
	    } else {
            /* x/4 */
	        qx = u32f32(ix - 0x01000000).f;
	    }
	    f32 hz = 0.5f*z - qx;
	    f32 a  = gOneF32 - qx;
	    return a - (hz - (z * r - x * y));
	}
}

internal f32
__kernel_sinf(f32 x, f32 y, b32 isYNonZero)
{
	u32 ix = u32f32(x).u;
	ix &= 0x7fffffff;

    /* |x| < 2**-27 */
	if (ix < 0x32000000)
    {
        if((s32)x == 0) {
            return x;		/* generate inexact */
        }
    }
	f32 z = x*x;
	f32 v = z*x;
	f32 r = gSin2F32 + z * (gSin3F32 + z * (gSin4F32 + z * (gSin5F32 + z * gSin6F32)));
	if (isYNonZero == 0) {
        return x + v * (gSin1F32 + z * r);
    } else {
        return x - ((z * (gHalfF32 * y - v * r) - y) - v * gSin1F32);
    }
}

/* __ieee754_rem_pio2f(x,y)
 *
 * return the remainder of x rem pi/2 in y[0]+y[1]
 */

internal s32
__ieee754_rem_pio2f(f32 x, f32 *y)
{
    f32 z,w,t,r,fn;
    f32 tx[3];
    s32 i,j,n;
    s32 e0, nx;

	s32 hx = (s32)u32f32(x).u;
	s32 ix = hx & 0x7fffffff;

    if (ix <= 0x3f490fd8)   /* |x| ~<= pi/4 , no need for reduction */
    {
        y[0] = x;
        y[1] = 0;
        return 0;
    }
    if (ix < 0x4016cbe4)
    {  /* |x| < 3pi/4, special case with n=+-1 */
        if (hx > 0)
        {
            z = x - gPiOver2F32_1;
            if ((ix & 0xfffffff0) != 0x3fc90fd0)
            { /* 24+24 bit pi OK */
                y[0] = z - gPiOver2F32_1t;
                y[1] = (z - y[0]) - gPiOver2F32_1t;
            }
            else
            {		/* near pi/2, use 24+24+24 bit pi */
                z -= gPiOver2F32_2;
                y[0] = z - gPiOver2F32_2t;
                y[1] = (z - y[0]) - gPiOver2F32_2t;
            }
            return 1;
        }
        else
        {	/* negative x */
            z = x + gPiOver2F32_1;
            if ((ix & 0xfffffff0) != 0x3fc90fd0)
            { /* 24+24 bit pi OK */
                y[0] = z + gPiOver2F32_1t;
                y[1] = (z - y[0]) + gPiOver2F32_1t;
            }
            else
            {		/* near pi/2, use 24+24+24 bit pi */
                z += gPiOver2F32_2;
                y[0] = z + gPiOver2F32_2t;
                y[1] = (z - y[0]) + gPiOver2F32_2t;
            }
            return -1;
        }
    }
    if (ix <= 0x43490f80)
    { /* |x| ~<= 2^7*(pi/2), medium size */
        t  = sun93_fabsf(x);
        n  = (s32)(t * g2OverPiF32 + gHalfF32);
        fn = (f32)n;
        r  = t - fn * gPiOver2F32_1;
        w  = fn * gPiOver2F32_1t;	/* 1st round good to 40 bit */
        if (n < 32 && (s32)(ix & 0xffffff00) != gNPiOver2_hw[n - 1])
        {
            y[0] = r - w;	/* quick check no cancellation */
        }
        else
        {
            j  = ix >> 23;
            y[0] = r - w;
            u32 high = u32f32(y[0]).u;
            i = j - ((high >> 23) & 0xff);
            if (i > 8)
            {  /* 2nd iteration needed, good to 57 */
                t  = r;
                w  = fn * gPiOver2F32_2;
                r  = t - w;
                w  = fn * gPiOver2F32_2t - ((t - r) - w);
                y[0] = r - w;
                high = u32f32(y[0]).u;
                i = j - ((high >> 23) & 0xff);
                if (i > 25)
                {	/* 3rd iteration need, 74 bits acc */
                    t  = r;	/* will cover all possible cases */
                    w  = fn * gPiOver2F32_3;
                    r  = t - w;
                    w  = fn * gPiOver2F32_3t - ((t - r) - w);
                    y[0] = r - w;
                }
            }
        }
        y[1] = (r - y[0]) - w;
        if (hx < 0) {
            y[0] = -y[0];
            y[1] = -y[1];
            return -n;
        } else {
            return n;
        }
    }
    /*
     * all other (large) arguments
     */
    if (!FLT_UWORD_IS_FINITE(ix)) {
        y[0] = y[1] = x-x;
        return 0;
    }
    /* set z = scalbn(|x|,ilogb(x)-7) */
    e0 	= (s32)((ix >> 23) - 134);	/* e0 = ilogb(z)-7; */
    z = u32f32((u32)(ix - ((s32)e0 << 23))).f;
    for (i = 0; i < 2; i++) {
        tx[i] = (f32)((s32)(z));
        z     = (z - tx[i]) * g2pow8F32;
    }
    tx[2] = z;
    nx = 3;
    while (tx[nx - 1] == gZeroF32) {
        nx--;	/* skip zero term */
    }
    n  =  __kernel_rem_pio2f(tx, y, e0, nx, 2, g2OverPiArr);
    if (hx < 0) {
        y[0] = -y[0];
        y[1] = -y[1];
        return -n;
    }
    return n;
}

internal f32
sun93_cosf(f32 x)
{
	f32 y[2], z = 0.0f;
	s32 n;

	s32 ix = (s32)u32f32(x).u;

    /* |x| ~< pi/4 */
	ix &= 0x7fffffff;
	if (ix <= 0x3f490fd8) {
        return __kernel_cosf(x, z);
    } else if (!FLT_UWORD_IS_FINITE(ix)) {
        /* cos(Inf or NaN) is NaN */
        return x-x;
    } else {
        /* argument reduction needed */
	    n = __ieee754_rem_pio2f(x,y);
	    switch(n & 3)
        {
            case 0:  return  __kernel_cosf(y[0], y[1]);
            case 1:  return -__kernel_sinf(y[0], y[1], 1);
            case 2:  return -__kernel_cosf(y[0], y[1]);
            default: return  __kernel_sinf(y[0], y[1], 1);
	    }
	}
}

internal f32
sun93_sinf(f32 x)
{
	f32 y[2], z = 0.0f;

	s32 ix = (s32)u32f32(x).u;

    /* |x| ~< pi/4 */
	ix &= 0x7fffffff;
	if (ix <= 0x3f490fd8) {
        return __kernel_sinf(x, z, 0);
    } else if (!FLT_UWORD_IS_FINITE(ix)) {
        /* sin(Inf or NaN) is NaN */
        return x-x;
    } else {
        /* argument reduction needed */
	    s32 n = __ieee754_rem_pio2f(x,y);
	    switch (n & 3)
        {
            case 0:  return  __kernel_sinf(y[0], y[1], 1);
            case 1:  return  __kernel_cosf(y[0], y[1]);
            case 2:  return -__kernel_sinf(y[0], y[1], 1);
            default: return -__kernel_cosf(y[0], y[1]);
	    }
	}
}
