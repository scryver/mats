
internal void
fft_normal2(u32 count, c32 *signal, c32 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < halfCount; index += 2)
    {
        u32 index0 = index + 0;
        u32 index1 = index + 1;
        u32 index2 = index0 + halfCount;
        u32 index3 = index1 + halfCount;

        u32 revIndex0 = reverse_bits32(index0, highBit.index);
        u32 revIndex1 = revIndex0 ^ halfCount;
        u32 revIndex2 = revIndex0 + 1;
        u32 revIndex3 = revIndex1 + 1;
        dest[index0] = signal[revIndex0];
        dest[index1] = signal[revIndex1];
        dest[index2] = signal[revIndex2];
        dest[index3] = signal[revIndex3];
    }

    // NOTE(michiel): w = e^(-i*2pi*j/m)
    SinCos32 cs4 = sincos32(-0.5f * F32_PI);
    c32 w4 = complex32(cs4.cos, cs4.sin);
    SinCos32 cs81 = sincos32(-0.25f * F32_PI);
    c32 w81 = complex32(cs81.cos, cs81.sin);
    SinCos32 cs82 = cs4;
    c32 w82 = complex32(cs82.cos, cs82.sin);
    SinCos32 cs83 = sincos32(-0.75f * F32_PI);
    c32 w83 = complex32(cs83.cos, cs83.sin);

    for (u32 k = 0; k < count; k += 8)
    {
        c32 *src0 = dest + k + 0;
        c32 *src1 = dest + k + 1;
        c32 *src2 = dest + k + 2;
        c32 *src3 = dest + k + 3;
        c32 *src4 = dest + k + 4;
        c32 *src5 = dest + k + 5;
        c32 *src6 = dest + k + 6;
        c32 *src7 = dest + k + 7;

        // NOTE(michiel): 2 -> 2
        c32 E20 = *src0;
        c32 O20 = *src1;
        *src0 = E20 + O20;
        *src1 = E20 - O20;
        c32 E21 = *src2;
        c32 O21 = *src3;
        *src2 = E21 + O21;
        *src3 = E21 - O21;
        c32 E22 = *src4;
        c32 O22 = *src5;
        *src4 = E22 + O22;
        *src5 = E22 - O22;
        c32 E23 = *src6;
        c32 O23 = *src7;
        *src6 = E23 + O23;
        *src7 = E23 - O23;

        // NOTE(michiel): 4 -> 4
        c32 E40 = *src0;
        c32 O40 = *src2;
        *src0 = E40 + O40;
        *src2 = E40 - O40;

        c32 E41 = *src1;
        c32 O41 = w4 * *src3;
        *src1 = E41 + O41;
        *src3 = E41 - O41;

        c32 E42 = *src4;
        c32 O42 = *src6;
        *src4 = E42 + O42;
        *src6 = E42 - O42;

        c32 E43 = *src5;
        c32 O43 = w4 * *src7;
        *src5 = E43 + O43;
        *src7 = E43 - O43;

        // NOTE(michiel): 8 -> 8
        c32 E80 = *src0;
        c32 O80 = *src4;
        *src0 = E80 + O80;
        *src4 = E80 - O80;

        c32 E81 = *src1;
        c32 O81 = w81 * *src5;
        *src1 = E81 + O81;
        *src5 = E81 - O81;

        c32 E82 = *src2;
        c32 O82 = w82 * *src6;
        *src2 = E82 + O82;
        *src6 = E82 - O82;

        c32 E83 = *src3;
        c32 O83 = w83 * *src7;
        *src3 = E83 + O83;
        *src7 = E83 - O83;
    }

    u32 halfM = 8;
    u32 m = 16;
    while (m <= count)
    {
        f32 oneOverM = (2.0f * F32_PI) / (f32)m;

        for (u32 k = 0; k < count; k += m)
        {
            c32 *src0 = dest + k;
            c32 *src1 = dest + k + halfM;
            for (u32 j = 0; j < halfM; j += 8)
            {
                // NOTE(michiel): w = e^(-i*2pi*j/m)
                SinCos32 cs0 = sincos32(-(f32)j * oneOverM);
                c32 w0 = complex32(cs0.cos, cs0.sin);
                SinCos32 cs1 = sincos32(-(f32)(j + 1) * oneOverM);
                c32 w1 = complex32(cs1.cos, cs1.sin);
                SinCos32 cs2 = sincos32(-(f32)(j + 2) * oneOverM);
                c32 w2 = complex32(cs2.cos, cs2.sin);
                SinCos32 cs3 = sincos32(-(f32)(j + 3) * oneOverM);
                c32 w3 = complex32(cs3.cos, cs3.sin);
                SinCos32 cs4 = sincos32(-(f32)(j + 4) * oneOverM);
                c32 w4 = complex32(cs4.cos, cs4.sin);
                SinCos32 cs5 = sincos32(-(f32)(j + 5) * oneOverM);
                c32 w5 = complex32(cs5.cos, cs5.sin);
                SinCos32 cs6 = sincos32(-(f32)(j + 6) * oneOverM);
                c32 w6 = complex32(cs6.cos, cs6.sin);
                SinCos32 cs7 = sincos32(-(f32)(j + 7) * oneOverM);
                c32 w7 = complex32(cs7.cos, cs7.sin);
                c32 E0 = *src0;
                c32 O0 = w0 * *src1;
                *src0++ = E0 + O0;
                *src1++ = E0 - O0;
                c32 E1 = *src0;
                c32 O1 = w1 * *src1;
                *src0++ = E1 + O1;
                *src1++ = E1 - O1;
                c32 E2 = *src0;
                c32 O2 = w2 * *src1;
                *src0++ = E2 + O2;
                *src1++ = E2 - O2;
                c32 E3 = *src0;
                c32 O3 = w3 * *src1;
                *src0++ = E3 + O3;
                *src1++ = E3 - O3;
                c32 E4 = *src0;
                c32 O4 = w4 * *src1;
                *src0++ = E4 + O4;
                *src1++ = E4 - O4;
                c32 E5 = *src0;
                c32 O5 = w5 * *src1;
                *src0++ = E5 + O5;
                *src1++ = E5 - O5;
                c32 E6 = *src0;
                c32 O6 = w6 * *src1;
                *src0++ = E6 + O6;
                *src1++ = E6 - O6;
                c32 E7 = *src0;
                c32 O7 = w7 * *src1;
                *src0++ = E7 + O7;
                *src1++ = E7 - O7;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_normal3(u32 count, c32 *signal, c32 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < halfCount; index += 2)
    {
        u32 index0 = index + 0;
        u32 index1 = index + 1;
        u32 index2 = index0 + halfCount;
        u32 index3 = index1 + halfCount;

        u32 revIndex0 = reverse_bits32(index0, highBit.index);
        u32 revIndex1 = revIndex0 ^ halfCount;
        u32 revIndex2 = revIndex0 + 1;
        u32 revIndex3 = revIndex1 + 1;
        dest[index0] = signal[revIndex0];
        dest[index1] = signal[revIndex1];
        dest[index2] = signal[revIndex2];
        dest[index3] = signal[revIndex3];
    }

    // NOTE(michiel): w = e^(-i*2pi*j/m)
    SinCos32_4x csBase = sincos32_4x(F32_4x(0.0f, -0.25f * F32_PI, -0.5f * F32_PI, -0.75f * F32_PI));
    f32_4x wReal = csBase.cos;
    f32_4x wImag = csBase.sin;

    for (u32 k = 0; k < count; k += 8)
    {
        c32 src0 = *(dest + k + 0);
        c32 src1 = *(dest + k + 1);
        c32 src2 = *(dest + k + 2);
        c32 src3 = *(dest + k + 3);
        c32 src4 = *(dest + k + 4);
        c32 src5 = *(dest + k + 5);
        c32 src6 = *(dest + k + 6);
        c32 src7 = *(dest + k + 7);

        // NOTE(michiel): 2 -> 2
        f32 E20r = src0.real;
        f32 E20i = src0.imag;
        f32 O20r = src1.real;
        f32 O20i = src1.imag;

        f32 E21r = src2.real;
        f32 E21i = src2.imag;
        f32 O21r = src3.real;
        f32 O21i = src3.imag;

        f32 E22r = src4.real;
        f32 E22i = src4.imag;
        f32 O22r = src5.real;
        f32 O22i = src5.imag;

        f32 E23r = src6.real;
        f32 E23i = src6.imag;
        f32 O23r = src7.real;
        f32 O23i = src7.imag;

        src0.real = E20r + O20r;
        src0.imag = E20i + O20i;
        src1.real = E20r - O20r;
        src1.imag = E20i - O20i;

        src2.real = E21r + O21r;
        src2.imag = E21i + O21i;
        src3.real = E21r - O21r;
        src3.imag = E21i - O21i;

        src4.real = E22r + O22r;
        src4.imag = E22i + O22i;
        src5.real = E22r - O22r;
        src5.imag = E22i - O22i;

        src6.real = E23r + O23r;
        src6.imag = E23i + O23i;
        src7.real = E23r - O23r;
        src7.imag = E23i - O23i;

        // NOTE(michiel): 4 -> 4
        c32 E40 = src0;
        c32 E41 = src1;
        c32 O40 = src2;
        c32 O41 = complex32(wReal.e[2], wImag.e[2]) * src3;
        c32 E42 = src4;
        c32 E43 = src5;
        c32 O42 = src6;
        c32 O43 = complex32(wReal.e[2], wImag.e[2]) * src7;

        src0 = E40 + O40;
        src1 = E41 + O41;
        src2 = E40 - O40;
        src3 = E41 - O41;
        src4 = E42 + O42;
        src5 = E43 + O43;
        src6 = E42 - O42;
        src7 = E43 - O43;

        // NOTE(michiel): 8 -> 8
        c32 E80 = src0;
        c32 E81 = src1;
        c32 E82 = src2;
        c32 E83 = src3;
        c32 O80 = src4;
        c32 O81 = complex32(wReal.e[1], wImag.e[1]) * src5;
        c32 O82 = complex32(wReal.e[2], wImag.e[2]) * src6;
        c32 O83 = complex32(wReal.e[3], wImag.e[3]) * src7;

        src0 = E80 + O80;
        src1 = E81 + O81;
        src2 = E82 + O82;
        src3 = E83 + O83;
        src4 = E80 - O80;
        src5 = E81 - O81;
        src6 = E82 - O82;
        src7 = E83 - O83;

        *(dest + k + 0) = src0;
        *(dest + k + 1) = src1;
        *(dest + k + 2) = src2;
        *(dest + k + 3) = src3;
        *(dest + k + 4) = src4;
        *(dest + k + 5) = src5;
        *(dest + k + 6) = src6;
        *(dest + k + 7) = src7;
    }

    u32 halfM = 8;
    u32 m = 16;
    while (m <= count)
    {
        f32 oneOverM = (2.0f * F32_PI) / (f32)m;

        for (u32 k = 0; k < count; k += m)
        {
            c32 *src0 = dest + k;
            c32 *src1 = dest + k + halfM;
            for (u32 j = 0; j < halfM; j += 8)
            {
                SinCos32_4x cs0 = sincos32_4x(F32_4x(-(f32)(j + 0) * oneOverM,
                                                     -(f32)(j + 1) * oneOverM,
                                                     -(f32)(j + 2) * oneOverM,
                                                     -(f32)(j + 3) * oneOverM));
                SinCos32_4x cs1 = sincos32_4x(F32_4x(-(f32)(j + 4) * oneOverM,
                                                     -(f32)(j + 5) * oneOverM,
                                                     -(f32)(j + 6) * oneOverM,
                                                     -(f32)(j + 7) * oneOverM));
                f32_4x wReal0 = cs0.cos;
                f32_4x wImag0 = cs0.sin;
                f32_4x wReal1 = cs1.cos;
                f32_4x wImag1 = cs1.sin;
                c32 w0 = complex32(wReal0.e[0], wImag0.e[0]);
                c32 w1 = complex32(wReal0.e[1], wImag0.e[1]);
                c32 w2 = complex32(wReal0.e[2], wImag0.e[2]);
                c32 w3 = complex32(wReal0.e[3], wImag0.e[3]);
                c32 w4 = complex32(wReal1.e[0], wImag1.e[0]);
                c32 w5 = complex32(wReal1.e[1], wImag1.e[1]);
                c32 w6 = complex32(wReal1.e[2], wImag1.e[2]);
                c32 w7 = complex32(wReal1.e[3], wImag1.e[3]);

                c32 src00 = *(src0 + 0);
                c32 src01 = *(src0 + 1);
                c32 src02 = *(src0 + 2);
                c32 src03 = *(src0 + 3);
                c32 src04 = *(src0 + 4);
                c32 src05 = *(src0 + 5);
                c32 src06 = *(src0 + 6);
                c32 src07 = *(src0 + 7);

                c32 src10 = *(src1 + 0);
                c32 src11 = *(src1 + 1);
                c32 src12 = *(src1 + 2);
                c32 src13 = *(src1 + 3);
                c32 src14 = *(src1 + 4);
                c32 src15 = *(src1 + 5);
                c32 src16 = *(src1 + 6);
                c32 src17 = *(src1 + 7);

                c32 E0 = src00;
                c32 E1 = src01;
                c32 E2 = src02;
                c32 E3 = src03;
                c32 E4 = src04;
                c32 E5 = src05;
                c32 E6 = src06;
                c32 E7 = src07;
                c32 O0 = w0 * src10;
                c32 O1 = w1 * src11;
                c32 O2 = w2 * src12;
                c32 O3 = w3 * src13;
                c32 O4 = w4 * src14;
                c32 O5 = w5 * src15;
                c32 O6 = w6 * src16;
                c32 O7 = w7 * src17;

                src00 = E0 + O0;
                src01 = E1 + O1;
                src02 = E2 + O2;
                src03 = E3 + O3;
                src04 = E4 + O4;
                src05 = E5 + O5;
                src06 = E6 + O6;
                src07 = E7 + O7;
                src10 = E0 - O0;
                src11 = E1 - O1;
                src12 = E2 - O2;
                src13 = E3 - O3;
                src14 = E4 - O4;
                src15 = E5 - O5;
                src16 = E6 - O6;
                src17 = E7 - O7;

                *(src0 + 0) = src00;
                *(src0 + 1) = src01;
                *(src0 + 2) = src02;
                *(src0 + 3) = src03;
                *(src0 + 4) = src04;
                *(src0 + 5) = src05;
                *(src0 + 6) = src06;
                *(src0 + 7) = src07;

                *(src1 + 0) = src10;
                *(src1 + 1) = src11;
                *(src1 + 2) = src12;
                *(src1 + 3) = src13;
                *(src1 + 4) = src14;
                *(src1 + 5) = src15;
                *(src1 + 6) = src16;
                *(src1 + 7) = src17;

                src0 += 8;
                src1 += 8;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_normal4(u32 count, c32 *signal, c32 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < halfCount; index += 2)
    {
        u32 index0 = index + 0;
        u32 index1 = index + 1;
        u32 index2 = index0 + halfCount;
        u32 index3 = index1 + halfCount;

        u32 revIndex0 = reverse_bits32(index0, highBit.index);
        u32 revIndex1 = revIndex0 ^ halfCount;
        u32 revIndex2 = revIndex0 + 1;
        u32 revIndex3 = revIndex1 + 1;
        dest[index0] = signal[revIndex0];
        dest[index1] = signal[revIndex1];
        dest[index2] = signal[revIndex2];
        dest[index3] = signal[revIndex3];
    }

    // NOTE(michiel): w = e^(-i*2pi*j/m)
    SinCos32_4x csBase = sincos32_4x(F32_4x(0.0f, -0.25f * F32_PI, -0.5f * F32_PI, -0.75f * F32_PI));
    f32_4x wReal = csBase.cos;
    f32_4x wImag = csBase.sin;

    for (u32 k = 0; k < count; k += 8)
    {
        c32 src0 = *(dest + k + 0);
        c32 src1 = *(dest + k + 1);
        c32 src2 = *(dest + k + 2);
        c32 src3 = *(dest + k + 3);
        c32 src4 = *(dest + k + 4);
        c32 src5 = *(dest + k + 5);
        c32 src6 = *(dest + k + 6);
        c32 src7 = *(dest + k + 7);

        // NOTE(michiel): 2 -> 2
        f32 E20r = src0.real;
        f32 E20i = src0.imag;
        f32 O20r = src1.real;
        f32 O20i = src1.imag;

        f32 E21r = src2.real;
        f32 E21i = src2.imag;
        f32 O21r = src3.real;
        f32 O21i = src3.imag;

        f32 E22r = src4.real;
        f32 E22i = src4.imag;
        f32 O22r = src5.real;
        f32 O22i = src5.imag;

        f32 E23r = src6.real;
        f32 E23i = src6.imag;
        f32 O23r = src7.real;
        f32 O23i = src7.imag;

        f32 o41hr = E21r - O21r;
        f32 o41hi = E21i - O21i;
        f32 o43hr = E23r - O23r;
        f32 o43hi = E23i - O23i;

        // NOTE(michiel): 4 -> 4
        f32 E40r = E20r + O20r;
        f32 E40i = E20i + O20i;
        f32 E41r = E20r - O20r;
        f32 E41i = E20i - O20i;

        f32 O40r = E21r + O21r;
        f32 O40i = E21i + O21i;
        f32 O41r = wReal.e[2] * o41hr - wImag.e[2] * o41hi;
        f32 O41i = wReal.e[2] * o41hi + wImag.e[2] * o41hr;

        f32 E42r = E22r + O22r;
        f32 E42i = E22i + O22i;
        f32 E43r = E22r - O22r;
        f32 E43i = E22i - O22i;

        f32 O42r = E23r + O23r;
        f32 O42i = E23i + O23i;
        f32 O43r = wReal.e[2] * o43hr - wImag.e[2] * o43hi;
        f32 O43i = wReal.e[2] * o43hi + wImag.e[2] * o43hr;

        f32 o80hr = E42r + O42r;
        f32 o80hi = E42i + O42i;
        f32 o81hr = E43r + O43r;
        f32 o81hi = E43i + O43i;
        f32 o82hr = E42r - O42r;
        f32 o82hi = E42i - O42i;
        f32 o83hr = E43r - O43r;
        f32 o83hi = E43i - O43i;

        // NOTE(michiel): 8 -> 8
        f32 E80r = E40r + O40r;
        f32 E80i = E40i + O40i;
        f32 E81r = E41r + O41r;
        f32 E81i = E41i + O41i;
        f32 E82r = E40r - O40r;
        f32 E82i = E40i - O40i;
        f32 E83r = E41r - O41r;
        f32 E83i = E41i - O41i;
        f32 O80r = wReal.e[0] * o80hr - wImag.e[0] * o80hi;
        f32 O80i = wReal.e[0] * o80hi + wImag.e[0] * o80hr;
        f32 O81r = wReal.e[1] * o81hr - wImag.e[1] * o81hi;
        f32 O81i = wReal.e[1] * o81hi + wImag.e[1] * o81hr;

        f32 O82r = wReal.e[2] * o82hr - wImag.e[2] * o82hi;
        f32 O82i = wReal.e[2] * o82hi + wImag.e[2] * o82hr;
        f32 O83r = wReal.e[3] * o83hr - wImag.e[3] * o83hi;
        f32 O83i = wReal.e[3] * o83hi + wImag.e[3] * o83hr;

        src0.real = E80r + O80r;
        src0.imag = E80i + O80i;
        src1.real = E81r + O81r;
        src1.imag = E81i + O81i;

        src2.real = E82r + O82r;
        src2.imag = E82i + O82i;
        src3.real = E83r + O83r;
        src3.imag = E83i + O83i;

        src4.real = E80r - O80r;
        src4.imag = E80i - O80i;
        src5.real = E81r - O81r;
        src5.imag = E81i - O81i;

        src6.real = E82r - O82r;
        src6.imag = E82i - O82i;
        src7.real = E83r - O83r;
        src7.imag = E83i - O83i;

        *(dest + k + 0) = src0;
        *(dest + k + 1) = src1;
        *(dest + k + 2) = src2;
        *(dest + k + 3) = src3;
        *(dest + k + 4) = src4;
        *(dest + k + 5) = src5;
        *(dest + k + 6) = src6;
        *(dest + k + 7) = src7;
    }

    u32 halfM = 8;
    u32 m = 16;
    while (m <= count)
    {
        f32 oneOverM = (2.0f * F32_PI) / (f32)m;

        for (u32 k = 0; k < count; k += m)
        {
            c32 *src0 = dest + k;
            c32 *src1 = dest + k + halfM;
            for (u32 j = 0; j < halfM; j += 8)
            {
                SinCos32_4x cs0 = sincos32_4x(F32_4x(-(f32)(j + 0) * oneOverM,
                                                     -(f32)(j + 1) * oneOverM,
                                                     -(f32)(j + 2) * oneOverM,
                                                     -(f32)(j + 3) * oneOverM));
                SinCos32_4x cs1 = sincos32_4x(F32_4x(-(f32)(j + 4) * oneOverM,
                                                     -(f32)(j + 5) * oneOverM,
                                                     -(f32)(j + 6) * oneOverM,
                                                     -(f32)(j + 7) * oneOverM));

                c32 src00 = *(src0 + 0);
                c32 src01 = *(src0 + 1);
                c32 src02 = *(src0 + 2);
                c32 src03 = *(src0 + 3);
                c32 src04 = *(src0 + 4);
                c32 src05 = *(src0 + 5);
                c32 src06 = *(src0 + 6);
                c32 src07 = *(src0 + 7);

                c32 src10 = *(src1 + 0);
                c32 src11 = *(src1 + 1);
                c32 src12 = *(src1 + 2);
                c32 src13 = *(src1 + 3);
                c32 src14 = *(src1 + 4);
                c32 src15 = *(src1 + 5);
                c32 src16 = *(src1 + 6);
                c32 src17 = *(src1 + 7);

                f32 E0r = src00.real;
                f32 E0i = src00.imag;
                f32 E1r = src01.real;
                f32 E1i = src01.imag;

                f32 E2r = src02.real;
                f32 E2i = src02.imag;
                f32 E3r = src03.real;
                f32 E3i = src03.imag;

                f32 E4r = src04.real;
                f32 E4i = src04.imag;
                f32 E5r = src05.real;
                f32 E5i = src05.imag;

                f32 E6r = src06.real;
                f32 E6i = src06.imag;
                f32 E7r = src07.real;
                f32 E7i = src07.imag;

                f32 O0r = cs0.cos.e[0] * src10.real - cs0.sin.e[0] * src10.imag;
                f32 O0i = cs0.cos.e[0] * src10.imag + cs0.sin.e[0] * src10.real;
                f32 O1r = cs0.cos.e[1] * src11.real - cs0.sin.e[1] * src11.imag;
                f32 O1i = cs0.cos.e[1] * src11.imag + cs0.sin.e[1] * src11.real;

                f32 O2r = cs0.cos.e[2] * src12.real - cs0.sin.e[2] * src12.imag;
                f32 O2i = cs0.cos.e[2] * src12.imag + cs0.sin.e[2] * src12.real;
                f32 O3r = cs0.cos.e[3] * src13.real - cs0.sin.e[3] * src13.imag;
                f32 O3i = cs0.cos.e[3] * src13.imag + cs0.sin.e[3] * src13.real;

                f32 O4r = cs1.cos.e[0] * src14.real - cs1.sin.e[0] * src14.imag;
                f32 O4i = cs1.cos.e[0] * src14.imag + cs1.sin.e[0] * src14.real;
                f32 O5r = cs1.cos.e[1] * src15.real - cs1.sin.e[1] * src15.imag;
                f32 O5i = cs1.cos.e[1] * src15.imag + cs1.sin.e[1] * src15.real;

                f32 O6r = cs1.cos.e[2] * src16.real - cs1.sin.e[2] * src16.imag;
                f32 O6i = cs1.cos.e[2] * src16.imag + cs1.sin.e[2] * src16.real;
                f32 O7r = cs1.cos.e[3] * src17.real - cs1.sin.e[3] * src17.imag;
                f32 O7i = cs1.cos.e[3] * src17.imag + cs1.sin.e[3] * src17.real;

                src00.real = E0r + O0r;
                src00.imag = E0i + O0i;
                src01.real = E1r + O1r;
                src01.imag = E1i + O1i;

                src02.real = E2r + O2r;
                src02.imag = E2i + O2i;
                src03.real = E3r + O3r;
                src03.imag = E3i + O3i;

                src04.real = E4r + O4r;
                src04.imag = E4i + O4i;
                src05.real = E5r + O5r;
                src05.imag = E5i + O5i;

                src06.real = E6r + O6r;
                src06.imag = E6i + O6i;
                src07.real = E7r + O7r;
                src07.imag = E7i + O7i;

                src10.real = E0r - O0r;
                src10.imag = E0i - O0i;
                src11.real = E1r - O1r;
                src11.imag = E1i - O1i;

                src12.real = E2r - O2r;
                src12.imag = E2i - O2i;
                src13.real = E3r - O3r;
                src13.imag = E3i - O3i;

                src14.real = E4r - O4r;
                src14.imag = E4i - O4i;
                src15.real = E5r - O5r;
                src15.imag = E5i - O5i;

                src16.real = E6r - O6r;
                src16.imag = E6i - O6i;
                src17.real = E7r - O7r;
                src17.imag = E7i - O7i;

                *(src0 + 0) = src00;
                *(src0 + 1) = src01;
                *(src0 + 2) = src02;
                *(src0 + 3) = src03;
                *(src0 + 4) = src04;
                *(src0 + 5) = src05;
                *(src0 + 6) = src06;
                *(src0 + 7) = src07;

                *(src1 + 0) = src10;
                *(src1 + 1) = src11;
                *(src1 + 2) = src12;
                *(src1 + 3) = src13;
                *(src1 + 4) = src14;
                *(src1 + 5) = src15;
                *(src1 + 6) = src16;
                *(src1 + 7) = src17;

                src0 += 8;
                src1 += 8;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_normal5(u32 count, c32 *signal, c32 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < halfCount; index += 2)
    {
        u32 index0 = index + 0;
        u32 index1 = index + 1;
        u32 index2 = index0 + halfCount;
        u32 index3 = index1 + halfCount;

        u32 revIndex0 = reverse_bits32(index0, highBit.index);
        u32 revIndex1 = revIndex0 ^ halfCount;
        u32 revIndex2 = revIndex0 + 1;
        u32 revIndex3 = revIndex1 + 1;
        dest[index0] = signal[revIndex0];
        dest[index1] = signal[revIndex1];
        dest[index2] = signal[revIndex2];
        dest[index3] = signal[revIndex3];
    }

    // NOTE(michiel): w = e^(-i*2pi*j/m)
    SinCos32_4x csBase = sincos32_4x(F32_4x(0.0f, -0.25f * F32_PI, -0.5f * F32_PI, -0.75f * F32_PI));
    f32_4x cos42; cos42.m = _mm_shuffle_ps(csBase.cos.m, csBase.cos.m, MULTILANE_SHUFFLE_MASK(2, 2, 2, 2));
    f32_4x sin42; sin42.m = _mm_shuffle_ps(csBase.sin.m, csBase.sin.m, MULTILANE_SHUFFLE_MASK(2, 2, 2, 2));
    f32_4x cos801; cos801.m = _mm_shuffle_ps(csBase.cos.m, csBase.cos.m, MULTILANE_SHUFFLE_MASK(0, 0, 1, 1));
    f32_4x sin801; sin801.m = _mm_shuffle_ps(csBase.sin.m, csBase.sin.m, MULTILANE_SHUFFLE_MASK(0, 0, 1, 1));
    f32_4x cos823; cos823.m = _mm_shuffle_ps(csBase.cos.m, csBase.cos.m, MULTILANE_SHUFFLE_MASK(2, 2, 3, 3));
    f32_4x sin823; sin823.m = _mm_shuffle_ps(csBase.sin.m, csBase.sin.m, MULTILANE_SHUFFLE_MASK(2, 2, 3, 3));

    for (u32 k = 0; k < count; k += 8)
    {
        f32_4x EO20ri = F32_4x((f32 *)(dest + k + 0));
        f32_4x EO21ri = F32_4x((f32 *)(dest + k + 2));
        f32_4x EO22ri = F32_4x((f32 *)(dest + k + 4));
        f32_4x EO23ri = F32_4x((f32 *)(dest + k + 6));

        f32_4x E202ri; E202ri.m = _mm_shuffle_ps(EO20ri.m, EO22ri.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
        f32_4x O202ri; O202ri.m = _mm_shuffle_ps(EO20ri.m, EO22ri.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3));
        f32_4x E213ri; E213ri.m = _mm_shuffle_ps(EO21ri.m, EO23ri.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
        f32_4x O213ri; O213ri.m = _mm_shuffle_ps(EO21ri.m, EO23ri.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3));

        f32_4x E402ri = E202ri + O202ri;
        f32_4x E413ri = E202ri - O202ri;

        f32_4x O402ri = E213ri + O213ri;
        f32_4x o413hri = E213ri - O213ri;
        f32_4x o413hir; o413hir.m = _mm_shuffle_ps(o413hri.m, o413hri.m, MULTILANE_SHUFFLE_MASK(1, 0, 3, 2));
        f32_4x O413ri; O413ri.m = _mm_addsub_ps((cos42 * o413hri).m, (sin42 * o413hir).m);

        f32_4x EpO402ri = E402ri + O402ri;
        f32_4x EmO402ri = E402ri - O402ri;
        f32_4x EpO413ri = E413ri + O413ri;
        f32_4x EmO413ri = E413ri - O413ri;

        f32_4x E801ri; E801ri.m = _mm_shuffle_ps(EpO402ri.m, EpO413ri.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
        f32_4x E823ri; E823ri.m = _mm_shuffle_ps(EmO402ri.m, EmO413ri.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));

        f32_4x o801hri; o801hri.m = _mm_shuffle_ps(EpO402ri.m, EpO413ri.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3));
        f32_4x o801hir; o801hir.m = _mm_shuffle_ps(EpO402ri.m, EpO413ri.m, MULTILANE_SHUFFLE_MASK(3, 2, 3, 2));
        f32_4x O801ri; O801ri.m = _mm_addsub_ps((cos801 * o801hri).m, (sin801 * o801hir).m);

        f32_4x o823hri; o823hri.m = _mm_shuffle_ps(EmO402ri.m, EmO413ri.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3));
        f32_4x o823hir; o823hir.m = _mm_shuffle_ps(EmO402ri.m, EmO413ri.m, MULTILANE_SHUFFLE_MASK(3, 2, 3, 2));
        f32_4x O823ri; O823ri.m = _mm_addsub_ps((cos823 * o823hri).m, (sin823 * o823hir).m);

        f32_4x EpO801 = E801ri + O801ri;
        f32_4x EpO823 = E823ri + O823ri;
        f32_4x EmO801 = E801ri - O801ri;
        f32_4x EmO823 = E823ri - O823ri;

        _mm_store_ps((f32 *)(dest + k + 0), EpO801.m);
        _mm_store_ps((f32 *)(dest + k + 2), EpO823.m);
        _mm_store_ps((f32 *)(dest + k + 4), EmO801.m);
        _mm_store_ps((f32 *)(dest + k + 6), EmO823.m);
    }

    u32 halfM = 8;
    u32 m = 16;
    while (m <= count)
    {
        f32 oneOverM = (2.0f * F32_PI) / (f32)m;

        for (u32 k = 0; k < count; k += m)
        {
            c32 *src0 = dest + k;
            c32 *src1 = dest + k + halfM;
            for (u32 j = 0; j < halfM; j += 8)
            {
                SinCos32_4x cs0 = sincos32_4x(F32_4x(-(f32)(j + 0) * oneOverM,
                                                     -(f32)(j + 1) * oneOverM,
                                                     -(f32)(j + 2) * oneOverM,
                                                     -(f32)(j + 3) * oneOverM));
                SinCos32_4x cs1 = sincos32_4x(F32_4x(-(f32)(j + 4) * oneOverM,
                                                     -(f32)(j + 5) * oneOverM,
                                                     -(f32)(j + 6) * oneOverM,
                                                     -(f32)(j + 7) * oneOverM));

                c32 src00 = *(src0 + 0);
                c32 src01 = *(src0 + 1);
                c32 src02 = *(src0 + 2);
                c32 src03 = *(src0 + 3);
                c32 src04 = *(src0 + 4);
                c32 src05 = *(src0 + 5);
                c32 src06 = *(src0 + 6);
                c32 src07 = *(src0 + 7);

                c32 src10 = *(src1 + 0);
                c32 src11 = *(src1 + 1);
                c32 src12 = *(src1 + 2);
                c32 src13 = *(src1 + 3);
                c32 src14 = *(src1 + 4);
                c32 src15 = *(src1 + 5);
                c32 src16 = *(src1 + 6);
                c32 src17 = *(src1 + 7);

                f32 E0r = src00.real;
                f32 E0i = src00.imag;
                f32 E1r = src01.real;
                f32 E1i = src01.imag;

                f32 E2r = src02.real;
                f32 E2i = src02.imag;
                f32 E3r = src03.real;
                f32 E3i = src03.imag;

                f32 E4r = src04.real;
                f32 E4i = src04.imag;
                f32 E5r = src05.real;
                f32 E5i = src05.imag;

                f32 E6r = src06.real;
                f32 E6i = src06.imag;
                f32 E7r = src07.real;
                f32 E7i = src07.imag;

                f32 O0r = cs0.cos.e[0] * src10.real - cs0.sin.e[0] * src10.imag;
                f32 O0i = cs0.cos.e[0] * src10.imag + cs0.sin.e[0] * src10.real;
                f32 O1r = cs0.cos.e[1] * src11.real - cs0.sin.e[1] * src11.imag;
                f32 O1i = cs0.cos.e[1] * src11.imag + cs0.sin.e[1] * src11.real;

                f32 O2r = cs0.cos.e[2] * src12.real - cs0.sin.e[2] * src12.imag;
                f32 O2i = cs0.cos.e[2] * src12.imag + cs0.sin.e[2] * src12.real;
                f32 O3r = cs0.cos.e[3] * src13.real - cs0.sin.e[3] * src13.imag;
                f32 O3i = cs0.cos.e[3] * src13.imag + cs0.sin.e[3] * src13.real;

                f32 O4r = cs1.cos.e[0] * src14.real - cs1.sin.e[0] * src14.imag;
                f32 O4i = cs1.cos.e[0] * src14.imag + cs1.sin.e[0] * src14.real;
                f32 O5r = cs1.cos.e[1] * src15.real - cs1.sin.e[1] * src15.imag;
                f32 O5i = cs1.cos.e[1] * src15.imag + cs1.sin.e[1] * src15.real;

                f32 O6r = cs1.cos.e[2] * src16.real - cs1.sin.e[2] * src16.imag;
                f32 O6i = cs1.cos.e[2] * src16.imag + cs1.sin.e[2] * src16.real;
                f32 O7r = cs1.cos.e[3] * src17.real - cs1.sin.e[3] * src17.imag;
                f32 O7i = cs1.cos.e[3] * src17.imag + cs1.sin.e[3] * src17.real;

                src00.real = E0r + O0r;
                src00.imag = E0i + O0i;
                src01.real = E1r + O1r;
                src01.imag = E1i + O1i;

                src02.real = E2r + O2r;
                src02.imag = E2i + O2i;
                src03.real = E3r + O3r;
                src03.imag = E3i + O3i;

                src04.real = E4r + O4r;
                src04.imag = E4i + O4i;
                src05.real = E5r + O5r;
                src05.imag = E5i + O5i;

                src06.real = E6r + O6r;
                src06.imag = E6i + O6i;
                src07.real = E7r + O7r;
                src07.imag = E7i + O7i;

                src10.real = E0r - O0r;
                src10.imag = E0i - O0i;
                src11.real = E1r - O1r;
                src11.imag = E1i - O1i;

                src12.real = E2r - O2r;
                src12.imag = E2i - O2i;
                src13.real = E3r - O3r;
                src13.imag = E3i - O3i;

                src14.real = E4r - O4r;
                src14.imag = E4i - O4i;
                src15.real = E5r - O5r;
                src15.imag = E5i - O5i;

                src16.real = E6r - O6r;
                src16.imag = E6i - O6i;
                src17.real = E7r - O7r;
                src17.imag = E7i - O7i;

                *(src0 + 0) = src00;
                *(src0 + 1) = src01;
                *(src0 + 2) = src02;
                *(src0 + 3) = src03;
                *(src0 + 4) = src04;
                *(src0 + 5) = src05;
                *(src0 + 6) = src06;
                *(src0 + 7) = src07;

                *(src1 + 0) = src10;
                *(src1 + 1) = src11;
                *(src1 + 2) = src12;
                *(src1 + 3) = src13;
                *(src1 + 4) = src14;
                *(src1 + 5) = src15;
                *(src1 + 6) = src16;
                *(src1 + 7) = src17;

                src0 += 8;
                src1 += 8;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_inplace6(u32 count, c32 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < halfCount; index += 2)
    {
        u32 index0 = index + 0;
        u32 index1 = index + 1;
        u32 index2 = index0 + halfCount;
        u32 index3 = index1 + halfCount;

        u32 revIndex0 = reverse_bits32(index0, highBit.index);
        u32 revIndex1 = revIndex0 ^ halfCount;
        u32 revIndex2 = revIndex0 + 1;
        u32 revIndex3 = revIndex1 + 1;
        if (revIndex0 > index0) {
            c32 temp = dest[index0];
            dest[index0] = dest[revIndex0];
            dest[revIndex0] = temp;
        }
        if (revIndex1 > index1) {
            c32 temp = dest[index1];
            dest[index1] = dest[revIndex1];
            dest[revIndex1] = temp;
        }
        if (revIndex2 > index2) {
            c32 temp = dest[index2];
            dest[index2] = dest[revIndex2];
            dest[revIndex2] = temp;
        }
        if (revIndex3 > index3) {
            c32 temp = dest[index3];
            dest[index3] = dest[revIndex3];
            dest[revIndex3] = temp;
        }
    }

    // NOTE(michiel): w = e^(-i*2pi*j/m)
    SinCos32_4x csBase = sincos32_4x(F32_4x(0.0f, -0.25f * F32_PI, -0.5f * F32_PI, -0.75f * F32_PI));
    f32_4x cos42; cos42.m = _mm_shuffle_ps(csBase.cos.m, csBase.cos.m, MULTILANE_SHUFFLE_MASK(2, 2, 2, 2));
    f32_4x sin42; sin42.m = _mm_shuffle_ps(csBase.sin.m, csBase.sin.m, MULTILANE_SHUFFLE_MASK(2, 2, 2, 2));
    f32_4x cos801; cos801.m = _mm_shuffle_ps(csBase.cos.m, csBase.cos.m, MULTILANE_SHUFFLE_MASK(0, 0, 1, 1));
    f32_4x sin801; sin801.m = _mm_shuffle_ps(csBase.sin.m, csBase.sin.m, MULTILANE_SHUFFLE_MASK(0, 0, 1, 1));
    f32_4x cos823; cos823.m = _mm_shuffle_ps(csBase.cos.m, csBase.cos.m, MULTILANE_SHUFFLE_MASK(2, 2, 3, 3));
    f32_4x sin823; sin823.m = _mm_shuffle_ps(csBase.sin.m, csBase.sin.m, MULTILANE_SHUFFLE_MASK(2, 2, 3, 3));

    for (u32 k = 0; k < count; k += 8)
    {
        f32_4x EO20ri = F32_4x((f32 *)(dest + k + 0));
        f32_4x EO21ri = F32_4x((f32 *)(dest + k + 2));
        f32_4x EO22ri = F32_4x((f32 *)(dest + k + 4));
        f32_4x EO23ri = F32_4x((f32 *)(dest + k + 6));

        f32_4x E202ri; E202ri.m = _mm_shuffle_ps(EO20ri.m, EO22ri.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
        f32_4x O202ri; O202ri.m = _mm_shuffle_ps(EO20ri.m, EO22ri.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3));
        f32_4x E213ri; E213ri.m = _mm_shuffle_ps(EO21ri.m, EO23ri.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
        f32_4x O213ri; O213ri.m = _mm_shuffle_ps(EO21ri.m, EO23ri.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3));

        f32_4x E402ri = E202ri + O202ri;
        f32_4x E413ri = E202ri - O202ri;

        f32_4x O402ri = E213ri + O213ri;
        f32_4x o413hri = E213ri - O213ri;
        f32_4x o413hir; o413hir.m = _mm_shuffle_ps(o413hri.m, o413hri.m, MULTILANE_SHUFFLE_MASK(1, 0, 3, 2));
        f32_4x O413ri; O413ri.m = _mm_addsub_ps((cos42 * o413hri).m, (sin42 * o413hir).m);

        f32_4x EpO402ri = E402ri + O402ri;
        f32_4x EmO402ri = E402ri - O402ri;
        f32_4x EpO413ri = E413ri + O413ri;
        f32_4x EmO413ri = E413ri - O413ri;

        f32_4x E801ri; E801ri.m = _mm_shuffle_ps(EpO402ri.m, EpO413ri.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
        f32_4x E823ri; E823ri.m = _mm_shuffle_ps(EmO402ri.m, EmO413ri.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));

        f32_4x o801hri; o801hri.m = _mm_shuffle_ps(EpO402ri.m, EpO413ri.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3));
        f32_4x o801hir; o801hir.m = _mm_shuffle_ps(EpO402ri.m, EpO413ri.m, MULTILANE_SHUFFLE_MASK(3, 2, 3, 2));
        f32_4x O801ri; O801ri.m = _mm_addsub_ps((cos801 * o801hri).m, (sin801 * o801hir).m);

        f32_4x o823hri; o823hri.m = _mm_shuffle_ps(EmO402ri.m, EmO413ri.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3));
        f32_4x o823hir; o823hir.m = _mm_shuffle_ps(EmO402ri.m, EmO413ri.m, MULTILANE_SHUFFLE_MASK(3, 2, 3, 2));
        f32_4x O823ri; O823ri.m = _mm_addsub_ps((cos823 * o823hri).m, (sin823 * o823hir).m);

        f32_4x EpO801 = E801ri + O801ri;
        f32_4x EpO823 = E823ri + O823ri;
        f32_4x EmO801 = E801ri - O801ri;
        f32_4x EmO823 = E823ri - O823ri;

        _mm_store_ps((f32 *)(dest + k + 0), EpO801.m);
        _mm_store_ps((f32 *)(dest + k + 2), EpO823.m);
        _mm_store_ps((f32 *)(dest + k + 4), EmO801.m);
        _mm_store_ps((f32 *)(dest + k + 6), EmO823.m);
    }

    u32 halfM = 8;
    u32 m = 16;
    f32_4x angleMod0 = F32_4x(0.0f, 1.0f, 2.0f, 3.0f);

    while (m <= count)
    {
        f32_4x oneOverM = F32_4x((2.0f * F32_PI) / (f32)m);
        f32_4x angleBase = -oneOverM * angleMod0;
        f32_4x angleStep = F32_4x(-4.0f) * oneOverM;

        for (u32 k = 0; k < count; k += m)
        {
            c32 *src0 = dest + k;
            c32 *src1 = dest + k + halfM;

            f32_4x angles = angleBase;

            f32 *EGrab = (f32 *)src0;
            f32 *OGrab = (f32 *)src1;
            f32 *EPut  = (f32 *)src0;
            f32 *OPut  = (f32 *)src1;

            for (u32 j = 0; j < halfM; j += 8)
            {
                f32_4x src101ri = F32_4x(OGrab);
                OGrab += 4;
                f32_4x src123ri = F32_4x(OGrab);
                OGrab += 4;

                f32_4x src103r; src103r.m = _mm_shuffle_ps(src101ri.m, src123ri.m, MULTILANE_SHUFFLE_MASK(0, 2, 0, 2));
                f32_4x src103i; src103i.m = _mm_shuffle_ps(src101ri.m, src123ri.m, MULTILANE_SHUFFLE_MASK(1, 3, 1, 3));

                SinCos32_4x cs0 = sincos32_4x(angles);
                angles = angles + angleStep;

                f32_4x O03r = cs0.cos * src103r - cs0.sin * src103i;
                f32_4x O03i = cs0.cos * src103i + cs0.sin * src103r;

                f32_4x E01ri = F32_4x(EGrab);
                EGrab += 4;
                f32_4x E23ri = F32_4x(EGrab);
                EGrab += 4;

                f32_4x O01ri; O01ri.m = _mm_unpacklo_ps(O03r.m, O03i.m);
                f32_4x O23ri; O23ri.m = _mm_unpackhi_ps(O03r.m, O03i.m);

                f32_4x EpO01ri = E01ri + O01ri;
                f32_4x EmO01ri = E01ri - O01ri;
                f32_4x EpO23ri = E23ri + O23ri;
                f32_4x EmO23ri = E23ri - O23ri;

                _mm_store_ps(EPut, EpO01ri.m);
                EPut += 4;
                _mm_store_ps(OPut, EmO01ri.m);
                OPut += 4;
                _mm_store_ps(EPut, EpO23ri.m);
                EPut += 4;
                _mm_store_ps(OPut, EmO23ri.m);
                OPut += 4;

                f32_4x src145ri = F32_4x(OGrab);
                OGrab += 4;
                f32_4x src167ri = F32_4x(OGrab);
                OGrab += 4;

                f32_4x src147r; src147r.m = _mm_shuffle_ps(src145ri.m, src167ri.m, MULTILANE_SHUFFLE_MASK(0, 2, 0, 2));
                f32_4x src147i; src147i.m = _mm_shuffle_ps(src145ri.m, src167ri.m, MULTILANE_SHUFFLE_MASK(1, 3, 1, 3));

                SinCos32_4x cs1 = sincos32_4x(angles);
                angles = angles + angleStep;

                f32_4x O47r = cs1.cos * src147r - cs1.sin * src147i;
                f32_4x O47i = cs1.cos * src147i + cs1.sin * src147r;

                f32_4x E45ri = F32_4x(EGrab);
                EGrab += 4;
                f32_4x E67ri = F32_4x(EGrab);
                EGrab += 4;

                f32_4x O45ri; O45ri.m = _mm_unpacklo_ps(O47r.m, O47i.m);
                f32_4x O67ri; O67ri.m = _mm_unpackhi_ps(O47r.m, O47i.m);

                f32_4x EpO45ri = E45ri + O45ri;
                f32_4x EmO45ri = E45ri - O45ri;
                f32_4x EpO67ri = E67ri + O67ri;
                f32_4x EmO67ri = E67ri - O67ri;

                _mm_store_ps(EPut, EpO45ri.m);
                EPut += 4;
                _mm_store_ps(OPut, EmO45ri.m);
                OPut += 4;
                _mm_store_ps(EPut, EpO67ri.m);
                EPut += 4;
                _mm_store_ps(OPut, EmO67ri.m);
                OPut += 4;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_normal6(u32 count, c32 *signal, c32 *dest)
{
    copy(count * sizeof(c32), signal, dest);
    fft_inplace6(count, dest);
}

internal void
fft_inplace_inexact6(u32 count, c32 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < halfCount; index += 2)
    {
        u32 index0 = index + 0;
        u32 index1 = index + 1;
        u32 index2 = index0 + halfCount;
        u32 index3 = index1 + halfCount;

        u32 revIndex0 = reverse_bits32(index0, highBit.index);
        u32 revIndex1 = revIndex0 ^ halfCount;
        u32 revIndex2 = revIndex0 + 1;
        u32 revIndex3 = revIndex1 + 1;
        if (revIndex0 > index0) {
            c32 temp = dest[index0];
            dest[index0] = dest[revIndex0];
            dest[revIndex0] = temp;
        }
        if (revIndex1 > index1) {
            c32 temp = dest[index1];
            dest[index1] = dest[revIndex1];
            dest[revIndex1] = temp;
        }
        if (revIndex2 > index2) {
            c32 temp = dest[index2];
            dest[index2] = dest[revIndex2];
            dest[revIndex2] = temp;
        }
        if (revIndex3 > index3) {
            c32 temp = dest[index3];
            dest[index3] = dest[revIndex3];
            dest[revIndex3] = temp;
        }
    }

    // NOTE(michiel): w = e^(-i*2pi*j/m)
    SinCos32_4x csBase = sincos32_4x(F32_4x(0.0f, -0.25f * F32_PI, -0.5f * F32_PI, -0.75f * F32_PI));
    f32_4x cos42; cos42.m = _mm_shuffle_ps(csBase.cos.m, csBase.cos.m, MULTILANE_SHUFFLE_MASK(2, 2, 2, 2));
    f32_4x sin42; sin42.m = _mm_shuffle_ps(csBase.sin.m, csBase.sin.m, MULTILANE_SHUFFLE_MASK(2, 2, 2, 2));
    f32_4x cos801; cos801.m = _mm_shuffle_ps(csBase.cos.m, csBase.cos.m, MULTILANE_SHUFFLE_MASK(0, 0, 1, 1));
    f32_4x sin801; sin801.m = _mm_shuffle_ps(csBase.sin.m, csBase.sin.m, MULTILANE_SHUFFLE_MASK(0, 0, 1, 1));
    f32_4x cos823; cos823.m = _mm_shuffle_ps(csBase.cos.m, csBase.cos.m, MULTILANE_SHUFFLE_MASK(2, 2, 3, 3));
    f32_4x sin823; sin823.m = _mm_shuffle_ps(csBase.sin.m, csBase.sin.m, MULTILANE_SHUFFLE_MASK(2, 2, 3, 3));

    for (u32 k = 0; k < count; k += 8)
    {
        f32_4x EO20ri = F32_4x((f32 *)(dest + k + 0));
        f32_4x EO21ri = F32_4x((f32 *)(dest + k + 2));
        f32_4x EO22ri = F32_4x((f32 *)(dest + k + 4));
        f32_4x EO23ri = F32_4x((f32 *)(dest + k + 6));

        f32_4x E202ri; E202ri.m = _mm_shuffle_ps(EO20ri.m, EO22ri.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
        f32_4x O202ri; O202ri.m = _mm_shuffle_ps(EO20ri.m, EO22ri.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3));
        f32_4x E213ri; E213ri.m = _mm_shuffle_ps(EO21ri.m, EO23ri.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
        f32_4x O213ri; O213ri.m = _mm_shuffle_ps(EO21ri.m, EO23ri.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3));

        f32_4x E402ri = E202ri + O202ri;
        f32_4x E413ri = E202ri - O202ri;

        f32_4x O402ri = E213ri + O213ri;
        f32_4x o413hri = E213ri - O213ri;
        f32_4x o413hir; o413hir.m = _mm_shuffle_ps(o413hri.m, o413hri.m, MULTILANE_SHUFFLE_MASK(1, 0, 3, 2));
        f32_4x O413ri; O413ri.m = _mm_addsub_ps((cos42 * o413hri).m, (sin42 * o413hir).m);

        f32_4x EpO402ri = E402ri + O402ri;
        f32_4x EmO402ri = E402ri - O402ri;
        f32_4x EpO413ri = E413ri + O413ri;
        f32_4x EmO413ri = E413ri - O413ri;

        f32_4x E801ri; E801ri.m = _mm_shuffle_ps(EpO402ri.m, EpO413ri.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));
        f32_4x E823ri; E823ri.m = _mm_shuffle_ps(EmO402ri.m, EmO413ri.m, MULTILANE_SHUFFLE_MASK(0, 1, 0, 1));

        f32_4x o801hri; o801hri.m = _mm_shuffle_ps(EpO402ri.m, EpO413ri.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3));
        f32_4x o801hir; o801hir.m = _mm_shuffle_ps(EpO402ri.m, EpO413ri.m, MULTILANE_SHUFFLE_MASK(3, 2, 3, 2));
        f32_4x O801ri; O801ri.m = _mm_addsub_ps((cos801 * o801hri).m, (sin801 * o801hir).m);

        f32_4x o823hri; o823hri.m = _mm_shuffle_ps(EmO402ri.m, EmO413ri.m, MULTILANE_SHUFFLE_MASK(2, 3, 2, 3));
        f32_4x o823hir; o823hir.m = _mm_shuffle_ps(EmO402ri.m, EmO413ri.m, MULTILANE_SHUFFLE_MASK(3, 2, 3, 2));
        f32_4x O823ri; O823ri.m = _mm_addsub_ps((cos823 * o823hri).m, (sin823 * o823hir).m);

        f32_4x EpO801 = E801ri + O801ri;
        f32_4x EpO823 = E823ri + O823ri;
        f32_4x EmO801 = E801ri - O801ri;
        f32_4x EmO823 = E823ri - O823ri;

        _mm_store_ps((f32 *)(dest + k + 0), EpO801.m);
        _mm_store_ps((f32 *)(dest + k + 2), EpO823.m);
        _mm_store_ps((f32 *)(dest + k + 4), EmO801.m);
        _mm_store_ps((f32 *)(dest + k + 6), EmO823.m);
    }

    u32 halfM = 8;
    u32 m = 16;

    f32 oneOverMpre = 4.0f * F32_PI / (f32)m;
    SinCos32_4x sinCosPre = sincos32_4x(F32_4x(-oneOverMpre, -oneOverMpre*2.0f, -oneOverMpre*3.0f, -oneOverMpre*4.0f));

    c32 wm1 = complex32(sinCosPre.cos.e[0], sinCosPre.sin.e[0]);
    c32 wm2 = complex32(sinCosPre.cos.e[1], sinCosPre.sin.e[1]);
    c32 wm3 = complex32(sinCosPre.cos.e[2], sinCosPre.sin.e[2]);
    c32 wm4 = complex32(sinCosPre.cos.e[3], sinCosPre.sin.e[3]);

    while (m <= count)
    {
        f32 oneOverM = 0.5f * oneOverMpre; // (2.0f * F32_PI) / (f32)m;
        oneOverMpre = oneOverM;

        c32 wm8 = wm4;
        c32 wm6 = wm3;
        wm4 = wm2;
        wm2 = wm1;

        SinCos32_4x sinCos = sincos32_4x(F32_4x(-oneOverM, -oneOverM*3.0f, -oneOverM*5.0f, -oneOverM*7.0f));
        wm1 = complex32(sinCos.cos.e[0], sinCos.sin.e[0]);
        wm3 = complex32(sinCos.cos.e[1], sinCos.sin.e[1]);
        c32 wm5 = complex32(sinCos.cos.e[2], sinCos.sin.e[2]);
        c32 wm7 = complex32(sinCos.cos.e[3], sinCos.sin.e[3]);

        f32_4x wStepr = F32_4x(wm8.real);
        f32_4x wStepi = F32_4x(wm8.imag);

        f32_4x w0Startr = F32_4x(1.0f, wm1.real, wm2.real, wm3.real);
        f32_4x w0Starti = F32_4x(0.0f, wm1.imag, wm2.imag, wm3.imag);
        f32_4x w1Startr = F32_4x(wm4.real, wm5.real, wm6.real, wm7.real);
        f32_4x w1Starti = F32_4x(wm4.imag, wm5.imag, wm6.imag, wm7.imag);

        for (u32 k = 0; k < count; k += m)
        {
            c32 *src0 = dest + k;
            c32 *src1 = dest + k + halfM;

            f32_4x w0r = w0Startr;
            f32_4x w0i = w0Starti;
            f32_4x w1r = w1Startr;
            f32_4x w1i = w1Starti;

            f32 *EGrab = (f32 *)src0;
            f32 *OGrab = (f32 *)src1;
            f32 *EPut  = (f32 *)src0;
            f32 *OPut  = (f32 *)src1;

            for (u32 j = 0; j < halfM; j += 8)
            {
                f32_4x src101ri = F32_4x(OGrab);
                OGrab += 4;
                f32_4x src123ri = F32_4x(OGrab);
                OGrab += 4;

                f32_4x src103r; src103r.m = _mm_shuffle_ps(src101ri.m, src123ri.m, MULTILANE_SHUFFLE_MASK(0, 2, 0, 2));
                f32_4x src103i; src103i.m = _mm_shuffle_ps(src101ri.m, src123ri.m, MULTILANE_SHUFFLE_MASK(1, 3, 1, 3));

                f32_4x O03r = w0r * src103r - w0i * src103i;
                f32_4x O03i = w0r * src103i + w0i * src103r;

                f32_4x tempW = w0r * wStepr - w0i * wStepi;
                w0i = w0r * wStepi + w0i * wStepr;
                w0r = tempW;

                f32_4x E01ri = F32_4x(EGrab);
                EGrab += 4;
                f32_4x E23ri = F32_4x(EGrab);
                EGrab += 4;

                f32_4x O01ri; O01ri.m = _mm_unpacklo_ps(O03r.m, O03i.m);
                f32_4x O23ri; O23ri.m = _mm_unpackhi_ps(O03r.m, O03i.m);

                f32_4x EpO01ri = E01ri + O01ri;
                f32_4x EmO01ri = E01ri - O01ri;
                f32_4x EpO23ri = E23ri + O23ri;
                f32_4x EmO23ri = E23ri - O23ri;

                _mm_store_ps(EPut, EpO01ri.m);
                EPut += 4;
                _mm_store_ps(OPut, EmO01ri.m);
                OPut += 4;
                _mm_store_ps(EPut, EpO23ri.m);
                EPut += 4;
                _mm_store_ps(OPut, EmO23ri.m);
                OPut += 4;

                f32_4x src145ri = F32_4x(OGrab);
                OGrab += 4;
                f32_4x src167ri = F32_4x(OGrab);
                OGrab += 4;

                f32_4x src147r; src147r.m = _mm_shuffle_ps(src145ri.m, src167ri.m, MULTILANE_SHUFFLE_MASK(0, 2, 0, 2));
                f32_4x src147i; src147i.m = _mm_shuffle_ps(src145ri.m, src167ri.m, MULTILANE_SHUFFLE_MASK(1, 3, 1, 3));

                f32_4x O47r = w1r * src147r - w1i * src147i;
                f32_4x O47i = w1r * src147i + w1i * src147r;

                f32_4x tempW2 = w1r * wStepr - w1i * wStepi;
                w1i = w1r * wStepi + w1i * wStepr;
                w1r = tempW2;

                f32_4x E45ri = F32_4x(EGrab);
                EGrab += 4;
                f32_4x E67ri = F32_4x(EGrab);
                EGrab += 4;

                f32_4x O45ri; O45ri.m = _mm_unpacklo_ps(O47r.m, O47i.m);
                f32_4x O67ri; O67ri.m = _mm_unpackhi_ps(O47r.m, O47i.m);

                f32_4x EpO45ri = E45ri + O45ri;
                f32_4x EmO45ri = E45ri - O45ri;
                f32_4x EpO67ri = E67ri + O67ri;
                f32_4x EmO67ri = E67ri - O67ri;

                _mm_store_ps(EPut, EpO45ri.m);
                EPut += 4;
                _mm_store_ps(OPut, EmO45ri.m);
                OPut += 4;
                _mm_store_ps(EPut, EpO67ri.m);
                EPut += 4;
                _mm_store_ps(OPut, EmO67ri.m);
                OPut += 4;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_inexact6(u32 count, c32 *signal, c32 *dest)
{
    copy(count * sizeof(c32), signal, dest);
    fft_inplace_inexact6(count, dest);
}

internal void
fft_normal2_64(u32 count, c64 *signal, c64 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < halfCount; index += 2)
    {
        u32 index0 = index + 0;
        u32 index1 = index + 1;
        u32 index2 = index0 + halfCount;
        u32 index3 = index1 + halfCount;

        u32 revIndex0 = reverse_bits32(index0, highBit.index);
        u32 revIndex1 = revIndex0 ^ halfCount;
        u32 revIndex2 = revIndex0 + 1;
        u32 revIndex3 = revIndex1 + 1;
        dest[index0] = signal[revIndex0];
        dest[index1] = signal[revIndex1];
        dest[index2] = signal[revIndex2];
        dest[index3] = signal[revIndex3];
    }

    // NOTE(michiel): w = e^(-i*2pi*j/m)
    SinCos64 cs4 = sincos64(-0.5 * F64_PI);
    c64 w4 = complex64(cs4.cos, cs4.sin);
    SinCos64 cs81 = sincos64(-0.25 * F64_PI);
    c64 w81 = complex64(cs81.cos, cs81.sin);
    SinCos64 cs82 = cs4;
    c64 w82 = complex64(cs82.cos, cs82.sin);
    SinCos64 cs83 = sincos64(-0.75 * F64_PI);
    c64 w83 = complex64(cs83.cos, cs83.sin);

    for (u32 k = 0; k < count; k += 8)
    {
        c64 *src0 = dest + k + 0;
        c64 *src1 = dest + k + 1;
        c64 *src2 = dest + k + 2;
        c64 *src3 = dest + k + 3;
        c64 *src4 = dest + k + 4;
        c64 *src5 = dest + k + 5;
        c64 *src6 = dest + k + 6;
        c64 *src7 = dest + k + 7;

        // NOTE(michiel): 2 -> 2
        c64 E20 = *src0;
        c64 O20 = *src1;
        *src0 = E20 + O20;
        *src1 = E20 - O20;
        c64 E21 = *src2;
        c64 O21 = *src3;
        *src2 = E21 + O21;
        *src3 = E21 - O21;
        c64 E22 = *src4;
        c64 O22 = *src5;
        *src4 = E22 + O22;
        *src5 = E22 - O22;
        c64 E23 = *src6;
        c64 O23 = *src7;
        *src6 = E23 + O23;
        *src7 = E23 - O23;

        // NOTE(michiel): 4 -> 4
        c64 E40 = *src0;
        c64 O40 = *src2;
        *src0 = E40 + O40;
        *src2 = E40 - O40;

        c64 E41 = *src1;
        c64 O41 = w4 * *src3;
        *src1 = E41 + O41;
        *src3 = E41 - O41;

        c64 E42 = *src4;
        c64 O42 = *src6;
        *src4 = E42 + O42;
        *src6 = E42 - O42;

        c64 E43 = *src5;
        c64 O43 = w4 * *src7;
        *src5 = E43 + O43;
        *src7 = E43 - O43;

        // NOTE(michiel): 8 -> 8
        c64 E80 = *src0;
        c64 O80 = *src4;
        *src0 = E80 + O80;
        *src4 = E80 - O80;

        c64 E81 = *src1;
        c64 O81 = w81 * *src5;
        *src1 = E81 + O81;
        *src5 = E81 - O81;

        c64 E82 = *src2;
        c64 O82 = w82 * *src6;
        *src2 = E82 + O82;
        *src6 = E82 - O82;

        c64 E83 = *src3;
        c64 O83 = w83 * *src7;
        *src3 = E83 + O83;
        *src7 = E83 - O83;
    }

    u32 halfM = 8;
    u32 m = 16;
    while (m <= count)
    {
        f64 oneOverM = (2.0 * F64_PI) / (f64)m;

        for (u32 k = 0; k < count; k += m)
        {
            c64 *src0 = dest + k;
            c64 *src1 = dest + k + halfM;
            for (u32 j = 0; j < halfM; j += 8)
            {
                // NOTE(michiel): w = e^(-i*2pi*j/m)
                SinCos64 cs0 = sincos64(-(f64)j * oneOverM);
                c64 w0 = complex64(cs0.cos, cs0.sin);
                SinCos64 cs1 = sincos64(-(f64)(j + 1) * oneOverM);
                c64 w1 = complex64(cs1.cos, cs1.sin);
                SinCos64 cs2 = sincos64(-(f64)(j + 2) * oneOverM);
                c64 w2 = complex64(cs2.cos, cs2.sin);
                SinCos64 cs3 = sincos64(-(f64)(j + 3) * oneOverM);
                c64 w3 = complex64(cs3.cos, cs3.sin);
                SinCos64 cs4 = sincos64(-(f64)(j + 4) * oneOverM);
                c64 w4 = complex64(cs4.cos, cs4.sin);
                SinCos64 cs5 = sincos64(-(f64)(j + 5) * oneOverM);
                c64 w5 = complex64(cs5.cos, cs5.sin);
                SinCos64 cs6 = sincos64(-(f64)(j + 6) * oneOverM);
                c64 w6 = complex64(cs6.cos, cs6.sin);
                SinCos64 cs7 = sincos64(-(f64)(j + 7) * oneOverM);
                c64 w7 = complex64(cs7.cos, cs7.sin);
                c64 E0 = *src0;
                c64 O0 = w0 * *src1;
                *src0++ = E0 + O0;
                *src1++ = E0 - O0;
                c64 E1 = *src0;
                c64 O1 = w1 * *src1;
                *src0++ = E1 + O1;
                *src1++ = E1 - O1;
                c64 E2 = *src0;
                c64 O2 = w2 * *src1;
                *src0++ = E2 + O2;
                *src1++ = E2 - O2;
                c64 E3 = *src0;
                c64 O3 = w3 * *src1;
                *src0++ = E3 + O3;
                *src1++ = E3 - O3;
                c64 E4 = *src0;
                c64 O4 = w4 * *src1;
                *src0++ = E4 + O4;
                *src1++ = E4 - O4;
                c64 E5 = *src0;
                c64 O5 = w5 * *src1;
                *src0++ = E5 + O5;
                *src1++ = E5 - O5;
                c64 E6 = *src0;
                c64 O6 = w6 * *src1;
                *src0++ = E6 + O6;
                *src1++ = E6 - O6;
                c64 E7 = *src0;
                c64 O7 = w7 * *src1;
                *src0++ = E7 + O7;
                *src1++ = E7 - O7;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_normal3_64(u32 count, c64 *signal, c64 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < halfCount; index += 2)
    {
        u32 index0 = index + 0;
        u32 index1 = index + 1;
        u32 index2 = index0 + halfCount;
        u32 index3 = index1 + halfCount;

        u32 revIndex0 = reverse_bits32(index0, highBit.index);
        u32 revIndex1 = revIndex0 ^ halfCount;
        u32 revIndex2 = revIndex0 + 1;
        u32 revIndex3 = revIndex1 + 1;
        dest[index0] = signal[revIndex0];
        dest[index1] = signal[revIndex1];
        dest[index2] = signal[revIndex2];
        dest[index3] = signal[revIndex3];
    }

    // NOTE(michiel): w = e^(-i*2pi*j/m)
    SinCos64 csBase0 = sincos64(0.0);
    SinCos64 csBase1 = sincos64(-0.25 * F64_PI);
    SinCos64 csBase2 = sincos64(-0.5 * F64_PI);
    SinCos64 csBase3 = sincos64(-0.75 * F64_PI);
    f64_2x wReal0 = F64_2x(csBase0.cos, csBase1.cos);
    f64_2x wReal1 = F64_2x(csBase2.cos, csBase3.cos);
    f64_2x wImag0 = F64_2x(csBase0.sin, csBase1.sin);
    f64_2x wImag1 = F64_2x(csBase2.sin, csBase3.sin);

    for (u32 k = 0; k < count; k += 8)
    {
        c64 src0 = *(dest + k + 0);
        c64 src1 = *(dest + k + 1);
        c64 src2 = *(dest + k + 2);
        c64 src3 = *(dest + k + 3);
        c64 src4 = *(dest + k + 4);
        c64 src5 = *(dest + k + 5);
        c64 src6 = *(dest + k + 6);
        c64 src7 = *(dest + k + 7);

        // NOTE(michiel): 2 -> 2
        f64 E20r = src0.real;
        f64 E20i = src0.imag;
        f64 O20r = src1.real;
        f64 O20i = src1.imag;

        f64 E21r = src2.real;
        f64 E21i = src2.imag;
        f64 O21r = src3.real;
        f64 O21i = src3.imag;

        f64 E22r = src4.real;
        f64 E22i = src4.imag;
        f64 O22r = src5.real;
        f64 O22i = src5.imag;

        f64 E23r = src6.real;
        f64 E23i = src6.imag;
        f64 O23r = src7.real;
        f64 O23i = src7.imag;

        src0.real = E20r + O20r;
        src0.imag = E20i + O20i;
        src1.real = E20r - O20r;
        src1.imag = E20i - O20i;

        src2.real = E21r + O21r;
        src2.imag = E21i + O21i;
        src3.real = E21r - O21r;
        src3.imag = E21i - O21i;

        src4.real = E22r + O22r;
        src4.imag = E22i + O22i;
        src5.real = E22r - O22r;
        src5.imag = E22i - O22i;

        src6.real = E23r + O23r;
        src6.imag = E23i + O23i;
        src7.real = E23r - O23r;
        src7.imag = E23i - O23i;

        // NOTE(michiel): 4 -> 4
        c64 E40 = src0;
        c64 E41 = src1;
        c64 O40 = src2;
        c64 O41 = complex64(wReal1.e[0], wImag1.e[0]) * src3;
        c64 E42 = src4;
        c64 E43 = src5;
        c64 O42 = src6;
        c64 O43 = complex64(wReal1.e[0], wImag1.e[0]) * src7;

        src0 = E40 + O40;
        src1 = E41 + O41;
        src2 = E40 - O40;
        src3 = E41 - O41;
        src4 = E42 + O42;
        src5 = E43 + O43;
        src6 = E42 - O42;
        src7 = E43 - O43;

        // NOTE(michiel): 8 -> 8
        c64 E80 = src0;
        c64 E81 = src1;
        c64 E82 = src2;
        c64 E83 = src3;
        c64 O80 = src4;
        c64 O81 = complex64(wReal0.e[1], wImag0.e[1]) * src5;
        c64 O82 = complex64(wReal1.e[0], wImag1.e[0]) * src6;
        c64 O83 = complex64(wReal1.e[1], wImag1.e[1]) * src7;

        src0 = E80 + O80;
        src1 = E81 + O81;
        src2 = E82 + O82;
        src3 = E83 + O83;
        src4 = E80 - O80;
        src5 = E81 - O81;
        src6 = E82 - O82;
        src7 = E83 - O83;

        *(dest + k + 0) = src0;
        *(dest + k + 1) = src1;
        *(dest + k + 2) = src2;
        *(dest + k + 3) = src3;
        *(dest + k + 4) = src4;
        *(dest + k + 5) = src5;
        *(dest + k + 6) = src6;
        *(dest + k + 7) = src7;
    }

    u32 halfM = 8;
    u32 m = 16;
    while (m <= count)
    {
        f64 oneOverM = (2.0 * F64_PI) / (f64)m;

        for (u32 k = 0; k < count; k += m)
        {
            c64 *src0 = dest + k;
            c64 *src1 = dest + k + halfM;
            for (u32 j = 0; j < halfM; j += 8)
            {
                SinCos64 cs00 = sincos64(-(f64)(j + 0) * oneOverM);
                SinCos64 cs01 = sincos64(-(f64)(j + 1) * oneOverM);
                SinCos64 cs02 = sincos64(-(f64)(j + 2) * oneOverM);
                SinCos64 cs03 = sincos64(-(f64)(j + 3) * oneOverM);
                SinCos64 cs10 = sincos64(-(f64)(j + 4) * oneOverM);
                SinCos64 cs11 = sincos64(-(f64)(j + 5) * oneOverM);
                SinCos64 cs12 = sincos64(-(f64)(j + 6) * oneOverM);
                SinCos64 cs13 = sincos64(-(f64)(j + 7) * oneOverM);

                f64_2x wReal00 = F64_2x(cs00.cos, cs01.cos);
                f64_2x wReal01 = F64_2x(cs02.cos, cs03.cos);
                f64_2x wImag00 = F64_2x(cs00.sin, cs01.sin);
                f64_2x wImag01 = F64_2x(cs02.sin, cs03.sin);
                f64_2x wReal10 = F64_2x(cs10.cos, cs11.cos);
                f64_2x wReal11 = F64_2x(cs12.cos, cs13.cos);
                f64_2x wImag10 = F64_2x(cs10.sin, cs11.sin);
                f64_2x wImag11 = F64_2x(cs12.sin, cs13.sin);
                c64 w0 = complex64(wReal00.e[0], wImag00.e[0]);
                c64 w1 = complex64(wReal00.e[1], wImag00.e[1]);
                c64 w2 = complex64(wReal01.e[0], wImag01.e[0]);
                c64 w3 = complex64(wReal01.e[1], wImag01.e[1]);
                c64 w4 = complex64(wReal10.e[0], wImag10.e[0]);
                c64 w5 = complex64(wReal10.e[1], wImag10.e[1]);
                c64 w6 = complex64(wReal11.e[0], wImag11.e[0]);
                c64 w7 = complex64(wReal11.e[1], wImag11.e[1]);

                c64 src00 = *(src0 + 0);
                c64 src01 = *(src0 + 1);
                c64 src02 = *(src0 + 2);
                c64 src03 = *(src0 + 3);
                c64 src04 = *(src0 + 4);
                c64 src05 = *(src0 + 5);
                c64 src06 = *(src0 + 6);
                c64 src07 = *(src0 + 7);

                c64 src10 = *(src1 + 0);
                c64 src11 = *(src1 + 1);
                c64 src12 = *(src1 + 2);
                c64 src13 = *(src1 + 3);
                c64 src14 = *(src1 + 4);
                c64 src15 = *(src1 + 5);
                c64 src16 = *(src1 + 6);
                c64 src17 = *(src1 + 7);

                c64 E0 = src00;
                c64 E1 = src01;
                c64 E2 = src02;
                c64 E3 = src03;
                c64 E4 = src04;
                c64 E5 = src05;
                c64 E6 = src06;
                c64 E7 = src07;
                c64 O0 = w0 * src10;
                c64 O1 = w1 * src11;
                c64 O2 = w2 * src12;
                c64 O3 = w3 * src13;
                c64 O4 = w4 * src14;
                c64 O5 = w5 * src15;
                c64 O6 = w6 * src16;
                c64 O7 = w7 * src17;

                src00 = E0 + O0;
                src01 = E1 + O1;
                src02 = E2 + O2;
                src03 = E3 + O3;
                src04 = E4 + O4;
                src05 = E5 + O5;
                src06 = E6 + O6;
                src07 = E7 + O7;
                src10 = E0 - O0;
                src11 = E1 - O1;
                src12 = E2 - O2;
                src13 = E3 - O3;
                src14 = E4 - O4;
                src15 = E5 - O5;
                src16 = E6 - O6;
                src17 = E7 - O7;

                *(src0 + 0) = src00;
                *(src0 + 1) = src01;
                *(src0 + 2) = src02;
                *(src0 + 3) = src03;
                *(src0 + 4) = src04;
                *(src0 + 5) = src05;
                *(src0 + 6) = src06;
                *(src0 + 7) = src07;

                *(src1 + 0) = src10;
                *(src1 + 1) = src11;
                *(src1 + 2) = src12;
                *(src1 + 3) = src13;
                *(src1 + 4) = src14;
                *(src1 + 5) = src15;
                *(src1 + 6) = src16;
                *(src1 + 7) = src17;

                src0 += 8;
                src1 += 8;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_normal4_64(u32 count, c64 *signal, c64 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < halfCount; index += 2)
    {
        u32 index0 = index + 0;
        u32 index1 = index + 1;
        u32 index2 = index0 + halfCount;
        u32 index3 = index1 + halfCount;

        u32 revIndex0 = reverse_bits32(index0, highBit.index);
        u32 revIndex1 = revIndex0 ^ halfCount;
        u32 revIndex2 = revIndex0 + 1;
        u32 revIndex3 = revIndex1 + 1;
        dest[index0] = signal[revIndex0];
        dest[index1] = signal[revIndex1];
        dest[index2] = signal[revIndex2];
        dest[index3] = signal[revIndex3];
    }

    // NOTE(michiel): w = e^(-i*2pi*j/m)
    SinCos64 csBase0 = sincos64(0.0);
    SinCos64 csBase1 = sincos64(-0.25 * F64_PI);
    SinCos64 csBase2 = sincos64(-0.5 * F64_PI);
    SinCos64 csBase3 = sincos64(-0.75 * F64_PI);
    f64_2x wReal0 = F64_2x(csBase0.cos, csBase1.cos);
    f64_2x wReal1 = F64_2x(csBase2.cos, csBase3.cos);
    f64_2x wImag0 = F64_2x(csBase0.sin, csBase1.sin);
    f64_2x wImag1 = F64_2x(csBase2.sin, csBase3.sin);

    for (u32 k = 0; k < count; k += 8)
    {
        c64 src0 = *(dest + k + 0);
        c64 src1 = *(dest + k + 1);
        c64 src2 = *(dest + k + 2);
        c64 src3 = *(dest + k + 3);
        c64 src4 = *(dest + k + 4);
        c64 src5 = *(dest + k + 5);
        c64 src6 = *(dest + k + 6);
        c64 src7 = *(dest + k + 7);

        // NOTE(michiel): 2 -> 2
        f64 E20r = src0.real;
        f64 E20i = src0.imag;
        f64 O20r = src1.real;
        f64 O20i = src1.imag;

        f64 E21r = src2.real;
        f64 E21i = src2.imag;
        f64 O21r = src3.real;
        f64 O21i = src3.imag;

        f64 E22r = src4.real;
        f64 E22i = src4.imag;
        f64 O22r = src5.real;
        f64 O22i = src5.imag;

        f64 E23r = src6.real;
        f64 E23i = src6.imag;
        f64 O23r = src7.real;
        f64 O23i = src7.imag;

        f64 o41hr = E21r - O21r;
        f64 o41hi = E21i - O21i;
        f64 o43hr = E23r - O23r;
        f64 o43hi = E23i - O23i;

        // NOTE(michiel): 4 -> 4
        f64 E40r = E20r + O20r;
        f64 E40i = E20i + O20i;
        f64 E41r = E20r - O20r;
        f64 E41i = E20i - O20i;

        f64 O40r = E21r + O21r;
        f64 O40i = E21i + O21i;
        f64 O41r = wReal1.e[0] * o41hr - wImag1.e[0] * o41hi;
        f64 O41i = wReal1.e[0] * o41hi + wImag1.e[0] * o41hr;

        f64 E42r = E22r + O22r;
        f64 E42i = E22i + O22i;
        f64 E43r = E22r - O22r;
        f64 E43i = E22i - O22i;

        f64 O42r = E23r + O23r;
        f64 O42i = E23i + O23i;
        f64 O43r = wReal1.e[0] * o43hr - wImag1.e[0] * o43hi;
        f64 O43i = wReal1.e[0] * o43hi + wImag1.e[0] * o43hr;

        f64 o80hr = E42r + O42r;
        f64 o80hi = E42i + O42i;
        f64 o81hr = E43r + O43r;
        f64 o81hi = E43i + O43i;
        f64 o82hr = E42r - O42r;
        f64 o82hi = E42i - O42i;
        f64 o83hr = E43r - O43r;
        f64 o83hi = E43i - O43i;

        // NOTE(michiel): 8 -> 8
        f64 E80r = E40r + O40r;
        f64 E80i = E40i + O40i;
        f64 E81r = E41r + O41r;
        f64 E81i = E41i + O41i;
        f64 E82r = E40r - O40r;
        f64 E82i = E40i - O40i;
        f64 E83r = E41r - O41r;
        f64 E83i = E41i - O41i;
        f64 O80r = wReal0.e[0] * o80hr - wImag0.e[0] * o80hi;
        f64 O80i = wReal0.e[0] * o80hi + wImag0.e[0] * o80hr;
        f64 O81r = wReal0.e[1] * o81hr - wImag0.e[1] * o81hi;
        f64 O81i = wReal0.e[1] * o81hi + wImag0.e[1] * o81hr;

        f64 O82r = wReal1.e[0] * o82hr - wImag1.e[0] * o82hi;
        f64 O82i = wReal1.e[0] * o82hi + wImag1.e[0] * o82hr;
        f64 O83r = wReal1.e[1] * o83hr - wImag1.e[1] * o83hi;
        f64 O83i = wReal1.e[1] * o83hi + wImag1.e[1] * o83hr;

        src0.real = E80r + O80r;
        src0.imag = E80i + O80i;
        src1.real = E81r + O81r;
        src1.imag = E81i + O81i;

        src2.real = E82r + O82r;
        src2.imag = E82i + O82i;
        src3.real = E83r + O83r;
        src3.imag = E83i + O83i;

        src4.real = E80r - O80r;
        src4.imag = E80i - O80i;
        src5.real = E81r - O81r;
        src5.imag = E81i - O81i;

        src6.real = E82r - O82r;
        src6.imag = E82i - O82i;
        src7.real = E83r - O83r;
        src7.imag = E83i - O83i;

        *(dest + k + 0) = src0;
        *(dest + k + 1) = src1;
        *(dest + k + 2) = src2;
        *(dest + k + 3) = src3;
        *(dest + k + 4) = src4;
        *(dest + k + 5) = src5;
        *(dest + k + 6) = src6;
        *(dest + k + 7) = src7;
    }

    u32 halfM = 8;
    u32 m = 16;
    while (m <= count)
    {
        f64 oneOverM = (2.0 * F64_PI) / (f64)m;

        for (u32 k = 0; k < count; k += m)
        {
            c64 *src0 = dest + k;
            c64 *src1 = dest + k + halfM;
            for (u32 j = 0; j < halfM; j += 8)
            {
                SinCos64 cs00 = sincos64(-(f64)(j + 0) * oneOverM);
                SinCos64 cs01 = sincos64(-(f64)(j + 1) * oneOverM);
                SinCos64 cs02 = sincos64(-(f64)(j + 2) * oneOverM);
                SinCos64 cs03 = sincos64(-(f64)(j + 3) * oneOverM);
                SinCos64 cs10 = sincos64(-(f64)(j + 4) * oneOverM);
                SinCos64 cs11 = sincos64(-(f64)(j + 5) * oneOverM);
                SinCos64 cs12 = sincos64(-(f64)(j + 6) * oneOverM);
                SinCos64 cs13 = sincos64(-(f64)(j + 7) * oneOverM);

                c64 src00 = *(src0 + 0);
                c64 src01 = *(src0 + 1);
                c64 src02 = *(src0 + 2);
                c64 src03 = *(src0 + 3);
                c64 src04 = *(src0 + 4);
                c64 src05 = *(src0 + 5);
                c64 src06 = *(src0 + 6);
                c64 src07 = *(src0 + 7);

                c64 src10 = *(src1 + 0);
                c64 src11 = *(src1 + 1);
                c64 src12 = *(src1 + 2);
                c64 src13 = *(src1 + 3);
                c64 src14 = *(src1 + 4);
                c64 src15 = *(src1 + 5);
                c64 src16 = *(src1 + 6);
                c64 src17 = *(src1 + 7);

                f64 E0r = src00.real;
                f64 E0i = src00.imag;
                f64 E1r = src01.real;
                f64 E1i = src01.imag;

                f64 E2r = src02.real;
                f64 E2i = src02.imag;
                f64 E3r = src03.real;
                f64 E3i = src03.imag;

                f64 E4r = src04.real;
                f64 E4i = src04.imag;
                f64 E5r = src05.real;
                f64 E5i = src05.imag;

                f64 E6r = src06.real;
                f64 E6i = src06.imag;
                f64 E7r = src07.real;
                f64 E7i = src07.imag;

                f64 O0r = cs00.cos * src10.real - cs00.sin * src10.imag;
                f64 O0i = cs00.cos * src10.imag + cs00.sin * src10.real;
                f64 O1r = cs01.cos * src11.real - cs01.sin * src11.imag;
                f64 O1i = cs01.cos * src11.imag + cs01.sin * src11.real;

                f64 O2r = cs02.cos * src12.real - cs02.sin * src12.imag;
                f64 O2i = cs02.cos * src12.imag + cs02.sin * src12.real;
                f64 O3r = cs03.cos * src13.real - cs03.sin * src13.imag;
                f64 O3i = cs03.cos * src13.imag + cs03.sin * src13.real;

                f64 O4r = cs10.cos * src14.real - cs10.sin * src14.imag;
                f64 O4i = cs10.cos * src14.imag + cs10.sin * src14.real;
                f64 O5r = cs11.cos * src15.real - cs11.sin * src15.imag;
                f64 O5i = cs11.cos * src15.imag + cs11.sin * src15.real;

                f64 O6r = cs12.cos * src16.real - cs12.sin * src16.imag;
                f64 O6i = cs12.cos * src16.imag + cs12.sin * src16.real;
                f64 O7r = cs13.cos * src17.real - cs13.sin * src17.imag;
                f64 O7i = cs13.cos * src17.imag + cs13.sin * src17.real;

                src00.real = E0r + O0r;
                src00.imag = E0i + O0i;
                src01.real = E1r + O1r;
                src01.imag = E1i + O1i;

                src02.real = E2r + O2r;
                src02.imag = E2i + O2i;
                src03.real = E3r + O3r;
                src03.imag = E3i + O3i;

                src04.real = E4r + O4r;
                src04.imag = E4i + O4i;
                src05.real = E5r + O5r;
                src05.imag = E5i + O5i;

                src06.real = E6r + O6r;
                src06.imag = E6i + O6i;
                src07.real = E7r + O7r;
                src07.imag = E7i + O7i;

                src10.real = E0r - O0r;
                src10.imag = E0i - O0i;
                src11.real = E1r - O1r;
                src11.imag = E1i - O1i;

                src12.real = E2r - O2r;
                src12.imag = E2i - O2i;
                src13.real = E3r - O3r;
                src13.imag = E3i - O3i;

                src14.real = E4r - O4r;
                src14.imag = E4i - O4i;
                src15.real = E5r - O5r;
                src15.imag = E5i - O5i;

                src16.real = E6r - O6r;
                src16.imag = E6i - O6i;
                src17.real = E7r - O7r;
                src17.imag = E7i - O7i;

                *(src0 + 0) = src00;
                *(src0 + 1) = src01;
                *(src0 + 2) = src02;
                *(src0 + 3) = src03;
                *(src0 + 4) = src04;
                *(src0 + 5) = src05;
                *(src0 + 6) = src06;
                *(src0 + 7) = src07;

                *(src1 + 0) = src10;
                *(src1 + 1) = src11;
                *(src1 + 2) = src12;
                *(src1 + 3) = src13;
                *(src1 + 4) = src14;
                *(src1 + 5) = src15;
                *(src1 + 6) = src16;
                *(src1 + 7) = src17;

                src0 += 8;
                src1 += 8;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_normal5_64(u32 count, c64 *signal, c64 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < halfCount; index += 2)
    {
        u32 index0 = index + 0;
        u32 index1 = index + 1;
        u32 index2 = index0 + halfCount;
        u32 index3 = index1 + halfCount;

        u32 revIndex0 = reverse_bits32(index0, highBit.index);
        u32 revIndex1 = revIndex0 ^ halfCount;
        u32 revIndex2 = revIndex0 + 1;
        u32 revIndex3 = revIndex1 + 1;
        dest[index0] = signal[revIndex0];
        dest[index1] = signal[revIndex1];
        dest[index2] = signal[revIndex2];
        dest[index3] = signal[revIndex3];
    }

    // NOTE(michiel): w = e^(-i*2pi*j/m)
    SinCos64 csBase0 = sincos64(0.0);
    SinCos64 csBase1 = sincos64(-0.25 * F64_PI);
    SinCos64 csBase2 = sincos64(-0.5 * F64_PI);
    SinCos64 csBase3 = sincos64(-0.75 * F64_PI);
    f64_2x cos42 = F64_2x(csBase2.cos, csBase2.cos);
    f64_2x sin42 = F64_2x(csBase2.sin, csBase2.sin);
    f64_2x cos80 = F64_2x(csBase0.cos, csBase0.cos);
    f64_2x cos81 = F64_2x(csBase1.cos, csBase1.cos);
    f64_2x cos82 = F64_2x(csBase2.cos, csBase2.cos);
    f64_2x cos83 = F64_2x(csBase3.cos, csBase3.cos);
    f64_2x sin80 = F64_2x(csBase0.sin, csBase0.sin);
    f64_2x sin81 = F64_2x(csBase1.sin, csBase1.sin);
    f64_2x sin82 = F64_2x(csBase2.sin, csBase2.sin);
    f64_2x sin83 = F64_2x(csBase3.sin, csBase3.sin);

    for (u32 k = 0; k < count; k += 8)
    {
        f64_2x E20ri = F64_2x((f64 *)(dest + k + 0));
        f64_2x O20ri = F64_2x((f64 *)(dest + k + 1));
        f64_2x E21ri = F64_2x((f64 *)(dest + k + 2));
        f64_2x O21ri = F64_2x((f64 *)(dest + k + 3));
        f64_2x E22ri = F64_2x((f64 *)(dest + k + 4));
        f64_2x O22ri = F64_2x((f64 *)(dest + k + 5));
        f64_2x E23ri = F64_2x((f64 *)(dest + k + 6));
        f64_2x O23ri = F64_2x((f64 *)(dest + k + 7));

        f64_2x E40ri = E20ri + O20ri;
        f64_2x E41ri = E20ri - O20ri;
        f64_2x E42ri = E22ri + O22ri;
        f64_2x E43ri = E22ri - O22ri;

        f64_2x O40ri = E21ri + O21ri;
        f64_2x o41hri = E21ri - O21ri;
        f64_2x o41hir; o41hir.md = _mm_shuffle_pd(o41hri.md, o41hri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O41ri; O41ri.md = _mm_addsub_pd((cos42 * o41hri).md, (sin42 * o41hir).md);
        f64_2x O42ri = E23ri + O23ri;
        f64_2x o43hri = E23ri - O23ri;
        f64_2x o43hir; o43hir.md = _mm_shuffle_pd(o43hri.md, o43hri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O43ri; O43ri.md = _mm_addsub_pd((cos42 * o43hri).md, (sin42 * o43hir).md);

        f64_2x EpO40ri = E40ri + O40ri;
        f64_2x EmO40ri = E40ri - O40ri;
        f64_2x EpO41ri = E41ri + O41ri;
        f64_2x EmO41ri = E41ri - O41ri;
        f64_2x EpO42ri = E42ri + O42ri;
        f64_2x EmO42ri = E42ri - O42ri;
        f64_2x EpO43ri = E43ri + O43ri;
        f64_2x EmO43ri = E43ri - O43ri;

        f64_2x E80ri = EpO40ri;
        f64_2x E81ri = EpO41ri;
        f64_2x E82ri = EmO40ri;
        f64_2x E83ri = EmO41ri;

        f64_2x o80hri = EpO42ri;
        f64_2x o80hir; o80hir.md = _mm_shuffle_pd(EpO42ri.md, EpO42ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O80ri; O80ri.md = _mm_addsub_pd((cos80 * o80hri).md, (sin80 * o80hir).md);
        f64_2x o81hri = EpO43ri;
        f64_2x o81hir; o81hir.md = _mm_shuffle_pd(EpO43ri.md, EpO43ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O81ri; O81ri.md = _mm_addsub_pd((cos81 * o81hri).md, (sin81 * o81hir).md);

        f64_2x o82hri = EmO42ri;
        f64_2x o82hir; o82hir.md = _mm_shuffle_pd(EmO42ri.md, EmO42ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O82ri; O82ri.md = _mm_addsub_pd((cos82 * o82hri).md, (sin82 * o82hir).md);
        f64_2x o83hri = EmO43ri;
        f64_2x o83hir; o83hir.md = _mm_shuffle_pd(EmO43ri.md, EmO43ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O83ri; O83ri.md = _mm_addsub_pd((cos83 * o83hri).md, (sin83 * o83hir).md);

        f64_2x EpO80 = E80ri + O80ri;
        f64_2x EpO81 = E81ri + O81ri;
        f64_2x EpO82 = E82ri + O82ri;
        f64_2x EpO83 = E83ri + O83ri;
        f64_2x EmO80 = E80ri - O80ri;
        f64_2x EmO81 = E81ri - O81ri;
        f64_2x EmO82 = E82ri - O82ri;
        f64_2x EmO83 = E83ri - O83ri;

        _mm_store_pd((f64 *)(dest + k + 0), EpO80.md);
        _mm_store_pd((f64 *)(dest + k + 1), EpO81.md);
        _mm_store_pd((f64 *)(dest + k + 2), EpO82.md);
        _mm_store_pd((f64 *)(dest + k + 3), EpO83.md);
        _mm_store_pd((f64 *)(dest + k + 4), EmO80.md);
        _mm_store_pd((f64 *)(dest + k + 5), EmO81.md);
        _mm_store_pd((f64 *)(dest + k + 6), EmO82.md);
        _mm_store_pd((f64 *)(dest + k + 7), EmO83.md);
    }

    u32 halfM = 8;
    u32 m = 16;
    while (m <= count)
    {
        f64 oneOverM = (2.0 * F64_PI) / (f64)m;

        for (u32 k = 0; k < count; k += m)
        {
            c64 *src0 = dest + k;
            c64 *src1 = dest + k + halfM;
            for (u32 j = 0; j < halfM; j += 8)
            {
                SinCos64 cs00 = sincos64(-(f64)(j + 0) * oneOverM);
                SinCos64 cs01 = sincos64(-(f64)(j + 1) * oneOverM);
                SinCos64 cs02 = sincos64(-(f64)(j + 2) * oneOverM);
                SinCos64 cs03 = sincos64(-(f64)(j + 3) * oneOverM);
                SinCos64 cs10 = sincos64(-(f64)(j + 4) * oneOverM);
                SinCos64 cs11 = sincos64(-(f64)(j + 5) * oneOverM);
                SinCos64 cs12 = sincos64(-(f64)(j + 6) * oneOverM);
                SinCos64 cs13 = sincos64(-(f64)(j + 7) * oneOverM);

                c64 src00 = *(src0 + 0);
                c64 src01 = *(src0 + 1);
                c64 src02 = *(src0 + 2);
                c64 src03 = *(src0 + 3);
                c64 src04 = *(src0 + 4);
                c64 src05 = *(src0 + 5);
                c64 src06 = *(src0 + 6);
                c64 src07 = *(src0 + 7);

                c64 src10 = *(src1 + 0);
                c64 src11 = *(src1 + 1);
                c64 src12 = *(src1 + 2);
                c64 src13 = *(src1 + 3);
                c64 src14 = *(src1 + 4);
                c64 src15 = *(src1 + 5);
                c64 src16 = *(src1 + 6);
                c64 src17 = *(src1 + 7);

                f64 E0r = src00.real;
                f64 E0i = src00.imag;
                f64 E1r = src01.real;
                f64 E1i = src01.imag;

                f64 E2r = src02.real;
                f64 E2i = src02.imag;
                f64 E3r = src03.real;
                f64 E3i = src03.imag;

                f64 E4r = src04.real;
                f64 E4i = src04.imag;
                f64 E5r = src05.real;
                f64 E5i = src05.imag;

                f64 E6r = src06.real;
                f64 E6i = src06.imag;
                f64 E7r = src07.real;
                f64 E7i = src07.imag;

                f64 O0r = cs00.cos * src10.real - cs00.sin * src10.imag;
                f64 O0i = cs00.cos * src10.imag + cs00.sin * src10.real;
                f64 O1r = cs01.cos * src11.real - cs01.sin * src11.imag;
                f64 O1i = cs01.cos * src11.imag + cs01.sin * src11.real;

                f64 O2r = cs02.cos * src12.real - cs02.sin * src12.imag;
                f64 O2i = cs02.cos * src12.imag + cs02.sin * src12.real;
                f64 O3r = cs03.cos * src13.real - cs03.sin * src13.imag;
                f64 O3i = cs03.cos * src13.imag + cs03.sin * src13.real;

                f64 O4r = cs10.cos * src14.real - cs10.sin * src14.imag;
                f64 O4i = cs10.cos * src14.imag + cs10.sin * src14.real;
                f64 O5r = cs11.cos * src15.real - cs11.sin * src15.imag;
                f64 O5i = cs11.cos * src15.imag + cs11.sin * src15.real;

                f64 O6r = cs12.cos * src16.real - cs12.sin * src16.imag;
                f64 O6i = cs12.cos * src16.imag + cs12.sin * src16.real;
                f64 O7r = cs13.cos * src17.real - cs13.sin * src17.imag;
                f64 O7i = cs13.cos * src17.imag + cs13.sin * src17.real;

                src00.real = E0r + O0r;
                src00.imag = E0i + O0i;
                src01.real = E1r + O1r;
                src01.imag = E1i + O1i;

                src02.real = E2r + O2r;
                src02.imag = E2i + O2i;
                src03.real = E3r + O3r;
                src03.imag = E3i + O3i;

                src04.real = E4r + O4r;
                src04.imag = E4i + O4i;
                src05.real = E5r + O5r;
                src05.imag = E5i + O5i;

                src06.real = E6r + O6r;
                src06.imag = E6i + O6i;
                src07.real = E7r + O7r;
                src07.imag = E7i + O7i;

                src10.real = E0r - O0r;
                src10.imag = E0i - O0i;
                src11.real = E1r - O1r;
                src11.imag = E1i - O1i;

                src12.real = E2r - O2r;
                src12.imag = E2i - O2i;
                src13.real = E3r - O3r;
                src13.imag = E3i - O3i;

                src14.real = E4r - O4r;
                src14.imag = E4i - O4i;
                src15.real = E5r - O5r;
                src15.imag = E5i - O5i;

                src16.real = E6r - O6r;
                src16.imag = E6i - O6i;
                src17.real = E7r - O7r;
                src17.imag = E7i - O7i;

                *(src0 + 0) = src00;
                *(src0 + 1) = src01;
                *(src0 + 2) = src02;
                *(src0 + 3) = src03;
                *(src0 + 4) = src04;
                *(src0 + 5) = src05;
                *(src0 + 6) = src06;
                *(src0 + 7) = src07;

                *(src1 + 0) = src10;
                *(src1 + 1) = src11;
                *(src1 + 2) = src12;
                *(src1 + 3) = src13;
                *(src1 + 4) = src14;
                *(src1 + 5) = src15;
                *(src1 + 6) = src16;
                *(src1 + 7) = src17;

                src0 += 8;
                src1 += 8;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_inplace6_64(u32 count, c64 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < halfCount; index += 2)
    {
        u32 index0 = index + 0;
        u32 index1 = index + 1;
        u32 index2 = index0 + halfCount;
        u32 index3 = index1 + halfCount;

        u32 revIndex0 = reverse_bits32(index0, highBit.index);
        u32 revIndex1 = revIndex0 ^ halfCount;
        u32 revIndex2 = revIndex0 + 1;
        u32 revIndex3 = revIndex1 + 1;

        if (revIndex0 > index0) {
            c64 temp = dest[index0];
            dest[index0] = dest[revIndex0];
            dest[revIndex0] = temp;
        }
        if (revIndex1 > index1) {
            c64 temp = dest[index1];
            dest[index1] = dest[revIndex1];
            dest[revIndex1] = temp;
        }
        if (revIndex2 > index2) {
            c64 temp = dest[index2];
            dest[index2] = dest[revIndex2];
            dest[revIndex2] = temp;
        }
        if (revIndex3 > index3) {
            c64 temp = dest[index3];
            dest[index3] = dest[revIndex3];
            dest[revIndex3] = temp;
        }
    }

    // NOTE(michiel): w = e^(-i*2pi*j/m)
    SinCos64 csBase0 = sincos64(0.0);
    SinCos64 csBase1 = sincos64(-0.25 * F64_PI);
    SinCos64 csBase2 = sincos64(-0.5 * F64_PI);
    SinCos64 csBase3 = sincos64(-0.75 * F64_PI);
    f64_2x cos42 = F64_2x(csBase2.cos, csBase2.cos);
    f64_2x sin42 = F64_2x(csBase2.sin, csBase2.sin);
    f64_2x cos80 = F64_2x(csBase0.cos, csBase0.cos);
    f64_2x cos81 = F64_2x(csBase1.cos, csBase1.cos);
    f64_2x cos82 = F64_2x(csBase2.cos, csBase2.cos);
    f64_2x cos83 = F64_2x(csBase3.cos, csBase3.cos);
    f64_2x sin80 = F64_2x(csBase0.sin, csBase0.sin);
    f64_2x sin81 = F64_2x(csBase1.sin, csBase1.sin);
    f64_2x sin82 = F64_2x(csBase2.sin, csBase2.sin);
    f64_2x sin83 = F64_2x(csBase3.sin, csBase3.sin);

    for (u32 k = 0; k < count; k += 8)
    {
        f64_2x E20ri = F64_2x((f64 *)(dest + k + 0));
        f64_2x O20ri = F64_2x((f64 *)(dest + k + 1));
        f64_2x E21ri = F64_2x((f64 *)(dest + k + 2));
        f64_2x O21ri = F64_2x((f64 *)(dest + k + 3));
        f64_2x E22ri = F64_2x((f64 *)(dest + k + 4));
        f64_2x O22ri = F64_2x((f64 *)(dest + k + 5));
        f64_2x E23ri = F64_2x((f64 *)(dest + k + 6));
        f64_2x O23ri = F64_2x((f64 *)(dest + k + 7));

        f64_2x E40ri = E20ri + O20ri;
        f64_2x E41ri = E20ri - O20ri;
        f64_2x E42ri = E22ri + O22ri;
        f64_2x E43ri = E22ri - O22ri;

        f64_2x O40ri = E21ri + O21ri;
        f64_2x o41hri = E21ri - O21ri;
        f64_2x o41hir; o41hir.md = _mm_shuffle_pd(o41hri.md, o41hri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O41ri; O41ri.md = _mm_addsub_pd((cos42 * o41hri).md, (sin42 * o41hir).md);
        f64_2x O42ri = E23ri + O23ri;
        f64_2x o43hri = E23ri - O23ri;
        f64_2x o43hir; o43hir.md = _mm_shuffle_pd(o43hri.md, o43hri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O43ri; O43ri.md = _mm_addsub_pd((cos42 * o43hri).md, (sin42 * o43hir).md);

        f64_2x EpO40ri = E40ri + O40ri;
        f64_2x EmO40ri = E40ri - O40ri;
        f64_2x EpO41ri = E41ri + O41ri;
        f64_2x EmO41ri = E41ri - O41ri;
        f64_2x EpO42ri = E42ri + O42ri;
        f64_2x EmO42ri = E42ri - O42ri;
        f64_2x EpO43ri = E43ri + O43ri;
        f64_2x EmO43ri = E43ri - O43ri;

        f64_2x o80hir; o80hir.md = _mm_shuffle_pd(EpO42ri.md, EpO42ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O80ri; O80ri.md = _mm_addsub_pd((cos80 * EpO42ri).md, (sin80 * o80hir).md);
        f64_2x o81hir; o81hir.md = _mm_shuffle_pd(EpO43ri.md, EpO43ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O81ri; O81ri.md = _mm_addsub_pd((cos81 * EpO43ri).md, (sin81 * o81hir).md);

        f64_2x o82hir; o82hir.md = _mm_shuffle_pd(EmO42ri.md, EmO42ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O82ri; O82ri.md = _mm_addsub_pd((cos82 * EmO42ri).md, (sin82 * o82hir).md);
        f64_2x o83hir; o83hir.md = _mm_shuffle_pd(EmO43ri.md, EmO43ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O83ri; O83ri.md = _mm_addsub_pd((cos83 * EmO43ri).md, (sin83 * o83hir).md);

        f64_2x EpO80 = EpO40ri + O80ri;
        f64_2x EpO81 = EpO41ri + O81ri;
        f64_2x EpO82 = EmO40ri + O82ri;
        f64_2x EpO83 = EmO41ri + O83ri;
        f64_2x EmO80 = EpO40ri - O80ri;
        f64_2x EmO81 = EpO41ri - O81ri;
        f64_2x EmO82 = EmO40ri - O82ri;
        f64_2x EmO83 = EmO41ri - O83ri;

        _mm_store_pd((f64 *)(dest + k + 0), EpO80.md);
        _mm_store_pd((f64 *)(dest + k + 1), EpO81.md);
        _mm_store_pd((f64 *)(dest + k + 2), EpO82.md);
        _mm_store_pd((f64 *)(dest + k + 3), EpO83.md);
        _mm_store_pd((f64 *)(dest + k + 4), EmO80.md);
        _mm_store_pd((f64 *)(dest + k + 5), EmO81.md);
        _mm_store_pd((f64 *)(dest + k + 6), EmO82.md);
        _mm_store_pd((f64 *)(dest + k + 7), EmO83.md);
    }

    u32 halfM = 8;
    u32 m = 16;
    while (m <= count)
    {
        f64 oneOverM = (2.0 * F64_PI) / (f64)m;

        for (u32 k = 0; k < count; k += m)
        {
            c64 *src0 = dest + k;
            c64 *src1 = dest + k + halfM;

            f64 *EGrab = (f64 *)src0;
            f64 *OGrab = (f64 *)src1;
            f64 *EPut  = (f64 *)src0;
            f64 *OPut  = (f64 *)src1;

            for (u32 j = 0; j < halfM; j += 8)
            {
                SinCos64 cs0 = sincos64(-(f64)(j + 0) * oneOverM);
                SinCos64 cs1 = sincos64(-(f64)(j + 1) * oneOverM);
                SinCos64 cs2 = sincos64(-(f64)(j + 2) * oneOverM);
                SinCos64 cs3 = sincos64(-(f64)(j + 3) * oneOverM);
                SinCos64 cs4 = sincos64(-(f64)(j + 4) * oneOverM);
                SinCos64 cs5 = sincos64(-(f64)(j + 5) * oneOverM);
                SinCos64 cs6 = sincos64(-(f64)(j + 6) * oneOverM);
                SinCos64 cs7 = sincos64(-(f64)(j + 7) * oneOverM);

                f64_2x oSrc0ri = F64_2x(OGrab);
                OGrab += 2;
                f64_2x oSrc1ri = F64_2x(OGrab);
                OGrab += 2;
                f64_2x oSrc2ri = F64_2x(OGrab);
                OGrab += 2;
                f64_2x oSrc3ri = F64_2x(OGrab);
                OGrab += 2;
                f64_2x oSrc4ri = F64_2x(OGrab);
                OGrab += 2;
                f64_2x oSrc5ri = F64_2x(OGrab);
                OGrab += 2;
                f64_2x oSrc6ri = F64_2x(OGrab);
                OGrab += 2;
                f64_2x oSrc7ri = F64_2x(OGrab);
                OGrab += 2;

                f64_2x oSrc0ir; oSrc0ir.md = _mm_shuffle_pd(oSrc0ri.md, oSrc0ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x ocs0ri = F64_2x(cs0.cos) * oSrc0ri;
                f64_2x osc0ir = F64_2x(cs0.sin) * oSrc0ir;
                f64_2x o0ri; o0ri.md = _mm_addsub_pd(ocs0ri.md, osc0ir.md);

                f64_2x oSrc1ir; oSrc1ir.md = _mm_shuffle_pd(oSrc1ri.md, oSrc1ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x ocs1ri = F64_2x(cs1.cos) * oSrc1ri;
                f64_2x osc1ir = F64_2x(cs1.sin) * oSrc1ir;
                f64_2x o1ri; o1ri.md = _mm_addsub_pd(ocs1ri.md, osc1ir.md);

                f64_2x oSrc2ir; oSrc2ir.md = _mm_shuffle_pd(oSrc2ri.md, oSrc2ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x ocs2ri = F64_2x(cs2.cos) * oSrc2ri;
                f64_2x osc2ir = F64_2x(cs2.sin) * oSrc2ir;
                f64_2x o2ri; o2ri.md = _mm_addsub_pd(ocs2ri.md, osc2ir.md);

                f64_2x oSrc3ir; oSrc3ir.md = _mm_shuffle_pd(oSrc3ri.md, oSrc3ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x ocs3ri = F64_2x(cs3.cos) * oSrc3ri;
                f64_2x osc3ir = F64_2x(cs3.sin) * oSrc3ir;
                f64_2x o3ri; o3ri.md = _mm_addsub_pd(ocs3ri.md, osc3ir.md);

                f64_2x oSrc4ir; oSrc4ir.md = _mm_shuffle_pd(oSrc4ri.md, oSrc4ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x ocs4ri = F64_2x(cs4.cos) * oSrc4ri;
                f64_2x osc4ir = F64_2x(cs4.sin) * oSrc4ir;
                f64_2x o4ri; o4ri.md = _mm_addsub_pd(ocs4ri.md, osc4ir.md);

                f64_2x oSrc5ir; oSrc5ir.md = _mm_shuffle_pd(oSrc5ri.md, oSrc5ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x ocs5ri = F64_2x(cs5.cos) * oSrc5ri;
                f64_2x osc5ir = F64_2x(cs5.sin) * oSrc5ir;
                f64_2x o5ri; o5ri.md = _mm_addsub_pd(ocs5ri.md, osc5ir.md);

                f64_2x oSrc6ir; oSrc6ir.md = _mm_shuffle_pd(oSrc6ri.md, oSrc6ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x ocs6ri = F64_2x(cs6.cos) * oSrc6ri;
                f64_2x osc6ir = F64_2x(cs6.sin) * oSrc6ir;
                f64_2x o6ri; o6ri.md = _mm_addsub_pd(ocs6ri.md, osc6ir.md);

                f64_2x oSrc7ir; oSrc7ir.md = _mm_shuffle_pd(oSrc7ri.md, oSrc7ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x ocs7ri = F64_2x(cs7.cos) * oSrc7ri;
                f64_2x osc7ir = F64_2x(cs7.sin) * oSrc7ir;
                f64_2x o7ri; o7ri.md = _mm_addsub_pd(ocs7ri.md, osc7ir.md);

                f64_2x e0ri = F64_2x(EGrab);
                EGrab += 2;
                f64_2x e1ri = F64_2x(EGrab);
                EGrab += 2;
                f64_2x e2ri = F64_2x(EGrab);
                EGrab += 2;
                f64_2x e3ri = F64_2x(EGrab);
                EGrab += 2;
                f64_2x e4ri = F64_2x(EGrab);
                EGrab += 2;
                f64_2x e5ri = F64_2x(EGrab);
                EGrab += 2;
                f64_2x e6ri = F64_2x(EGrab);
                EGrab += 2;
                f64_2x e7ri = F64_2x(EGrab);
                EGrab += 2;

                f64_2x s00ri = e0ri + o0ri;
                f64_2x s10ri = e0ri - o0ri;

                f64_2x s01ri = e1ri + o1ri;
                f64_2x s11ri = e1ri - o1ri;

                f64_2x s02ri = e2ri + o2ri;
                f64_2x s12ri = e2ri - o2ri;

                f64_2x s03ri = e3ri + o3ri;
                f64_2x s13ri = e3ri - o3ri;

                f64_2x s04ri = e4ri + o4ri;
                f64_2x s14ri = e4ri - o4ri;

                f64_2x s05ri = e5ri + o5ri;
                f64_2x s15ri = e5ri - o5ri;

                f64_2x s06ri = e6ri + o6ri;
                f64_2x s16ri = e6ri - o6ri;

                f64_2x s07ri = e7ri + o7ri;
                f64_2x s17ri = e7ri - o7ri;

                _mm_store_pd(EPut, s00ri.md);
                EPut += 2;
                _mm_store_pd(OPut, s10ri.md);
                OPut += 2;
                _mm_store_pd(EPut, s01ri.md);
                EPut += 2;
                _mm_store_pd(OPut, s11ri.md);
                OPut += 2;
                _mm_store_pd(EPut, s02ri.md);
                EPut += 2;
                _mm_store_pd(OPut, s12ri.md);
                OPut += 2;
                _mm_store_pd(EPut, s03ri.md);
                EPut += 2;
                _mm_store_pd(OPut, s13ri.md);
                OPut += 2;
                _mm_store_pd(EPut, s04ri.md);
                EPut += 2;
                _mm_store_pd(OPut, s14ri.md);
                OPut += 2;
                _mm_store_pd(EPut, s05ri.md);
                EPut += 2;
                _mm_store_pd(OPut, s15ri.md);
                OPut += 2;
                _mm_store_pd(EPut, s06ri.md);
                EPut += 2;
                _mm_store_pd(OPut, s16ri.md);
                OPut += 2;
                _mm_store_pd(EPut, s07ri.md);
                EPut += 2;
                _mm_store_pd(OPut, s17ri.md);
                OPut += 2;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_inplace_inexact6_64(u32 count, c64 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < halfCount; index += 2)
    {
        u32 index0 = index + 0;
        u32 index1 = index + 1;
        u32 index2 = index0 + halfCount;
        u32 index3 = index1 + halfCount;

        u32 revIndex0 = reverse_bits32(index0, highBit.index);
        u32 revIndex1 = revIndex0 ^ halfCount;
        u32 revIndex2 = revIndex0 + 1;
        u32 revIndex3 = revIndex1 + 1;

        if (revIndex0 > index0) {
            c64 temp = dest[index0];
            dest[index0] = dest[revIndex0];
            dest[revIndex0] = temp;
        }
        if (revIndex1 > index1) {
            c64 temp = dest[index1];
            dest[index1] = dest[revIndex1];
            dest[revIndex1] = temp;
        }
        if (revIndex2 > index2) {
            c64 temp = dest[index2];
            dest[index2] = dest[revIndex2];
            dest[revIndex2] = temp;
        }
        if (revIndex3 > index3) {
            c64 temp = dest[index3];
            dest[index3] = dest[revIndex3];
            dest[revIndex3] = temp;
        }
    }

    // NOTE(michiel): w = e^(-i*2pi*j/m)
    SinCos64 csBase0 = sincos64(0.0);
    SinCos64 csBase1 = sincos64(-0.25 * F64_PI);
    SinCos64 csBase2 = sincos64(-0.5 * F64_PI);
    SinCos64 csBase3 = sincos64(-0.75 * F64_PI);
    f64_2x cos42 = F64_2x(csBase2.cos, csBase2.cos);
    f64_2x sin42 = F64_2x(csBase2.sin, csBase2.sin);
    f64_2x cos80 = F64_2x(csBase0.cos, csBase0.cos);
    f64_2x cos81 = F64_2x(csBase1.cos, csBase1.cos);
    f64_2x cos82 = F64_2x(csBase2.cos, csBase2.cos);
    f64_2x cos83 = F64_2x(csBase3.cos, csBase3.cos);
    f64_2x sin80 = F64_2x(csBase0.sin, csBase0.sin);
    f64_2x sin81 = F64_2x(csBase1.sin, csBase1.sin);
    f64_2x sin82 = F64_2x(csBase2.sin, csBase2.sin);
    f64_2x sin83 = F64_2x(csBase3.sin, csBase3.sin);

    for (u32 k = 0; k < count; k += 8)
    {
        f64_2x E20ri = F64_2x((f64 *)(dest + k + 0));
        f64_2x O20ri = F64_2x((f64 *)(dest + k + 1));
        f64_2x E21ri = F64_2x((f64 *)(dest + k + 2));
        f64_2x O21ri = F64_2x((f64 *)(dest + k + 3));
        f64_2x E22ri = F64_2x((f64 *)(dest + k + 4));
        f64_2x O22ri = F64_2x((f64 *)(dest + k + 5));
        f64_2x E23ri = F64_2x((f64 *)(dest + k + 6));
        f64_2x O23ri = F64_2x((f64 *)(dest + k + 7));

        f64_2x E40ri = E20ri + O20ri;
        f64_2x E41ri = E20ri - O20ri;
        f64_2x E42ri = E22ri + O22ri;
        f64_2x E43ri = E22ri - O22ri;

        f64_2x O40ri = E21ri + O21ri;
        f64_2x o41hri = E21ri - O21ri;
        f64_2x o41hir; o41hir.md = _mm_shuffle_pd(o41hri.md, o41hri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O41ri; O41ri.md = _mm_addsub_pd((cos42 * o41hri).md, (sin42 * o41hir).md);
        f64_2x O42ri = E23ri + O23ri;
        f64_2x o43hri = E23ri - O23ri;
        f64_2x o43hir; o43hir.md = _mm_shuffle_pd(o43hri.md, o43hri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O43ri; O43ri.md = _mm_addsub_pd((cos42 * o43hri).md, (sin42 * o43hir).md);

        f64_2x EpO40ri = E40ri + O40ri;
        f64_2x EmO40ri = E40ri - O40ri;
        f64_2x EpO41ri = E41ri + O41ri;
        f64_2x EmO41ri = E41ri - O41ri;
        f64_2x EpO42ri = E42ri + O42ri;
        f64_2x EmO42ri = E42ri - O42ri;
        f64_2x EpO43ri = E43ri + O43ri;
        f64_2x EmO43ri = E43ri - O43ri;

        f64_2x o80hir; o80hir.md = _mm_shuffle_pd(EpO42ri.md, EpO42ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O80ri; O80ri.md = _mm_addsub_pd((cos80 * EpO42ri).md, (sin80 * o80hir).md);
        f64_2x o81hir; o81hir.md = _mm_shuffle_pd(EpO43ri.md, EpO43ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O81ri; O81ri.md = _mm_addsub_pd((cos81 * EpO43ri).md, (sin81 * o81hir).md);

        f64_2x o82hir; o82hir.md = _mm_shuffle_pd(EmO42ri.md, EmO42ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O82ri; O82ri.md = _mm_addsub_pd((cos82 * EmO42ri).md, (sin82 * o82hir).md);
        f64_2x o83hir; o83hir.md = _mm_shuffle_pd(EmO43ri.md, EmO43ri.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
        f64_2x O83ri; O83ri.md = _mm_addsub_pd((cos83 * EmO43ri).md, (sin83 * o83hir).md);

        f64_2x EpO80 = EpO40ri + O80ri;
        f64_2x EpO81 = EpO41ri + O81ri;
        f64_2x EpO82 = EmO40ri + O82ri;
        f64_2x EpO83 = EmO41ri + O83ri;
        f64_2x EmO80 = EpO40ri - O80ri;
        f64_2x EmO81 = EpO41ri - O81ri;
        f64_2x EmO82 = EmO40ri - O82ri;
        f64_2x EmO83 = EmO41ri - O83ri;

        _mm_store_pd((f64 *)(dest + k + 0), EpO80.md);
        _mm_store_pd((f64 *)(dest + k + 1), EpO81.md);
        _mm_store_pd((f64 *)(dest + k + 2), EpO82.md);
        _mm_store_pd((f64 *)(dest + k + 3), EpO83.md);
        _mm_store_pd((f64 *)(dest + k + 4), EmO80.md);
        _mm_store_pd((f64 *)(dest + k + 5), EmO81.md);
        _mm_store_pd((f64 *)(dest + k + 6), EmO82.md);
        _mm_store_pd((f64 *)(dest + k + 7), EmO83.md);
    }

    u32 halfM = 8;
    u32 m = 16;

    f64 oneOverMpre = (4.0 * F64_PI) / (f64)m;
    SinCos64 csPre0 = sincos64(-oneOverMpre);
    SinCos64 csPre1 = sincos64(-oneOverMpre*2.0);
    SinCos64 csPre2 = sincos64(-oneOverMpre*3.0);
    SinCos64 csPre3 = sincos64(-oneOverMpre*4.0);

    c64 wm1 = complex64(csPre0.cos, csPre0.sin);
    c64 wm2 = complex64(csPre1.cos, csPre1.sin);
    c64 wm3 = complex64(csPre2.cos, csPre2.sin);
    c64 wm4 = complex64(csPre3.cos, csPre3.sin);

    while (m <= count)
    {
        f64 oneOverM = 0.5 * oneOverMpre;
        oneOverMpre = oneOverM;

        c64 wm8 = wm4;
        c64 wm6 = wm3;
        wm4 = wm2;
        wm2 = wm1;

        SinCos64 csLoop0 = sincos64(-oneOverM);
        SinCos64 csLoop1 = sincos64(-oneOverM*3.0);
        SinCos64 csLoop2 = sincos64(-oneOverM*5.0);
        SinCos64 csLoop3 = sincos64(-oneOverM*7.0);

        wm1 = complex64(csLoop0.cos, csLoop0.sin);
        wm3 = complex64(csLoop1.cos, csLoop1.sin);
        c64 wm5 = complex64(csLoop2.cos, csLoop2.sin);
        c64 wm7 = complex64(csLoop3.cos, csLoop3.sin);

        f64_2x wStepr = F64_2x(wm8.real);
        f64_2x wStepi = F64_2x(wm8.imag);

        f64_2x w0Startr = F64_2x(1.0, wm1.real);
        f64_2x w0Starti = F64_2x(0.0, wm1.imag);
        f64_2x w1Startr = F64_2x(wm2.real, wm3.real);
        f64_2x w1Starti = F64_2x(wm2.imag, wm3.imag);
        f64_2x w2Startr = F64_2x(wm4.real, wm5.real);
        f64_2x w2Starti = F64_2x(wm4.imag, wm5.imag);
        f64_2x w3Startr = F64_2x(wm6.real, wm7.real);
        f64_2x w3Starti = F64_2x(wm6.imag, wm7.imag);

        for (u32 k = 0; k < count; k += m)
        {
            c64 *src0 = dest + k;
            c64 *src1 = dest + k + halfM;

            f64_2x w0r = w0Startr;
            f64_2x w0i = w0Starti;
            f64_2x w1r = w1Startr;
            f64_2x w1i = w1Starti;
            f64_2x w2r = w2Startr;
            f64_2x w2i = w2Starti;
            f64_2x w3r = w3Startr;
            f64_2x w3i = w3Starti;

            f64 *EGrab = (f64 *)src0;
            f64 *OGrab = (f64 *)src1;
            f64 *EPut  = (f64 *)src0;
            f64 *OPut  = (f64 *)src1;

            for (u32 j = 0; j < halfM; j += 8)
            {
                f64_2x oSrc01 = F64_2x(OGrab);
                f64_2x eSrc01 = F64_2x(EGrab);
                OGrab += 2;
                EGrab += 2;

                f64_2x oSrc23 = F64_2x(OGrab);
                f64_2x eSrc23 = F64_2x(EGrab);
                OGrab += 2;
                EGrab += 2;

                f64_2x w0r0; w0r0.md = _mm_shuffle_pd(w0r.md, w0r.md, MULTILANE_SHUFFLE_MASK_D(0, 0));
                f64_2x w0i0; w0i0.md = _mm_shuffle_pd(w0i.md, w0i.md, MULTILANE_SHUFFLE_MASK_D(0, 0));
                f64_2x src10; src10.md = _mm_shuffle_pd(oSrc01.md, oSrc01.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x O01ri; O01ri.md = _mm_addsub_pd((oSrc01*w0r0).md, (src10*w0i0).md);

                f64_2x w0r1; w0r1.md = _mm_shuffle_pd(w0r.md, w0r.md, MULTILANE_SHUFFLE_MASK_D(1, 1));
                f64_2x w0i1; w0i1.md = _mm_shuffle_pd(w0i.md, w0i.md, MULTILANE_SHUFFLE_MASK_D(1, 1));
                f64_2x src32; src32.md = _mm_shuffle_pd(oSrc23.md, oSrc23.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x O23ri; O23ri.md = _mm_addsub_pd((oSrc23*w0r1).md, (src32*w0i1).md);

                _mm_store_pd(EPut, (eSrc01 + O01ri).md);
                EPut += 2;
                _mm_store_pd(OPut, (eSrc01 - O01ri).md);
                OPut += 2;

                _mm_store_pd(EPut, (eSrc23 + O23ri).md);
                EPut += 2;
                _mm_store_pd(OPut, (eSrc23 - O23ri).md);
                OPut += 2;

                f64_2x tempW0 = w0r * wStepr - w0i * wStepi;
                w0i = w0r * wStepi + w0i * wStepr;
                w0r = tempW0;

                f64_2x oSrc45 = F64_2x(OGrab);
                f64_2x eSrc45 = F64_2x(EGrab);
                OGrab += 2;
                EGrab += 2;

                f64_2x oSrc67 = F64_2x(OGrab);
                f64_2x eSrc67 = F64_2x(EGrab);
                OGrab += 2;
                EGrab += 2;

                f64_2x w1r0; w1r0.md = _mm_shuffle_pd(w1r.md, w1r.md, MULTILANE_SHUFFLE_MASK_D(0, 0));
                f64_2x w1i0; w1i0.md = _mm_shuffle_pd(w1i.md, w1i.md, MULTILANE_SHUFFLE_MASK_D(0, 0));
                f64_2x src54; src54.md = _mm_shuffle_pd(oSrc45.md, oSrc45.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x O45ri; O45ri.md = _mm_addsub_pd((oSrc45*w1r0).md, (src54*w1i0).md);

                f64_2x w1r1; w1r1.md = _mm_shuffle_pd(w1r.md, w1r.md, MULTILANE_SHUFFLE_MASK_D(1, 1));
                f64_2x w1i1; w1i1.md = _mm_shuffle_pd(w1i.md, w1i.md, MULTILANE_SHUFFLE_MASK_D(1, 1));
                f64_2x src76; src76.md = _mm_shuffle_pd(oSrc67.md, oSrc67.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x O67ri; O67ri.md = _mm_addsub_pd((oSrc67*w1r1).md, (src76*w1i1).md);

                _mm_store_pd(EPut, (eSrc45 + O45ri).md);
                EPut += 2;
                _mm_store_pd(OPut, (eSrc45 - O45ri).md);
                OPut += 2;

                _mm_store_pd(EPut, (eSrc67 + O67ri).md);
                EPut += 2;
                _mm_store_pd(OPut, (eSrc67 - O67ri).md);
                OPut += 2;

                f64_2x tempW1 = w1r * wStepr - w1i * wStepi;
                w1i = w1r * wStepi + w1i * wStepr;
                w1r = tempW1;

                f64_2x oSrc89 = F64_2x(OGrab);
                f64_2x eSrc89 = F64_2x(EGrab);
                OGrab += 2;
                EGrab += 2;

                f64_2x oSrc1011 = F64_2x(OGrab);
                f64_2x eSrc1011 = F64_2x(EGrab);
                OGrab += 2;
                EGrab += 2;

                f64_2x w2r0; w2r0.md = _mm_shuffle_pd(w2r.md, w2r.md, MULTILANE_SHUFFLE_MASK_D(0, 0));
                f64_2x w2i0; w2i0.md = _mm_shuffle_pd(w2i.md, w2i.md, MULTILANE_SHUFFLE_MASK_D(0, 0));
                f64_2x src98; src98.md = _mm_shuffle_pd(oSrc89.md, oSrc89.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x O89ri; O89ri.md = _mm_addsub_pd((oSrc89*w2r0).md, (src98*w2i0).md);

                f64_2x w2r1; w2r1.md = _mm_shuffle_pd(w2r.md, w2r.md, MULTILANE_SHUFFLE_MASK_D(1, 1));
                f64_2x w2i1; w2i1.md = _mm_shuffle_pd(w2i.md, w2i.md, MULTILANE_SHUFFLE_MASK_D(1, 1));
                f64_2x src1110; src1110.md = _mm_shuffle_pd(oSrc1011.md, oSrc1011.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x O1011ri; O1011ri.md = _mm_addsub_pd((oSrc1011*w2r1).md, (src1110*w2i1).md);

                _mm_store_pd(EPut, (eSrc89 + O89ri).md);
                EPut += 2;
                _mm_store_pd(OPut, (eSrc89 - O89ri).md);
                OPut += 2;

                _mm_store_pd(EPut, (eSrc1011 + O1011ri).md);
                EPut += 2;
                _mm_store_pd(OPut, (eSrc1011 - O1011ri).md);
                OPut += 2;

                f64_2x tempW2 = w2r * wStepr - w2i * wStepi;
                w2i = w2r * wStepi + w2i * wStepr;
                w2r = tempW2;

                f64_2x oSrc1213 = F64_2x(OGrab);
                f64_2x eSrc1213 = F64_2x(EGrab);
                OGrab += 2;
                EGrab += 2;

                f64_2x oSrc1415 = F64_2x(OGrab);
                f64_2x eSrc1415 = F64_2x(EGrab);
                OGrab += 2;
                EGrab += 2;

                f64_2x w3r0; w3r0.md = _mm_shuffle_pd(w3r.md, w3r.md, MULTILANE_SHUFFLE_MASK_D(0, 0));
                f64_2x w3i0; w3i0.md = _mm_shuffle_pd(w3i.md, w3i.md, MULTILANE_SHUFFLE_MASK_D(0, 0));
                f64_2x src1312; src1312.md = _mm_shuffle_pd(oSrc1213.md, oSrc1213.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x O1213ri; O1213ri.md = _mm_addsub_pd((oSrc1213*w3r0).md, (src1312*w3i0).md);

                f64_2x w3r1; w3r1.md = _mm_shuffle_pd(w3r.md, w3r.md, MULTILANE_SHUFFLE_MASK_D(1, 1));
                f64_2x w3i1; w3i1.md = _mm_shuffle_pd(w3i.md, w3i.md, MULTILANE_SHUFFLE_MASK_D(1, 1));
                f64_2x src1514; src1514.md = _mm_shuffle_pd(oSrc1415.md, oSrc1415.md, MULTILANE_SHUFFLE_MASK_D(1, 0));
                f64_2x O1415ri; O1415ri.md = _mm_addsub_pd((oSrc1415*w3r1).md, (src1514*w3i1).md);

                _mm_store_pd(EPut, (eSrc1213 + O1213ri).md);
                EPut += 2;
                _mm_store_pd(OPut, (eSrc1213 - O1213ri).md);
                OPut += 2;

                _mm_store_pd(EPut, (eSrc1415 + O1415ri).md);
                EPut += 2;
                _mm_store_pd(OPut, (eSrc1415 - O1415ri).md);
                OPut += 2;

                f64_2x tempW3 = w3r * wStepr - w3i * wStepi;
                w3i = w3r * wStepi + w3i * wStepr;
                w3r = tempW3;
            }
        }
        halfM = m;
        m <<= 1;
    }
}
