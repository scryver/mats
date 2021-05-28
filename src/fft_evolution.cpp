
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
    SinCos_4x csBase = sincos32_4x(F32_4x(0.0f, -0.25f * F32_PI, -0.5f * F32_PI, -0.75f * F32_PI));
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
                SinCos_4x cs0 = sincos32_4x(F32_4x(-(f32)(j + 0) * oneOverM,
                                                   -(f32)(j + 1) * oneOverM,
                                                   -(f32)(j + 2) * oneOverM,
                                                   -(f32)(j + 3) * oneOverM));
                SinCos_4x cs1 = sincos32_4x(F32_4x(-(f32)(j + 4) * oneOverM,
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
    SinCos_4x csBase = sincos32_4x(F32_4x(0.0f, -0.25f * F32_PI, -0.5f * F32_PI, -0.75f * F32_PI));
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
                SinCos_4x cs0 = sincos32_4x(F32_4x(-(f32)(j + 0) * oneOverM,
                                                   -(f32)(j + 1) * oneOverM,
                                                   -(f32)(j + 2) * oneOverM,
                                                   -(f32)(j + 3) * oneOverM));
                SinCos_4x cs1 = sincos32_4x(F32_4x(-(f32)(j + 4) * oneOverM,
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
    SinCos_4x csBase = sincos32_4x(F32_4x(0.0f, -0.25f * F32_PI, -0.5f * F32_PI, -0.75f * F32_PI));
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
                SinCos_4x cs0 = sincos32_4x(F32_4x(-(f32)(j + 0) * oneOverM,
                                                   -(f32)(j + 1) * oneOverM,
                                                   -(f32)(j + 2) * oneOverM,
                                                   -(f32)(j + 3) * oneOverM));
                SinCos_4x cs1 = sincos32_4x(F32_4x(-(f32)(j + 4) * oneOverM,
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
    SinCos_4x csBase = sincos32_4x(F32_4x(0.0f, -0.25f * F32_PI, -0.5f * F32_PI, -0.75f * F32_PI));
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

                SinCos_4x cs0 = sincos32_4x(angles);
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

                SinCos_4x cs1 = sincos32_4x(angles);
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
    SinCos_4x csBase = sincos32_4x(F32_4x(0.0f, -0.25f * F32_PI, -0.5f * F32_PI, -0.75f * F32_PI));
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
    SinCos_4x sinCosPre = sincos32_4x(F32_4x(-oneOverMpre, -oneOverMpre*2.0f, -oneOverMpre*3.0f, -oneOverMpre*4.0f));

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

        SinCos_4x sinCos = sincos32_4x(F32_4x(-oneOverM, -oneOverM*3.0f, -oneOverM*5.0f, -oneOverM*7.0f));
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
