#include "../libberdip/src/common.h"
#include "../libberdip/src/intrinsics.h"
#include "../libberdip/src/multilane.h"

#include <unistd.h>
#include <time.h>

#define MATS_USE_SSE4 1
#include "../mats/mats.h"
#include "../mats/mats_wide.h"

#include "../mats/matc.h"
#include "../mats/mats_fft.h"

internal f32
square(f32 x)
{
    return x * x;
}

internal f64
square(f64 x)
{
    return x * x;
}

internal f32
cos_pi(f32 x)
{
    return cos32(x);
}

internal f32
sin_pi(f32 x)
{
    return sin32(x);
}

typedef c32 Complex32;
typedef c64 Complex64;
#include "../libberdip/src/fft.cpp"
#include "../mats/mats_fft.cpp"

internal void
fft_normal0(u32 count, c32 *signal, c32 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 2);

    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < count; ++index)
    {
        u32 reversedIndex = reverse_bits32(index, highBit.index);
        dest[index] = signal[reversedIndex];
    }

    u32 halfM = 1;
    u32 m = 2;
    while (m <= count)
    {
        f32 oneOverM = (2.0f * F32_PI) / (f32)m;

        for (u32 k = 0; k < count; k += m)
        {
            c32 *src0 = dest + k;
            c32 *src1 = dest + k + halfM;
            for (u32 j = 0; j < halfM; ++j)
            {
                // NOTE(michiel): w = e^(-i*2pi*j/m)
                SinCos32 cs = sincos32(-(f32)j * oneOverM);
                c32 w = complex32(cs.cos, cs.sin);
                c32 E = *src0;
                c32 O = w * *src1;
                *src0++ = E + O;
                *src1++ = E - O;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_inplace0(u32 count, c32 *signal)
{
    i_expect(is_pow2(count));
    i_expect(count > 2);

    BitScanResult highBit = find_most_significant_set_bit(count);
    for (u32 index = 0; index < count; ++index)
    {
        //u32 reversedIndex = bit_reverse(index, highBit.index);
        //reversedIndex >>= (31 - highBit.index);
        u32 reversedIndex = reverse_bits32(index, highBit.index);
        if (reversedIndex > index)
        {
            c32 temp = signal[index];
            signal[index] = signal[reversedIndex];
            signal[reversedIndex] = temp;
        }
    }

    u32 halfM = 1;
    u32 m = 2;
    while (m <= count)
    {
        f32 oneOverM = (2.0f * F32_PI) / (f32)m;

        for (u32 k = 0; k < count; k += m)
        {
            c32 *src0 = signal + k;
            c32 *src1 = signal + k + halfM;
            for (u32 j = 0; j < halfM; ++j)
            {
                // NOTE(michiel): w = e^(-i*2pi*j/m)
                SinCos32 cs = sincos32(-(f32)j * oneOverM);
                c32 w = complex32(cs.cos, cs.sin);
                c32 E = *src0;
                c32 O = w * *src1;
                *src0++ = E + O;
                *src1++ = E - O;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_normal0_64(u32 count, c64 *signal, c64 *dest)
{
    i_expect(is_pow2(count));
    i_expect(count > 2);

    BitScanResult highBit = find_most_significant_set_bit(count);
    i_expect(highBit.found);
    for (u32 index = 0; index < count; ++index)
    {
        u32 reversedIndex = reverse_bits32(index, highBit.index);
        dest[index] = signal[reversedIndex];
    }

    u32 halfM = 1;
    u32 m = 2;
    while (m <= count)
    {
        f64 oneOverM = (2.0 * F64_PI) / (f64)m;

        for (u32 k = 0; k < count; k += m)
        {
            c64 *src0 = dest + k;
            c64 *src1 = dest + k + halfM;
            for (u32 j = 0; j < halfM; ++j)
            {
                // NOTE(michiel): w = e^(-i*2pi*j/m)
                SinCos64 cs = sincos64(-(f64)j * oneOverM);
                c64 w = complex64(cs.cos, cs.sin);
                c64 E = *src0;
                c64 O = w * *src1;
                *src0++ = E + O;
                *src1++ = E - O;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

internal void
fft_inplace0_64(u32 count, c64 *signal)
{
    i_expect(is_pow2(count));
    i_expect(count > 2);

    BitScanResult highBit = find_most_significant_set_bit(count);
    for (u32 index = 0; index < count; ++index)
    {
        //u32 reversedIndex = bit_reverse(index, highBit.index);
        //reversedIndex >>= (31 - highBit.index);
        u32 reversedIndex = reverse_bits32(index, highBit.index);
        if (reversedIndex > index)
        {
            c64 temp = signal[index];
            signal[index] = signal[reversedIndex];
            signal[reversedIndex] = temp;
        }
    }

    u32 halfM = 1;
    u32 m = 2;
    while (m <= count)
    {
        f64 oneOverM = (2.0 * F64_PI) / (f64)m;

        for (u32 k = 0; k < count; k += m)
        {
            c64 *src0 = signal + k;
            c64 *src1 = signal + k + halfM;
            for (u32 j = 0; j < halfM; ++j)
            {
                // NOTE(michiel): w = e^(-i*2pi*j/m)
                SinCos64 cs = sincos64(-(f64)j * oneOverM);
                c64 w = complex64(cs.cos, cs.sin);
                c64 E = *src0;
                c64 O = w * *src1;
                *src0++ = E + O;
                *src1++ = E - O;
            }
        }
        halfM = m;
        m <<= 1;
    }
}

#include "../src/fft_evolution.cpp"

internal struct timespec
linux_get_wall_clock()
{
    struct timespec clock;
    clock_gettime(CLOCK_MONOTONIC, &clock);
    return clock;
}

internal f32
linux_get_seconds_elapsed(struct timespec start, struct timespec end)
{
    return ((f32)(end.tv_sec - start.tv_sec)
            + ((f32)(end.tv_nsec - start.tv_nsec) * 1e-9f));
}

s32 main(s32 argc, char **argv)
{
#define FFT_COUNT 64
    u32 fftCount = FFT_COUNT;
    c32 inputSignal[FFT_COUNT] = {};
    c32 outputSignal[FFT_COUNT];
    c32 outpu2Signal[FFT_COUNT];
    c32 outpu3Signal[FFT_COUNT];
    for (u32 index = 0; index < FFT_COUNT; ++index)
    {
        //inputSignal[index].real = 0.5f*cos32(5.0f * 2.0f * F32_PI * (f32)index / (f32)FFT_COUNT);
        //inputSignal[index].imag = 0.5f*sin32(5.0f * 2.0f * F32_PI * (f32)index / (f32)FFT_COUNT);
        inputSignal[index].real = sin32(5.0f * 2.0f * F32_PI * (f32)index / (f32)FFT_COUNT);
    }
    fft_normal0(fftCount, inputSignal, outputSignal);
    //fft_normal2(fftCount, inputSignal, outpu2Signal);
    fft_normal(fftCount, inputSignal, outpu2Signal);
    fft_fast(fftCount, inputSignal, outpu3Signal);
    //copy(FFT_COUNT * sizeof(c32), inputSignal, outpu3Signal);
    //fft(fftCount, outpu3Signal);
    fft_inplace0(fftCount, inputSignal);
    for (u32 index = 0; index < FFT_COUNT; ++index)
    {
        fprintf(stdout, "FFT[%3u]: ", index);

        {
            f32 magnitude = absolute32(inputSignal[index]) / (0.5f * (f32)fftCount);
            f32 phase     = argument32(inputSignal[index]);
            f32 magIndB   = 20.0f * log10_32(magnitude);
            fprintf(stdout, "%-15a % 7.3f [% 7.2f dB]", magnitude, phase, magIndB);
        }
        fprintf(stdout, " | ");
        {
            f32 magnitude = absolute32(outputSignal[index]) / (0.5f * (f32)fftCount);
            f32 phase     = argument32(outputSignal[index]);
            f32 magIndB   = 20.0f * log10_32(magnitude);
            fprintf(stdout, "%-15a % 7.3f [% 7.2f dB]", magnitude, phase, magIndB);
        }
        fprintf(stdout, " | ");
        {
            f32 magnitude = absolute32(outpu2Signal[index]) / (0.5f * (f32)fftCount);
            f32 phase     = argument32(outpu2Signal[index]);
            f32 magIndB   = 20.0f * log10_32(magnitude);
            fprintf(stdout, "%-15a % 7.3f [% 7.2f dB]", magnitude, phase, magIndB);
        }
        fprintf(stdout, " | ");
        {
            f32 magnitude = absolute32(outpu3Signal[index]) / (0.5f * (f32)fftCount);
            f32 phase     = argument32(outpu3Signal[index]);
            f32 magIndB   = 20.0f * log10_32(magnitude);
            fprintf(stdout, "%-15a % 7.3f [% 7.2f dB]", magnitude, phase, magIndB);
        }

        fprintf(stdout, "\n");
    }

    //
    // NOTE(michiel): FFT -> IFFT
    //

    for (u32 index = 0; index < FFT_COUNT; ++index)
    {
        inputSignal[index].real = sin32(5.0f * 2.0f * F32_PI * (f32)index / (f32)FFT_COUNT);
    }
    c32 reverseSignal[FFT_COUNT];
    fft_normal(FFT_COUNT, inputSignal, outputSignal);
    ifft_normal(FFT_COUNT, outputSignal, reverseSignal);

#if 0
    for (u32 index = 0; index < FFT_COUNT; ++index)
    {
        fprintf(stdout, "%3d: %7.5f -> %7.5f\n", index, inputSignal[index].real, reverseSignal[index].real);
    }
#endif

    fprintf(stdout, "%f,%f -> %f,%f; %f,%f -> %f,%f; %f,%f -> %f,%f\n",
            inputSignal[4].real, inputSignal[4].imag,
            reverseSignal[4].real, reverseSignal[4].imag,
            inputSignal[15].real, inputSignal[15].imag,
            reverseSignal[15].real, reverseSignal[15].imag,
            inputSignal[34].real, inputSignal[34].imag,
            reverseSignal[34].real, reverseSignal[34].imag);

    //
    // NOTE(michiel): Large arrays
    //

    u32 largeCount = megabytes(16);
    c32 *inputL = (c32 *)malloc(sizeof(c32) * largeCount);
    c32 *outputL = (c32 *)malloc(sizeof(c32) * largeCount);
    for (u32 index = 0; index < largeCount; ++index)
    {
        inputL[index].real = 0.5f*cos32(5.0f * 2.0f * F32_PI * (f32)index / (f32)largeCount);
        inputL[index].imag = 0.5f*sin32(5.0f * 2.0f * F32_PI * (f32)index / (f32)largeCount);;
    }

    struct timespec start = linux_get_wall_clock();
    fft_normal0(largeCount, inputL, outputL);
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    start = linux_get_wall_clock();
    fft_normal2(largeCount, inputL, outputL);
    f32 seconds2 = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    start = linux_get_wall_clock();
    fft_normal6(largeCount, inputL, outputL);
    f32 seconds3 = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    start = linux_get_wall_clock();
    fft_inexact6(largeCount, inputL, outputL);
    f32 seconds4 = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    start = linux_get_wall_clock();
    fft_normal(largeCount, inputL, outputL);
    f32 seconds5 = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    start = linux_get_wall_clock();
    fft_fast(largeCount, inputL, outputL);
    f32 seconds6 = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    start = linux_get_wall_clock();
    copy(largeCount * sizeof(c32), inputL, outputL);
    fft(largeCount, outputL);
    f32 secondsRef = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    fprintf(stdout, "Speed: (S = samples)\n");
    fprintf(stdout, "%u kS in %6.3f seconds, (1: %6.3fMS/sec) (2: %6.3fMS/sec) (6a: %6.3fMS/sec) (6b: %6.3fMS/sec) (mfft: %6.3fMS/sec) (mffft: %6.3fMS/sec) [target: %6.3fMS/sec]\n", largeCount / 1024, seconds,
            (f64)(largeCount / 1024 / 1024) / seconds,
            (f64)(largeCount / 1024 / 1024) / seconds2,
            (f64)(largeCount / 1024 / 1024) / seconds3,
            (f64)(largeCount / 1024 / 1024) / seconds4,
            (f64)(largeCount / 1024 / 1024) / seconds5,
            (f64)(largeCount / 1024 / 1024) / seconds6,
            (f64)(largeCount / 1024 / 1024) / secondsRef);

    c64 inputSignal64[FFT_COUNT] = {};
    c64 outputSignal64[FFT_COUNT];
    c64 outpu2Signal64[FFT_COUNT];
    c64 outpu3Signal64[FFT_COUNT];
    for (u32 index = 0; index < FFT_COUNT; ++index)
    {
        //inputSignal[index].real = 0.5*cos64(5.0 * 2.0 * F64_PI * (f64)index / (f64)FFT_COUNT);
        //inputSignal[index].imag = 0.5*sin64(5.0 * 2.0 * F64_PI * (f64)index / (f64)FFT_COUNT);
        inputSignal64[index].real = sin64(2.0 * 2.0 * F64_PI * (f64)index / (f64)FFT_COUNT);
    }
    fft_normal0_64(fftCount, inputSignal64, outputSignal64);
    fft_normal2_64(fftCount, inputSignal64, outpu2Signal64);
    //fft_normal5_64(fftCount, inputSignal64, outpu3Signal64);
    //fft_normal(fftCount, inputSignal, outpu2Signal);
    //fft_fast(fftCount, inputSignal, outpu3Signal);
    //copy(FFT_COUNT * sizeof(c32), inputSignal, outpu3Signal);
    //fft(fftCount, outpu3Signal);
    //fft_inplace0_64(fftCount, inputSignal64);
    copy(FFT_COUNT * sizeof(c64), inputSignal64, outpu3Signal64);
    fft_inplace6_64(fftCount, outpu3Signal64);

    f64 magnitude64[FFT_COUNT];
    f64 phase64[FFT_COUNT];
    db_from_fft64(fftCount, outpu3Signal64, magnitude64);
    unwrapped_phase_from_fft64(fftCount, outpu3Signal64, phase64);

    for (u32 index = 0; index < FFT_COUNT; ++index)
    {
        fprintf(stdout, "FFT64[%3u]: ", index);

        {
            f64 magnitude = absolute64(inputSignal64[index]) / (0.5 * (f64)fftCount);
            f64 phase     = argument64(inputSignal64[index]);
            f64 magIndB   = 20.0 * log10_64(magnitude);
            fprintf(stdout, "%-22a % 7.3f [% 7.2f dB]", magnitude, phase, magIndB);
        }
        fprintf(stdout, " | ");
        {
            f64 magnitude = absolute64(outputSignal64[index]) / (0.5 * (f64)fftCount);
            f64 phase     = argument64(outputSignal64[index]);
            f64 magIndB   = 20.0 * log10_64(magnitude);
            fprintf(stdout, "%-22a % 7.3f [% 7.2f dB]", magnitude, phase, magIndB);
        }
        fprintf(stdout, " | ");
        {
            f64 magnitude = absolute64(outpu2Signal64[index]) / (0.5 * (f64)fftCount);
            f64 phase     = argument64(outpu2Signal64[index]);
            f64 magIndB   = 20.0 * log10_64(magnitude);
            fprintf(stdout, "%-22a % 7.3f [% 7.2f dB]", magnitude, phase, magIndB);
        }
        fprintf(stdout, " | ");
        {
            fprintf(stdout, "%-22a % 7.3f [% 7.2f dB]", pow64(10.0, magnitude64[index] / 20.0), phase64[index], magnitude64[index]);
        }
#if 0
        {
            f64 magnitude = absolute64(outpu3Signal64[index]) / (0.5 * (f64)fftCount);
            f64 phase     = argument64(outpu3Signal64[index]);
            f64 magIndB   = 20.0 * log10_64(magnitude);
            fprintf(stdout, "%-22a % 7.3f [% 7.2f dB]", magnitude, phase, magIndB);
        }
#endif

        fprintf(stdout, "\n");
    }

#if 0
    //
    // NOTE(michiel): FFT -> IFFT
    //

    for (u32 index = 0; index < FFT_COUNT; ++index)
    {
        inputSignal[index].real = sin32(5.0f * 2.0f * F32_PI * (f32)index / (f32)FFT_COUNT);
    }
    c32 reverseSignal[FFT_COUNT];
    fft_normal(FFT_COUNT, inputSignal, outputSignal);
    ifft_normal(FFT_COUNT, outputSignal, reverseSignal);

#if 0
    for (u32 index = 0; index < FFT_COUNT; ++index)
    {
        fprintf(stdout, "%3d: %7.5f -> %7.5f\n", index, inputSignal[index].real, reverseSignal[index].real);
    }
#endif

    fprintf(stdout, "%f,%f -> %f,%f; %f,%f -> %f,%f; %f,%f -> %f,%f\n",
            inputSignal[4].real, inputSignal[4].imag,
            reverseSignal[4].real, reverseSignal[4].imag,
            inputSignal[15].real, inputSignal[15].imag,
            reverseSignal[15].real, reverseSignal[15].imag,
            inputSignal[34].real, inputSignal[34].imag,
            reverseSignal[34].real, reverseSignal[34].imag);
#endif

    //
    // NOTE(michiel): Large arrays
    //

    c64 *inputL64 = (c64 *)malloc(sizeof(c64) * largeCount);
    c64 *outputL64 = (c64 *)malloc(sizeof(c64) * largeCount);
    for (u32 index = 0; index < largeCount; ++index)
    {
        inputL64[index].real = 0.5*cos64(5.0 * 2.0 * F64_PI * (f64)index / (f64)largeCount);
        inputL64[index].imag = 0.5*sin64(5.0 * 2.0 * F64_PI * (f64)index / (f64)largeCount);;
    }

    start = linux_get_wall_clock();
    fft_normal0_64(largeCount, inputL64, outputL64);
    seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    start = linux_get_wall_clock();
    fft_normal2_64(largeCount, inputL64, outputL64);
    seconds2 = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    start = linux_get_wall_clock();
    fft_normal5_64(largeCount, inputL64, outputL64);
    seconds3 = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    start = linux_get_wall_clock();
    copy(largeCount * sizeof(c64), inputL64, outputL64);
    fft_inplace6_64(largeCount, outputL64);
    seconds4 = linux_get_seconds_elapsed(start, linux_get_wall_clock());

#if 0
    start = linux_get_wall_clock();
    fft_inexact6(largeCount, inputL, outputL);
    f32 seconds4 = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    start = linux_get_wall_clock();
    fft_normal(largeCount, inputL, outputL);
    f32 seconds5 = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    start = linux_get_wall_clock();
    fft_fast(largeCount, inputL, outputL);
    f32 seconds6 = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    start = linux_get_wall_clock();
    copy(largeCount * sizeof(c32), inputL, outputL);
    fft(largeCount, outputL);
    f32 secondsRef = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    fprintf(stdout, "Speed: (S = samples)\n");
    fprintf(stdout, "%u kS in %6.3f seconds, (1: %6.3fMS/sec) (2: %6.3fMS/sec) (6a: %6.3fMS/sec) (6b: %6.3fMS/sec) (mfft: %6.3fMS/sec) (mffft: %6.3fMS/sec) [target: %6.3fMS/sec]\n", largeCount / 1024, seconds,
            (f64)(largeCount / 1024 / 1024) / seconds,
            (f64)(largeCount / 1024 / 1024) / seconds2,
            (f64)(largeCount / 1024 / 1024) / seconds3,
            (f64)(largeCount / 1024 / 1024) / seconds4,
            (f64)(largeCount / 1024 / 1024) / seconds5,
            (f64)(largeCount / 1024 / 1024) / seconds6,
            (f64)(largeCount / 1024 / 1024) / secondsRef);
#else
    fprintf(stdout, "Speed: (S = samples)\n");
    fprintf(stdout, "%u kS in %6.3f seconds, (1: %6.3fMS/sec) (2: %6.3fMS/sec) (5: %6.3fMS/sec) (6: %6.3fMS/sec)\n", largeCount / 1024, seconds,
            (f64)(largeCount / 1024 / 1024) / seconds,
            (f64)(largeCount / 1024 / 1024) / seconds2,
            (f64)(largeCount / 1024 / 1024) / seconds3,
            (f64)(largeCount / 1024 / 1024) / seconds4);
#endif

    return 0;
}
