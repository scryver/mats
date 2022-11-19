
// NOTE(michiel): The fast version are slightly less exact, the difference is in the calculation of e^(-i*2pi*j/m).
//   The fast version do an adding of the angle instead of calling sincos multiple times, this makes it run 1.8 to 3 times as fast.

internal void fft_inplace(u32 count, c32 *data);
internal void fft_inplace_fast(u32 count, c32 *data);

internal void
fft_normal(u32 count, c32 *signal, c32 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = signal[index];
    }
    fft_inplace(count, dest);
}

internal void
fft_fast(u32 count, c32 *signal, c32 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = signal[index];
    }
    fft_inplace_fast(count, dest);
}

internal void
fft_real(u32 count, f32 *signal, c32 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = complex32(signal[index], 0.0f);
    }
    fft_inplace(count, dest);
}

internal void
fft_real_fast(u32 count, f32 *signal, c32 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = complex32(signal[index], 0.0f);
    }
    fft_inplace_fast(count, dest);
}

internal void ifft_inplace(u32 count, c32 *data);
internal void ifft_inplace_fast(u32 count, c32 *data);

internal void
ifft_normal(u32 count, c32 *signal, c32 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = signal[index];
    }
    ifft_inplace(count, dest);
}

internal void
ifft_fast(u32 count, c32 *signal, c32 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = signal[index];
    }
    ifft_inplace_fast(count, dest);
}

internal void
ifft_real(u32 count, f32 *signal, c32 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = complex32(signal[index], 0.0f);
    }
    ifft_inplace(count, dest);
}

internal void
ifft_real_fast(u32 count, f32 *signal, c32 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = complex32(signal[index], 0.0f);
    }
    ifft_inplace_fast(count, dest);
}

//
// NOTE(michiel): 64-bit version
//

internal void fft_inplace64(u32 count, c64 *data);
internal void fft_inplace_fast64(u32 count, c64 *data);

internal void
fft_normal64(u32 count, c64 *signal, c64 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = signal[index];
    }
    fft_inplace64(count, dest);
}

internal void
fft_fast64(u32 count, c64 *signal, c64 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = signal[index];
    }
    fft_inplace_fast64(count, dest);
}

internal void
fft_real64(u32 count, f64 *signal, c64 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = complex64(signal[index], 0.0);
    }
    fft_inplace64(count, dest);
}

internal void
fft_real_fast64(u32 count, f64 *signal, c64 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = complex64(signal[index], 0.0);
    }
    fft_inplace_fast64(count, dest);
}

internal void ifft_inplace64(u32 count, c64 *data);
internal void ifft_inplace_fast64(u32 count, c64 *data);

internal void
ifft_normal64(u32 count, c64 *signal, c64 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = signal[index];
    }
    ifft_inplace64(count, dest);
}

internal void
ifft_fast64(u32 count, c64 *signal, c64 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = signal[index];
    }
    ifft_inplace_fast64(count, dest);
}

internal void
ifft_real64(u32 count, f64 *signal, c64 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = complex64(signal[index], 0.0);
    }
    ifft_inplace64(count, dest);
}

internal void
ifft_real_fast64(u32 count, f64 *signal, c64 *dest)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = complex64(signal[index], 0.0);
    }
    ifft_inplace_fast64(count, dest);
}
