
// NOTE(michiel): The fast version are slightly less exact, the difference is in the calculation of e^(-i*2pi*j/m).
//   The fast version do an adding of the angle instead of calling sincos multiple times, this makes it run 1.6 times as fast.
internal void fft_inplace(u32 count, c32 *data);
internal void fft_inplace_fast(u32 count, c32 *data);

internal void
fft_normal(u32 count, c32 *signal, c32 *dest)
{
    copy(count * sizeof(c32), signal, dest);
    fft_inplace(count, dest);
}

internal void
fft_fast(u32 count, c32 *signal, c32 *dest)
{
    copy(count * sizeof(c32), signal, dest);
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
