
internal void fft_shift32(u32 count, c32 *data);
internal void fft_shift64(u32 count, c64 *data);

internal void magnitude_from_fft32(u32 count, c32 *source, f32 *dest);
internal void magnitude_from_fft64(u32 count, c64 *source, f64 *dest);
internal void db_from_fft32(u32 count, c32 *source, f32 *dest);
internal void db_from_fft64(u32 count, c64 *source, f64 *dest);
internal void phase_from_fft32(u32 count, c32 *source, f32 *dest);
internal void phase_from_fft64(u32 count, c64 *source, f64 *dest);
internal void unwrapped_phase_from_fft32(u32 count, c32 *source, f32 *dest);
internal void unwrapped_phase_from_fft64(u32 count, c64 *source, f64 *dest);

internal u32
fft_table_count(u32 count)
{
    // NOTE(michiel): Amount of f32 or f64 to reserve for the sincos table
    u32 result = 2 * (count * 2 - 8);
    return result;
}

// NOTE(michiel): The fast version are slightly less exact, the difference is in the calculation of e^(-i*2pi*j/m).
//   The fast version do an adding of the angle instead of calling sincos multiple times, this makes it run 1.8 to 3 times as fast.

internal void fft_inplace(u32 count, c32 *data);
internal void fft_inplace_fast(u32 count, c32 *data);

internal void fft_build_table_f32(u32 count, f32 *dest);
internal void fft_inplace_table(u32 count, c32 *dest, f32 *sinCosTable);

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
fft_table(u32 count, c32 *signal, c32 *dest, f32 *sinCosTable)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = signal[index];
    }
    fft_inplace_table(count, dest, sinCosTable);
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

internal void
fft_real_table(u32 count, f32 *signal, c32 *dest, f32 *sinCosTable)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = complex32(signal[index], 0.0f);
    }
    fft_inplace_table(count, dest, sinCosTable);
}

internal void ifft_inplace(u32 count, c32 *data);
internal void ifft_inplace_fast(u32 count, c32 *data);

internal void ifft_build_table_f32(u32 count, f32 *dest);
internal void ifft_inplace_table(u32 count, c32 *dest, f32 *sinCosTable);

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
ifft_table(u32 count, c32 *signal, c32 *dest, f32 *sinCosTable)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = signal[index];
    }
    ifft_inplace_table(count, dest, sinCosTable);
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

internal void
ifft_real_table(u32 count, f32 *signal, c32 *dest, f32 *sinCosTable)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = complex32(signal[index], 0.0f);
    }
    ifft_inplace_table(count, dest, sinCosTable);
}

//
// NOTE(michiel): 64-bit version
//

internal void fft_inplace64(u32 count, c64 *data);
internal void fft_inplace_fast64(u32 count, c64 *data);

internal void fft_build_table_f64(u32 count, f64 *dest);
internal void fft_inplace_table64(u32 count, c64 *dest, f64 *sinCosTable);

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
fft_table64(u32 count, c64 *signal, c64 *dest, f64 *sinCosTable)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = signal[index];
    }
    fft_inplace_table64(count, dest, sinCosTable);
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

internal void
fft_real_table64(u32 count, f64 *signal, c64 *dest, f64 *sinCosTable)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = complex64(signal[index], 0.0);
    }
    fft_inplace_table64(count, dest, sinCosTable);
}

internal void ifft_inplace64(u32 count, c64 *data);
internal void ifft_inplace_fast64(u32 count, c64 *data);

internal void ifft_build_table_f64(u32 count, f64 *dest);
internal void ifft_inplace_table64(u32 count, c64 *dest, f64 *sinCosTable);

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
ifft_table64(u32 count, c64 *signal, c64 *dest, f64 *sinCosTable)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = signal[index];
    }
    ifft_inplace_table64(count, dest, sinCosTable);
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

internal void
ifft_real_table64(u32 count, f64 *signal, c64 *dest, f64 *sinCosTable)
{
    for (u32 index = 0; index < count; ++index) {
        dest[index] = complex64(signal[index], 0.0);
    }
    ifft_inplace_table64(count, dest, sinCosTable);
}
