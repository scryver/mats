
internal void
fft_inplace(u32 count, c32 *dest)
{
    MATS_ASSERT(is_pow2(count));
    MATS_ASSERT(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    MATS_ASSERT(highBit.found);
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
fft_inplace_fast(u32 count, c32 *dest)
{
    MATS_ASSERT(is_pow2(count));
    MATS_ASSERT(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    MATS_ASSERT(highBit.found);
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
ifft_inplace(u32 count, c32 *dest)
{
    MATS_ASSERT(is_pow2(count));
    MATS_ASSERT(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    MATS_ASSERT(highBit.found);
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

    // NOTE(michiel): w = e^(i*2pi*j/m)
    SinCos32_4x csBase = sincos32_4x(F32_4x(0.0f, 0.25f * F32_PI, 0.5f * F32_PI, 0.75f * F32_PI));
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
        f32_4x angleBase = oneOverM * angleMod0;
        f32_4x angleStep = F32_4x(4.0f) * oneOverM;

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

    f32 scale = 1.0f / (f32)count;
    for (u32 index = 0; index < count; ++index)
    {
        dest[index] *= scale;
    }
}

internal void
ifft_inplace_fast(u32 count, c32 *dest)
{
    MATS_ASSERT(is_pow2(count));
    MATS_ASSERT(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    MATS_ASSERT(highBit.found);
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
    SinCos32_4x csBase = sincos32_4x(F32_4x(0.0f, 0.25f * F32_PI, 0.5f * F32_PI, 0.75f * F32_PI));
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
    SinCos32_4x sinCosPre = sincos32_4x(F32_4x(oneOverMpre, oneOverMpre*2.0f, oneOverMpre*3.0f, oneOverMpre*4.0f));

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

        SinCos32_4x sinCos = sincos32_4x(F32_4x(oneOverM, oneOverM*3.0f, oneOverM*5.0f, oneOverM*7.0f));
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

    f32 scale = 1.0f / (f32)count;
    for (u32 index = 0; index < count; ++index)
    {
        dest[index] *= scale;
    }
}

//
// NOTE(michiel): 64-bit version
//

internal void
fft_inplace64(u32 count, c64 *dest)
{
    MATS_ASSERT(is_pow2(count));
    MATS_ASSERT(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    MATS_ASSERT(highBit.found);
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
ifft_inplace64(u32 count, c64 *dest)
{
    MATS_ASSERT(is_pow2(count));
    MATS_ASSERT(count > 8);

    u32 halfCount = count / 2;
    BitScanResult highBit = find_most_significant_set_bit(count);
    MATS_ASSERT(highBit.found);
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

    // NOTE(michiel): w = e^(i*2pi*j/m)
    SinCos64 csBase0 = sincos64(0.0);
    SinCos64 csBase1 = sincos64(0.25 * F64_PI);
    SinCos64 csBase2 = sincos64(0.5 * F64_PI);
    SinCos64 csBase3 = sincos64(0.75 * F64_PI);
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
                SinCos64 cs0 = sincos64((f64)(j + 0) * oneOverM);
                SinCos64 cs1 = sincos64((f64)(j + 1) * oneOverM);
                SinCos64 cs2 = sincos64((f64)(j + 2) * oneOverM);
                SinCos64 cs3 = sincos64((f64)(j + 3) * oneOverM);
                SinCos64 cs4 = sincos64((f64)(j + 4) * oneOverM);
                SinCos64 cs5 = sincos64((f64)(j + 5) * oneOverM);
                SinCos64 cs6 = sincos64((f64)(j + 6) * oneOverM);
                SinCos64 cs7 = sincos64((f64)(j + 7) * oneOverM);

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

    f64 scale = 1.0 / (f64)count;
    for (u32 index = 0; index < count; ++index)
    {
        dest[index] *= scale;
    }
}
