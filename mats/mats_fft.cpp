
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
