
global const f32 gZeroF32          =  0.0f;
global const f32 gOneF32           =  1.0000000000e+00f;  /* 0x3f800000 */
global const f32 gHalfF32          =  5.0000000000e-01f;  /* 0x3f000000 */
global const f32 g2pow8F32         =  2.5600000000e+02f;  /* 0x43800000 = 2^8 */
global const f32 g2pow25F32        =  3.355443200e+07f;   /* 0x4c000000 = 2^25 */
global const f32 g2powMin8F32      =  3.9062500000e-03f;  /* 0x3b800000 = 2^-8 */
global const f32 g2powMin25F32     =  2.9802322388e-08f;  /* 0x33000000 = 2^-25 */
global const f32 g2powMin100F32    =  7.8886090522e-31f;  /* 0x0d800000 = 2^-100 */
global const f32 gLn2F32           =  6.9314718246e-01f;  /* 0x3f317218 = ln(2) */
global const f32 gInvLn2F32        =  1.4426950216e+00f;  /* 0x3fb8aa3b = 1/ln(2) */
global const f32 gInvLn10F32       =  4.3429449201e-01f;  /* 0x3ede5bd9 = 1/ln(10) */
global const f32 gHugeF32          =  1.0e+30f;
global const f32 gTinyF32          =  1.0e-30f;

global const f64 gZeroF64          =  0.0;
global const f64 gOneF64           =  1.0;
global const f64 gHalfF64          =  0.5;
//global const f64 g2pow8F64         =  2.5600000000e+02;  /* 0x43800000 = 2^8 */
global const f64 g2pow24F64        =  1.67772160000000000000e+07;  /* 0x4170000000000000 = 2^24 */
global const f64 g2pow54F64        =  1.80143985094819840000e+16;  /* 0x4350000000000000 = 2^54 */
//global const f64 g2powMin8F64      =  3.9062500000e-03;  /* 0x3b800000 = 2^-8 */
global const f64 g2powMin24F64     =  5.96046447753906250000e-08;  /* 0x3E70000000000000 = 2^-24 */
global const f64 g2powMin54F64     =  5.55111512312578270212e-17;  /* 0x3C90000000000000 = 2^-54 */
global const f64 g2powMin1000F64   =  9.33263618503218878990e-302; /* 0x0170000000000000 = 2^-1000 */
global const f64 gLn2F64           =  6.93147180559945286227e-01;  /* 0x3FE62E42FEFA39EF = ln(2) */
global const f64 gInvLn2F64        =  1.44269504088896338700e+00;  /* 0x3ff71547652b82fe = 1/ln(2) */
global const f64 gInvLn10F64       =  4.34294481903251816668e-01;  /* 0x3FDBCB7B1526E50E = 1/ln(10) */
global const f64 gHugeF64          =  1.0e+300;
global const f64 gTinyF64          =  1.0e-300;

global const f32 gPiF32            =  3.1415927410e+00f; /* 0x40490fdb, 24 bits of pi */
global const f32 gPiF32_lo         = -8.7422776573e-08f; /* 0x40490fdb, pi - gPiF32 */
global const f32 gPiOver2F32       =  1.5707963705e+00f; /* 0x3fc90fdb, 24 bits of pi / 2 */
global const f32 gPiOver4F32       =  7.8539818525e-01f; /* 0x3f490fdb, 24 bits of pi / 4 */
global const f32 g2OverPiF32       =  6.3661980629e-01f; /* 0x3f22f984, 24 bits of 2 / pi */
global const f32 gPiOver4F32_hi    =  7.8539812565e-01f; /* 0x3F490FDA, 24 bits of pi / 4 */
global const f32 gPiOver4F32_lo    =  3.7748947079e-08f; /* 0x33222168, pi / 4 - gPiOver4F32_hi */
global const f32 gPiOver2F32_1     =  1.5707855225e+00f; /* 0x3fc90f80, first 17 bit of pi/2 */
global const f32 gPiOver2F32_1t    =  1.0804334124e-05f; /* 0x37354443, pi/2 - pio2_1 */
global const f32 gPiOver2F32_2     =  1.0804273188e-05f; /* 0x37354400, second 17 bit of pi/2 */
global const f32 gPiOver2F32_2t    =  6.0770999344e-11f; /* 0x2e85a308, pi/2 - (pio2_1+pio2_2) */
global const f32 gPiOver2F32_3     =  6.0770943833e-11f; /* 0x2e85a300, third  17 bit of pi/2 */
global const f32 gPiOver2F32_3t    =  6.1232342629e-17f; /* 0x248d3132, pi/2 - (pio2_1+pio2_2+pio2_3) */
global const f32 gLog10F32_2_hi    =  3.0102920532e-01f; /* 0x3e9a2080, log10(2) */
global const f32 gLog10F32_2_lo    =  7.9034151668e-07f; /* 0x355427db, log10(2) - gLog10_2_hi */

global const f64 gPiF64            =  3.14159265358979311600e+00; /* 0x400921FB54442D18, 52 bits of pi */
global const f64 gPiF64_lo         =  1.22464679914735317720e-16; /* 0x3CA1A62633145C07, pi - gPiF64 */
global const f64 gPiOver2F64       =  1.57079632679489655800e+00; /* 0x3FF921FB54442D18, 52 bits of pi / 2 */
global const f64 gPiOver2F64_lo    =  6.12323399573676603587e-17; /* 0x3C91A62633145C07, pi / 2 - gPiOver2F64 */
global const f64 gPiOver4F64       =  7.85398163397448278999e-01; /* 0x3FE921FB54442D18, 52 bits of pi / 4 */
global const f64 gPiOver4F64_lo    =  3.06161699786838301793e-17; /* 0x3C81A62633145C07, pi / 4 - gPiOver4F64 */
global const f64 g2OverPiF64       =  6.36619772367581382433e-01; /* 0x3FE45F306DC9C883, 52 bits of 2 / pi */
global const f64 gPiOver2F64_1     =  1.57079632673412561417e+00; /* 0x3FF921FB54400000, first bits of pi/2 */
global const f64 gPiOver2F64_1t    =  6.07710050650619224932e-11; /* 0x3DD0B4611A626331, pi/2 - pio2_1 */
global const f64 gPiOver2F64_2     =  6.07710050630396597660e-11; /* 0x3DD0B4611A600000, second bits of pi/2 */
global const f64 gPiOver2F64_2t    =  2.02226624879595063154e-21; /* 0x3BA3198A2E037073, pi/2 - (pio2_1 + pio2_2) */
global const f64 gPiOver2F64_3     =  2.02226624871116645580e-21; /* 0x3BA3198A2E000000, third bits of pi/2 */
global const f64 gPiOver2F64_3t    =  8.47842766036889956997e-32; /* 0x397B839A252049C1, pi/2 - (pio2_1 + pio2_2 + pio2_3) */

global const f64 gLog10F64_2_hi    =  3.01029995663611771306e-01; /* 0x3FD34413509F6000, log10(2) */
global const f64 gLog10F64_2_lo    =  3.69423907715893078616e-13; /* 0x3D59FEF311F12B36, log10(2) - gLog10_2_hi */

global const f32 gHalfSignF32s[2]  = {0.5f, -0.5f};
global const f32 gLn2HighF32s[2]   = {6.9313812256e-01, -6.9313812256e-01}; /* 0x3f317180, 0xbf317180 */
global const f32 gLn2LowF32s[2]    = {9.0580006145e-06, -9.0580006145e-06}; /* 0x3717f7d1, 0xb717f7d1 */

global const f64 gHalfSignF64s[2]  = {0.5, -0.5};
global const f64 gLn2HighF64s[2]   = {6.93147180369123816490e-01, -6.93147180369123816490e-01}; // 0x3fe62e42fee00000, 0xbfe62e42fee00000
global const f64 gLn2LowF64s[2]    = {1.90821492927058770002e-10, -1.90821492927058770002e-10}; // 0x3dea39ef35793c76. 0xbdea39ef35793c76

// NOTE(michiel): Pi
global const s32 g2OverPiF64_Table[] = {
    0x00A2F983, 0x006E4E44, 0x001529FC, 0x002757D1, 0x00F534DD, 0x00C0DB62,
    0x0095993C, 0x00439041, 0x00FE5163, 0x00ABDEBB, 0x00C561B7, 0x00246E3A,
    0x00424DD2, 0x00E00649, 0x002EEA09, 0x00D1921C, 0x00FE1DEB, 0x001CB129,
    0x00A73EE8, 0x008235F5, 0x002EBB44, 0x0084E99C, 0x007026B4, 0x005F7E41,
    0x003991D6, 0x00398353, 0x0039F49C, 0x00845F8B, 0x00BDF928, 0x003B1FF8,
    0x0097FFDE, 0x0005980F, 0x00EF2F11, 0x008B5A0A, 0x006D1F6D, 0x00367ECF,
    0x0027CB09, 0x00B74F46, 0x003F669E, 0x005FEA2D, 0x007527BA, 0x00C7EBE5,
    0x00F17B3D, 0x000739F7, 0x008A5292, 0x00EA6BFB, 0x005FB11F, 0x008D5D08,
    0x00560330, 0x0046FC7B, 0x006BABF0, 0x00CFBC20, 0x009AF436, 0x001DA9E3,
    0x0091615E, 0x00E61B08, 0x00659985, 0x005F14A0, 0x0068408D, 0x00FFD880,
    0x004D7327, 0x00310606, 0x001556CA, 0x0073A8C9, 0x0060E27B, 0x00C08C6B,
};
global const s32 gNPiOver2F64_hw[] = {
    0x3FF921FB, 0x400921FB, 0x4012D97C, 0x401921FB, 0x401F6A7A, 0x4022D97C,
    0x4025FDBB, 0x402921FB, 0x402C463A, 0x402F6A7A, 0x4031475C, 0x4032D97C,
    0x40346B9C, 0x4035FDBB, 0x40378FDB, 0x403921FB, 0x403AB41B, 0x403C463A,
    0x403DD85A, 0x403F6A7A, 0x40407E4C, 0x4041475C, 0x4042106C, 0x4042D97C,
    0x4043A28C, 0x40446B9C, 0x404534AC, 0x4045FDBB, 0x4046C6CB, 0x40478FDB,
    0x404858EB, 0x404921FB,
};
global const f64 gPiOver2F64_Table[] = {
    1.57079625129699707031e+00, /* 0x3FF921FB40000000 */
    7.54978941586159635335e-08, /* 0x3E74442D00000000 */
    5.39030252995776476554e-15, /* 0x3CF8469880000000 */
    3.28200341580791294123e-22, /* 0x3B78CC5160000000 */
    1.27065575308067607349e-29, /* 0x39F01B8380000000 */
    1.22933308981111328932e-36, /* 0x387A252040000000 */
    2.73370053816464559624e-44, /* 0x36E3822280000000 */
    2.16741683877804819444e-51, /* 0x3569F31D00000000 */
};

// NOTE(michiel): Atan
global const f32 gAtanHiF32[] = {
    4.6364760399e-01, /* atan(0.5)hi 0x3eed6338 */
    7.8539812565e-01, /* atan(1.0)hi 0x3f490fda */
    9.8279368877e-01, /* atan(1.5)hi 0x3f7b985e */
    1.5707962513e+00, /* atan(inf)hi 0x3fc90fda */
};
global const f32 gAtanLoF32[] = {
    5.0121582440e-09, /* atan(0.5)lo 0x31ac3769 */
    3.7748947079e-08, /* atan(1.0)lo 0x33222168 */
    3.4473217170e-08, /* atan(1.5)lo 0x33140fb4 */
    7.5497894159e-08, /* atan(inf)lo 0x33a22168 */
};

global const f64 gAtanHiF64[] = {
    4.63647609000806093515e-01, /* atan(0.5)hi 0x3FDDAC670561BB4F */
    gPiOver4F64,                /* atan(1.0)hi 0x3FE921FB54442D18 */
    9.82793723247329054082e-01, /* atan(1.5)hi 0x3FEF730BD281F69B */
    gPiOver2F64,                /* atan(inf)hi 0x3FF921FB54442D18 */
};
global const f64 gAtanLoF64[] = {
    2.26987774529616870924e-17, /* atan(0.5)lo 0x3C7A2B7F222F65E2 */
    gPiOver4F64_lo,             /* atan(1.0)lo 0x3C81A62633145C07 */
    1.39033110312309984516e-17, /* atan(1.5)lo 0x3C7007887AF0CBBD */
    gPiOver2F64_lo,             /* atan(inf)lo 0x3C91A62633145C07 */
};


// NOTE(michiel): Exp
#define EXP2_32_TABLE_BITS 5
#define EXP2_32_POLY_ORDER 3
#define EXP2_32_N  ((f64)(1 << EXP2_32_TABLE_BITS))
global const f64 gExp2F32_Shift = 0x1.8p+52;
global const f64 gExp2F32_ShiftScaled = 0x1.8p+52 / EXP2_32_N;
global const f64 gExp2F32_InvLn2Scaled = 0x1.71547652b82fep+0 * EXP2_32_N;
global const f64 gExp2F32_Poly[EXP2_32_POLY_ORDER] = { 0x1.c6af84b912394p-5, 0x1.ebfce50fac4f3p-3, 0x1.62e42ff0c52d6p-1 };
global const f64 gExp2F32_PolyScaled[EXP2_32_POLY_ORDER] = {
    0x1.c6af84b912394p-5 / EXP2_32_N / EXP2_32_N / EXP2_32_N,
    0x1.ebfce50fac4f3p-3 / EXP2_32_N / EXP2_32_N,
    0x1.62e42ff0c52d6p-1 / EXP2_32_N,
};
global const u64 gExp2F32_Table[1 << EXP2_32_TABLE_BITS] = {
    0x3ff0000000000000ULL, 0x3fefd9b0d3158574ULL, 0x3fefb5586cf9890fULL, 0x3fef9301d0125b51ULL,
    0x3fef72b83c7d517bULL, 0x3fef54873168b9aaULL, 0x3fef387a6e756238ULL, 0x3fef1e9df51fdee1ULL,
    0x3fef06fe0a31b715ULL, 0x3feef1a7373aa9cbULL, 0x3feedea64c123422ULL, 0x3feece086061892dULL,
    0x3feebfdad5362a27ULL, 0x3feeb42b569d4f82ULL, 0x3feeab07dd485429ULL, 0x3feea47eb03a5585ULL,
    0x3feea09e667f3bcdULL, 0x3fee9f75e8ec5f74ULL, 0x3feea11473eb0187ULL, 0x3feea589994cce13ULL,
    0x3feeace5422aa0dbULL, 0x3feeb737b0cdc5e5ULL, 0x3feec49182a3f090ULL, 0x3feed503b23e255dULL,
    0x3feee89f995ad3adULL, 0x3feeff76f2fb5e47ULL, 0x3fef199bdd85529cULL, 0x3fef3720dcef9069ULL,
    0x3fef5818dcfba487ULL, 0x3fef7c97337b9b5fULL, 0x3fefa4afa2a490daULL, 0x3fefd0765b6e4540ULL,
};
#undef EXP2_32_N

// NOTE(michiel): Log
#define LOG32_TABLE_BITS 4
#define LOG32_POLY_ORDER 4
struct LogTable
{
    f64 invC;
    f64 logC;
};
global const f64 gLogF32_Ln2 = 0x1.62e42fefa39efp-1;
global const f64 gLogF32_Poly[LOG32_POLY_ORDER - 1] = {
    -0x1.00ea348b88334p-2,
    0x1.5575b0be00b6ap-2,
    -0x1.ffffef20a4123p-2,
};
global const LogTable gLogF32_Table[1 << LOG32_TABLE_BITS] = {
    { 0x1.661ec79f8f3bep+0, -0x1.57bf7808caadep-2 },
    { 0x1.571ed4aaf883dp+0, -0x1.2bef0a7c06ddbp-2 },
    { 0x1.49539f0f010bp+0, -0x1.01eae7f513a67p-2 },
    { 0x1.3c995b0b80385p+0, -0x1.b31d8a68224e9p-3 },
    { 0x1.30d190c8864a5p+0, -0x1.6574f0ac07758p-3 },
    { 0x1.25e227b0b8eap+0, -0x1.1aa2bc79c81p-3 },
    { 0x1.1bb4a4a1a343fp+0, -0x1.a4e76ce8c0e5ep-4 },
    { 0x1.12358f08ae5bap+0, -0x1.1973c5a611cccp-4 },
    { 0x1.0953f419900a7p+0, -0x1.252f438e10c1ep-5 },
    { 0x1p+0, 0x0p+0 },
    { 0x1.e608cfd9a47acp-1, 0x1.aa5aa5df25984p-5 },
    { 0x1.ca4b31f026aap-1, 0x1.c5e53aa362eb4p-4 },
    { 0x1.b2036576afce6p-1, 0x1.526e57720db08p-3 },
    { 0x1.9c2d163a1aa2dp-1, 0x1.bc2860d22477p-3 },
    { 0x1.886e6037841edp-1, 0x1.1058bc8a07ee1p-2 },
    { 0x1.767dcf5534862p-1, 0x1.4043057b6ee09p-2 },
};

// NOTE(michiel): Log2
#define LOG2_32_TABLE_BITS 4
#define LOG2_32_POLY_ORDER 4
global const f64 gLog2F32_Poly[LOG2_32_POLY_ORDER] = {
    -0x1.712b6f70a7e4dp-2,
    0x1.ecabf496832ep-2,
    -0x1.715479ffae3dep-1,
    0x1.715475f35c8b8p0,
};
global const LogTable gLog2F32_Table[1 << LOG2_32_TABLE_BITS] = {
    { 0x1.661ec79f8f3bep+0, -0x1.efec65b963019p-2 },
    { 0x1.571ed4aaf883dp+0, -0x1.b0b6832d4fca4p-2 },
    { 0x1.49539f0f010bp+0, -0x1.7418b0a1fb77bp-2 },
    { 0x1.3c995b0b80385p+0, -0x1.39de91a6dcf7bp-2 },
    { 0x1.30d190c8864a5p+0, -0x1.01d9bf3f2b631p-2 },
    { 0x1.25e227b0b8eap+0, -0x1.97c1d1b3b7afp-3 },
    { 0x1.1bb4a4a1a343fp+0, -0x1.2f9e393af3c9fp-3 },
    { 0x1.12358f08ae5bap+0, -0x1.960cbbf788d5cp-4 },
    { 0x1.0953f419900a7p+0, -0x1.a6f9db6475fcep-5 },
    { 0x1p+0, 0x0p+0 },
    { 0x1.e608cfd9a47acp-1, 0x1.338ca9f24f53dp-4 },
    { 0x1.ca4b31f026aap-1, 0x1.476a9543891bap-3 },
    { 0x1.b2036576afce6p-1, 0x1.e840b4ac4e4d2p-3 },
    { 0x1.9c2d163a1aa2dp-1, 0x1.40645f0c6651cp-2 },
    { 0x1.886e6037841edp-1, 0x1.88e9c2c1b9ff8p-2 },
    { 0x1.767dcf5534862p-1, 0x1.ce0a44eb17bccp-2 },
};

// NOTE(michiel): Pow
#define POW32_LOG2_TABLE_BITS 4
#define POW32_LOG2_POLY_ORDER 5
#define POW32_SCALE ((f64)(1 << EXP2_32_TABLE_BITS))
global const f64 gPowF32_Log2PolyScaled[POW32_LOG2_POLY_ORDER] = {
    0x1.27616c9496e0bp-2 * POW32_SCALE,
    -0x1.71969a075c67ap-2 * POW32_SCALE,
    0x1.ec70a6ca7baddp-2 * POW32_SCALE,
    -0x1.7154748bef6c8p-1 * POW32_SCALE,
    0x1.71547652ab82bp0 * POW32_SCALE,
}; // __powf_log2_data.poly
global const LogTable gPowF32_Log2TableScaled[1 << POW32_LOG2_TABLE_BITS] = {
    { 0x1.661ec79f8f3bep+0, -0x1.efec65b963019p-2 * POW32_SCALE },
    { 0x1.571ed4aaf883dp+0, -0x1.b0b6832d4fca4p-2 * POW32_SCALE },
    { 0x1.49539f0f010bp+0, -0x1.7418b0a1fb77bp-2 * POW32_SCALE },
    { 0x1.3c995b0b80385p+0, -0x1.39de91a6dcf7bp-2 * POW32_SCALE },
    { 0x1.30d190c8864a5p+0, -0x1.01d9bf3f2b631p-2 * POW32_SCALE },
    { 0x1.25e227b0b8eap+0, -0x1.97c1d1b3b7afp-3 * POW32_SCALE },
    { 0x1.1bb4a4a1a343fp+0, -0x1.2f9e393af3c9fp-3 * POW32_SCALE },
    { 0x1.12358f08ae5bap+0, -0x1.960cbbf788d5cp-4 * POW32_SCALE },
    { 0x1.0953f419900a7p+0, -0x1.a6f9db6475fcep-5 * POW32_SCALE },
    { 0x1p+0, 0x0p+0 * POW32_SCALE },
    { 0x1.e608cfd9a47acp-1, 0x1.338ca9f24f53dp-4 * POW32_SCALE },
    { 0x1.ca4b31f026aap-1, 0x1.476a9543891bap-3 * POW32_SCALE },
    { 0x1.b2036576afce6p-1, 0x1.e840b4ac4e4d2p-3 * POW32_SCALE },
    { 0x1.9c2d163a1aa2dp-1, 0x1.40645f0c6651cp-2 * POW32_SCALE },
    { 0x1.886e6037841edp-1, 0x1.88e9c2c1b9ff8p-2 * POW32_SCALE },
    { 0x1.767dcf5534862p-1, 0x1.ce0a44eb17bccp-2 * POW32_SCALE },
}; // __powf_log2_data.tab

global const f64 gPowF32_Log2Poly[POW32_LOG2_POLY_ORDER] = {
    0x1.27616c9496e0bp-2,
    -0x1.71969a075c67ap-2,
    0x1.ec70a6ca7baddp-2,
    -0x1.7154748bef6c8p-1,
    0x1.71547652ab82bp0,
};
global const LogTable gPowF32_Log2Table[1 << POW32_LOG2_TABLE_BITS] = {
    { 0x1.661ec79f8f3bep+0, -0x1.efec65b963019p-2 },
    { 0x1.571ed4aaf883dp+0, -0x1.b0b6832d4fca4p-2 },
    { 0x1.49539f0f010bp+0, -0x1.7418b0a1fb77bp-2 },
    { 0x1.3c995b0b80385p+0, -0x1.39de91a6dcf7bp-2 },
    { 0x1.30d190c8864a5p+0, -0x1.01d9bf3f2b631p-2 },
    { 0x1.25e227b0b8eap+0, -0x1.97c1d1b3b7afp-3 },
    { 0x1.1bb4a4a1a343fp+0, -0x1.2f9e393af3c9fp-3 },
    { 0x1.12358f08ae5bap+0, -0x1.960cbbf788d5cp-4 },
    { 0x1.0953f419900a7p+0, -0x1.a6f9db6475fcep-5 },
    { 0x1p+0, 0x0p+0 },
    { 0x1.e608cfd9a47acp-1, 0x1.338ca9f24f53dp-4 },
    { 0x1.ca4b31f026aap-1, 0x1.476a9543891bap-3 },
    { 0x1.b2036576afce6p-1, 0x1.e840b4ac4e4d2p-3 },
    { 0x1.9c2d163a1aa2dp-1, 0x1.40645f0c6651cp-2 },
    { 0x1.886e6037841edp-1, 0x1.88e9c2c1b9ff8p-2 },
    { 0x1.767dcf5534862p-1, 0x1.ce0a44eb17bccp-2 },
};

#define EXP64_TABLE_BITS 7
#define EXP64_POLY_ORDER 5
/* Use polynomial that is optimized for a wider input range.  This may be
   needed for good precision in non-nearest rounding and !TOINT_INTRINSICS.  */
#define EXP64_POLY_WIDE 0
/* Use close to nearest rounding toint when !TOINT_INTRINSICS.  This may be
   needed for good precision in non-nearest rouning and !EXP_POLY_WIDE.  */
#define EXP64_USE_TOINT_NARROW 0
#define EXP2_64_POLY_ORDER 5
#define EXP2_64_POLY_WIDE 0
struct ExpData
{
    f64 invln2N;
    f64 shift;
    f64 negln2hiN;
    f64 negln2loN;
    f64 poly[4]; /* Last four coefficients.  */
    f64 exp2_shift;
    f64 exp2_poly[EXP2_64_POLY_ORDER];
    u64 tab[2*(1 << EXP64_TABLE_BITS)];
};

#define LOG64_TABLE_BITS 7
#define LOG64_POLY_ORDER 6
#define LOG64_POLY1_ORDER 12
struct LogData
{
    f64 ln2hi;
    f64 ln2lo;
    f64 poly[LOG64_POLY_ORDER - 1]; /* First coefficient is 1.  */
    f64 poly1[LOG64_POLY1_ORDER - 1];
    struct {f64 invc, logc;} tab[1 << LOG64_TABLE_BITS];
#if !MATS_HAVE_FAST_FMA
    struct {f64 chi, clo;} tab2[1 << LOG64_TABLE_BITS];
#endif
};

#define LOG2_64_TABLE_BITS 6
#define LOG2_64_POLY_ORDER 7
#define LOG2_64_POLY1_ORDER 11
struct Log2Data
{
    f64 invln2hi;
    f64 invln2lo;
    f64 poly[LOG2_64_POLY_ORDER - 1];
    f64 poly1[LOG2_64_POLY1_ORDER - 1];
    struct {f64 invc, logc;} tab[1 << LOG2_64_TABLE_BITS];
#if !MATS_HAVE_FAST_FMA
    struct {f64 chi, clo;} tab2[1 << LOG2_64_TABLE_BITS];
#endif
};

#define POW64_LOG_TABLE_BITS 7
#define POW64_LOG_POLY_ORDER 8
struct PowLogData
{
    f64 ln2hi;
    f64 ln2lo;
    f64 poly[POW64_LOG_POLY_ORDER - 1]; // NOTE(michiel): First coefficient is 1.
    // NOTE(michiel): the pad field is unused, but allows slightly faster indexing.
    struct {f64 invc, pad, logc, logctail;} tab[1 << POW64_LOG_TABLE_BITS];
};

#define EXP64_N  (1 << EXP64_TABLE_BITS)
global const ExpData gExp64Data = {
    // N/ln2
    .invln2N = 0x1.71547652b82fep0 * EXP64_N,
    // Used for rounding when !TOINT_INTRINSICS
#if EXP64_USE_TOINT_NARROW
    .shift = 0x1800000000.8p0,
#else
    .shift = 0x1.8p52,
#endif

    // -ln2/N
#if EXP64_N == 64
    .negln2hiN = -0x1.62e42fefa0000p-7,
    .negln2loN = -0x1.cf79abc9e3b3ap-46,
#elif EXP64_N == 128
    .negln2hiN = -0x1.62e42fefa0000p-8,
    .negln2loN = -0x1.cf79abc9e3b3ap-47,
#elif EXP64_N == 256
    .negln2hiN = -0x1.62e42fefc0000p-9,
    .negln2loN = 0x1.c610ca86c3899p-45,
#endif
    // exp polynomial coefficients.
    .poly = {
#if EXP64_N == 64 && EXP64_POLY_ORDER == 5 && !EXP64_POLY_WIDE
        // abs error: 1.5543*2^-60
        // ulp error: 0.529 (0.533 without fma)
        // if |x| < ln2/128+eps
        // abs error if |x| < ln2/64: 1.7157*2^-50
        0x1.fffffffffdbcdp-2,
        0x1.555555555444cp-3,
        0x1.555573c6a9f7dp-5,
        0x1.1111266d28935p-7,
#elif EXP64_N == 64 && EXP64_POLY_ORDER == 6 && EXP64_POLY_WIDE
        // abs error: 1.6735*2^-64
        // ulp error: 0.518 (0.522 without fma)
        // if |x| < ln2/64
        0x1.5555555548f9ap-3,
        0x1.555555554bf5dp-5,
        0x1.11115b75f0f4dp-7,
        0x1.6c171a6b6303ep-10,
#elif EXP64_N == 128 && EXP64_POLY_ORDER == 5 && !EXP64_POLY_WIDE
        // abs error: 1.555*2^-66
        // ulp error: 0.509 (0.511 without fma)
        // if |x| < ln2/256+eps
        // abs error if |x| < ln2/256+0x1p-15: 1.09*2^-65
        // abs error if |x| < ln2/128: 1.7145*2^-56
        0x1.ffffffffffdbdp-2,
        0x1.555555555543cp-3,
        0x1.55555cf172b91p-5,
        0x1.1111167a4d017p-7,
#elif EXP64_N == 128 && EXP64_POLY_ORDER == 5 && EXP64_POLY_WIDE
        // abs error: 1.5542*2^-60
        // ulp error: 0.521 (0.523 without fma)
        // if |x| < ln2/128
        0x1.fffffffffdbcep-2,
        0x1.55555555543c2p-3,
        0x1.555573c64f2e3p-5,
        0x1.111126b4eff73p-7,
#elif EXP64_N == 128 && EXP64_POLY_ORDER == 6 && EXP64_POLY_WIDE
        // abs error: 1.6861*2^-71
        // ulp error: 0.509 (0.511 without fma)
        // if |x| < ln2/128
        0x1.55555555548fdp-3,
        0x1.555555555658fp-5,
        0x1.111123a859bb6p-7,
        0x1.6c16ba6920cabp-10,
#elif EXP64_N == 256 && EXP64_POLY_ORDER == 4 && !EXP64_POLY_WIDE
        // abs error: 1.43*2^-58
        // ulp error: 0.549 (0.550 without fma)
        // if |x| < ln2/512
        0x1p0, // unused
        0x1.fffffffffffd4p-2,
        0x1.5555571d6ef9p-3,
        0x1.5555576a5adcep-5,
#elif EXP64_N == 256 && EXP64_POLY_ORDER == 5 && EXP64_POLY_WIDE
        // abs error: 1.5547*2^-66
        // ulp error: 0.505 (0.506 without fma)
        // if |x| < ln2/256
        0x1.ffffffffffdbdp-2,
        0x1.555555555543cp-3,
        0x1.55555cf16e1edp-5,
        0x1.1111167a4b553p-7,
#endif
    },
    .exp2_shift = 0x1.8p52 / EXP64_N,
    // exp2 polynomial coefficients.
    .exp2_poly = {
#if EXP64_N == 64 && EXP2_64_POLY_ORDER == 6 && EXP2_64_POLY_WIDE
        // abs error: 1.3054*2^-63
        // ulp error: 0.515
        // if |x| < 1/64
        0x1.62e42fefa39efp-1,
        0x1.ebfbdff82c58fp-3,
        0x1.c6b08d7045cf1p-5,
        0x1.3b2ab6fb8fd0ep-7,
        0x1.5d884afec48d7p-10,
        0x1.43097dc684ae1p-13,
#elif EXP64_N == 128 && EXP2_64_POLY_ORDER == 5 && !EXP2_64_POLY_WIDE
        // abs error: 1.2195*2^-65
        // ulp error: 0.507 (0.511 without fma)
        // if |x| < 1/256
        // abs error if |x| < 1/128: 1.9941*2^-56
        0x1.62e42fefa39efp-1,
        0x1.ebfbdff82c424p-3,
        0x1.c6b08d70cf4b5p-5,
        0x1.3b2abd24650ccp-7,
        0x1.5d7e09b4e3a84p-10,
#elif EXP64_N == 256 && EXP2_64_POLY_ORDER == 5 && EXP2_64_POLY_WIDE
        // abs error: 1.2195*2^-65
        // ulp error: 0.504 (0.508 without fma)
        // if |x| < 1/256
        0x1.62e42fefa39efp-1,
        0x1.ebfbdff82c424p-3,
        0x1.c6b08d70cf4b5p-5,
        0x1.3b2abd24650ccp-7,
        0x1.5d7e09b4e3a84p-10,
#endif
    },
    // 2^(k/N) ~= H[k]*(1 + T[k]) for int k in [0,N)
    // tab[2*k] = asuint64(T[k])
    // tab[2*k+1] = asuint64(H[k]) - (k << 52)/N
    .tab = {
#if EXP64_N == 64
        0x0, 0x3ff0000000000000,
        0xbc7160139cd8dc5d, 0x3fefec9a3e778061,
        0x3c8cd2523567f613, 0x3fefd9b0d3158574,
        0x3c60f74e61e6c861, 0x3fefc74518759bc8,
        0x3c979aa65d837b6d, 0x3fefb5586cf9890f,
        0x3c3ebe3d702f9cd1, 0x3fefa3ec32d3d1a2,
        0xbc9556522a2fbd0e, 0x3fef9301d0125b51,
        0xbc91c923b9d5f416, 0x3fef829aaea92de0,
        0xbc801b15eaa59348, 0x3fef72b83c7d517b,
        0x3c8b898c3f1353bf, 0x3fef635beb6fcb75,
        0x3c9aecf73e3a2f60, 0x3fef54873168b9aa,
        0x3c8a6f4144a6c38d, 0x3fef463b88628cd6,
        0x3c968efde3a8a894, 0x3fef387a6e756238,
        0x3c80472b981fe7f2, 0x3fef2b4565e27cdd,
        0x3c82f7e16d09ab31, 0x3fef1e9df51fdee1,
        0x3c8b3782720c0ab4, 0x3fef1285a6e4030b,
        0x3c834d754db0abb6, 0x3fef06fe0a31b715,
        0x3c8fdd395dd3f84a, 0x3feefc08b26416ff,
        0xbc924aedcc4b5068, 0x3feef1a7373aa9cb,
        0xbc71d1e83e9436d2, 0x3feee7db34e59ff7,
        0x3c859f48a72a4c6d, 0x3feedea64c123422,
        0xbc58a78f4817895b, 0x3feed60a21f72e2a,
        0x3c4363ed60c2ac11, 0x3feece086061892d,
        0x3c6ecce1daa10379, 0x3feec6a2b5c13cd0,
        0x3c7690cebb7aafb0, 0x3feebfdad5362a27,
        0xbc8f94340071a38e, 0x3feeb9b2769d2ca7,
        0xbc78dec6bd0f385f, 0x3feeb42b569d4f82,
        0x3c93350518fdd78e, 0x3feeaf4736b527da,
        0x3c9063e1e21c5409, 0x3feeab07dd485429,
        0x3c9432e62b64c035, 0x3feea76f15ad2148,
        0xbc8c33c53bef4da8, 0x3feea47eb03a5585,
        0xbc93cedd78565858, 0x3feea23882552225,
        0xbc93b3efbf5e2228, 0x3feea09e667f3bcd,
        0xbc6367efb86da9ee, 0x3fee9fb23c651a2f,
        0xbc781f647e5a3ecf, 0x3fee9f75e8ec5f74,
        0xbc8619321e55e68a, 0x3fee9feb564267c9,
        0xbc7b32dcb94da51d, 0x3feea11473eb0187,
        0x3c65ebe1abd66c55, 0x3feea2f336cf4e62,
        0xbc9369b6f13b3734, 0x3feea589994cce13,
        0xbc94d450d872576e, 0x3feea8d99b4492ed,
        0x3c8db72fc1f0eab4, 0x3feeace5422aa0db,
        0x3c7bf68359f35f44, 0x3feeb1ae99157736,
        0xbc5da9b88b6c1e29, 0x3feeb737b0cdc5e5,
        0xbc92434322f4f9aa, 0x3feebd829fde4e50,
        0x3c71affc2b91ce27, 0x3feec49182a3f090,
        0xbc87c50422622263, 0x3feecc667b5de565,
        0xbc91bbd1d3bcbb15, 0x3feed503b23e255d,
        0x3c8469846e735ab3, 0x3feede6b5579fdbf,
        0x3c8c1a7792cb3387, 0x3feee89f995ad3ad,
        0xbc55c3d956dcaeba, 0x3feef3a2b84f15fb,
        0xbc68d6f438ad9334, 0x3feeff76f2fb5e47,
        0x3c74ffd70a5fddcd, 0x3fef0c1e904bc1d2,
        0x3c736eae30af0cb3, 0x3fef199bdd85529c,
        0x3c84e08fd10959ac, 0x3fef27f12e57d14b,
        0x3c676b2c6c921968, 0x3fef3720dcef9069,
        0xbc8fad5d3ffffa6f, 0x3fef472d4a07897c,
        0x3c74a385a63d07a7, 0x3fef5818dcfba487,
        0x3c8e5a50d5c192ac, 0x3fef69e603db3285,
        0xbc82d52107b43e1f, 0x3fef7c97337b9b5f,
        0x3c74b604603a88d3, 0x3fef902ee78b3ff6,
        0xbc8ff7128fd391f0, 0x3fefa4afa2a490da,
        0x3c8ec3bc41aa2008, 0x3fefba1bee615a27,
        0x3c8a64a931d185ee, 0x3fefd0765b6e4540,
        0x3c77893b4d91cd9d, 0x3fefe7c1819e90d8,
#elif EXP64_N == 128
        0x0, 0x3ff0000000000000,
        0x3c9b3b4f1a88bf6e, 0x3feff63da9fb3335,
        0xbc7160139cd8dc5d, 0x3fefec9a3e778061,
        0xbc905e7a108766d1, 0x3fefe315e86e7f85,
        0x3c8cd2523567f613, 0x3fefd9b0d3158574,
        0xbc8bce8023f98efa, 0x3fefd06b29ddf6de,
        0x3c60f74e61e6c861, 0x3fefc74518759bc8,
        0x3c90a3e45b33d399, 0x3fefbe3ecac6f383,
        0x3c979aa65d837b6d, 0x3fefb5586cf9890f,
        0x3c8eb51a92fdeffc, 0x3fefac922b7247f7,
        0x3c3ebe3d702f9cd1, 0x3fefa3ec32d3d1a2,
        0xbc6a033489906e0b, 0x3fef9b66affed31b,
        0xbc9556522a2fbd0e, 0x3fef9301d0125b51,
        0xbc5080ef8c4eea55, 0x3fef8abdc06c31cc,
        0xbc91c923b9d5f416, 0x3fef829aaea92de0,
        0x3c80d3e3e95c55af, 0x3fef7a98c8a58e51,
        0xbc801b15eaa59348, 0x3fef72b83c7d517b,
        0xbc8f1ff055de323d, 0x3fef6af9388c8dea,
        0x3c8b898c3f1353bf, 0x3fef635beb6fcb75,
        0xbc96d99c7611eb26, 0x3fef5be084045cd4,
        0x3c9aecf73e3a2f60, 0x3fef54873168b9aa,
        0xbc8fe782cb86389d, 0x3fef4d5022fcd91d,
        0x3c8a6f4144a6c38d, 0x3fef463b88628cd6,
        0x3c807a05b0e4047d, 0x3fef3f49917ddc96,
        0x3c968efde3a8a894, 0x3fef387a6e756238,
        0x3c875e18f274487d, 0x3fef31ce4fb2a63f,
        0x3c80472b981fe7f2, 0x3fef2b4565e27cdd,
        0xbc96b87b3f71085e, 0x3fef24dfe1f56381,
        0x3c82f7e16d09ab31, 0x3fef1e9df51fdee1,
        0xbc3d219b1a6fbffa, 0x3fef187fd0dad990,
        0x3c8b3782720c0ab4, 0x3fef1285a6e4030b,
        0x3c6e149289cecb8f, 0x3fef0cafa93e2f56,
        0x3c834d754db0abb6, 0x3fef06fe0a31b715,
        0x3c864201e2ac744c, 0x3fef0170fc4cd831,
        0x3c8fdd395dd3f84a, 0x3feefc08b26416ff,
        0xbc86a3803b8e5b04, 0x3feef6c55f929ff1,
        0xbc924aedcc4b5068, 0x3feef1a7373aa9cb,
        0xbc9907f81b512d8e, 0x3feeecae6d05d866,
        0xbc71d1e83e9436d2, 0x3feee7db34e59ff7,
        0xbc991919b3ce1b15, 0x3feee32dc313a8e5,
        0x3c859f48a72a4c6d, 0x3feedea64c123422,
        0xbc9312607a28698a, 0x3feeda4504ac801c,
        0xbc58a78f4817895b, 0x3feed60a21f72e2a,
        0xbc7c2c9b67499a1b, 0x3feed1f5d950a897,
        0x3c4363ed60c2ac11, 0x3feece086061892d,
        0x3c9666093b0664ef, 0x3feeca41ed1d0057,
        0x3c6ecce1daa10379, 0x3feec6a2b5c13cd0,
        0x3c93ff8e3f0f1230, 0x3feec32af0d7d3de,
        0x3c7690cebb7aafb0, 0x3feebfdad5362a27,
        0x3c931dbdeb54e077, 0x3feebcb299fddd0d,
        0xbc8f94340071a38e, 0x3feeb9b2769d2ca7,
        0xbc87deccdc93a349, 0x3feeb6daa2cf6642,
        0xbc78dec6bd0f385f, 0x3feeb42b569d4f82,
        0xbc861246ec7b5cf6, 0x3feeb1a4ca5d920f,
        0x3c93350518fdd78e, 0x3feeaf4736b527da,
        0x3c7b98b72f8a9b05, 0x3feead12d497c7fd,
        0x3c9063e1e21c5409, 0x3feeab07dd485429,
        0x3c34c7855019c6ea, 0x3feea9268a5946b7,
        0x3c9432e62b64c035, 0x3feea76f15ad2148,
        0xbc8ce44a6199769f, 0x3feea5e1b976dc09,
        0xbc8c33c53bef4da8, 0x3feea47eb03a5585,
        0xbc845378892be9ae, 0x3feea34634ccc320,
        0xbc93cedd78565858, 0x3feea23882552225,
        0x3c5710aa807e1964, 0x3feea155d44ca973,
        0xbc93b3efbf5e2228, 0x3feea09e667f3bcd,
        0xbc6a12ad8734b982, 0x3feea012750bdabf,
        0xbc6367efb86da9ee, 0x3fee9fb23c651a2f,
        0xbc80dc3d54e08851, 0x3fee9f7df9519484,
        0xbc781f647e5a3ecf, 0x3fee9f75e8ec5f74,
        0xbc86ee4ac08b7db0, 0x3fee9f9a48a58174,
        0xbc8619321e55e68a, 0x3fee9feb564267c9,
        0x3c909ccb5e09d4d3, 0x3feea0694fde5d3f,
        0xbc7b32dcb94da51d, 0x3feea11473eb0187,
        0x3c94ecfd5467c06b, 0x3feea1ed0130c132,
        0x3c65ebe1abd66c55, 0x3feea2f336cf4e62,
        0xbc88a1c52fb3cf42, 0x3feea427543e1a12,
        0xbc9369b6f13b3734, 0x3feea589994cce13,
        0xbc805e843a19ff1e, 0x3feea71a4623c7ad,
        0xbc94d450d872576e, 0x3feea8d99b4492ed,
        0x3c90ad675b0e8a00, 0x3feeaac7d98a6699,
        0x3c8db72fc1f0eab4, 0x3feeace5422aa0db,
        0xbc65b6609cc5e7ff, 0x3feeaf3216b5448c,
        0x3c7bf68359f35f44, 0x3feeb1ae99157736,
        0xbc93091fa71e3d83, 0x3feeb45b0b91ffc6,
        0xbc5da9b88b6c1e29, 0x3feeb737b0cdc5e5,
        0xbc6c23f97c90b959, 0x3feeba44cbc8520f,
        0xbc92434322f4f9aa, 0x3feebd829fde4e50,
        0xbc85ca6cd7668e4b, 0x3feec0f170ca07ba,
        0x3c71affc2b91ce27, 0x3feec49182a3f090,
        0x3c6dd235e10a73bb, 0x3feec86319e32323,
        0xbc87c50422622263, 0x3feecc667b5de565,
        0x3c8b1c86e3e231d5, 0x3feed09bec4a2d33,
        0xbc91bbd1d3bcbb15, 0x3feed503b23e255d,
        0x3c90cc319cee31d2, 0x3feed99e1330b358,
        0x3c8469846e735ab3, 0x3feede6b5579fdbf,
        0xbc82dfcd978e9db4, 0x3feee36bbfd3f37a,
        0x3c8c1a7792cb3387, 0x3feee89f995ad3ad,
        0xbc907b8f4ad1d9fa, 0x3feeee07298db666,
        0xbc55c3d956dcaeba, 0x3feef3a2b84f15fb,
        0xbc90a40e3da6f640, 0x3feef9728de5593a,
        0xbc68d6f438ad9334, 0x3feeff76f2fb5e47,
        0xbc91eee26b588a35, 0x3fef05b030a1064a,
        0x3c74ffd70a5fddcd, 0x3fef0c1e904bc1d2,
        0xbc91bdfbfa9298ac, 0x3fef12c25bd71e09,
        0x3c736eae30af0cb3, 0x3fef199bdd85529c,
        0x3c8ee3325c9ffd94, 0x3fef20ab5fffd07a,
        0x3c84e08fd10959ac, 0x3fef27f12e57d14b,
        0x3c63cdaf384e1a67, 0x3fef2f6d9406e7b5,
        0x3c676b2c6c921968, 0x3fef3720dcef9069,
        0xbc808a1883ccb5d2, 0x3fef3f0b555dc3fa,
        0xbc8fad5d3ffffa6f, 0x3fef472d4a07897c,
        0xbc900dae3875a949, 0x3fef4f87080d89f2,
        0x3c74a385a63d07a7, 0x3fef5818dcfba487,
        0xbc82919e2040220f, 0x3fef60e316c98398,
        0x3c8e5a50d5c192ac, 0x3fef69e603db3285,
        0x3c843a59ac016b4b, 0x3fef7321f301b460,
        0xbc82d52107b43e1f, 0x3fef7c97337b9b5f,
        0xbc892ab93b470dc9, 0x3fef864614f5a129,
        0x3c74b604603a88d3, 0x3fef902ee78b3ff6,
        0x3c83c5ec519d7271, 0x3fef9a51fbc74c83,
        0xbc8ff7128fd391f0, 0x3fefa4afa2a490da,
        0xbc8dae98e223747d, 0x3fefaf482d8e67f1,
        0x3c8ec3bc41aa2008, 0x3fefba1bee615a27,
        0x3c842b94c3a9eb32, 0x3fefc52b376bba97,
        0x3c8a64a931d185ee, 0x3fefd0765b6e4540,
        0xbc8e37bae43be3ed, 0x3fefdbfdad9cbe14,
        0x3c77893b4d91cd9d, 0x3fefe7c1819e90d8,
        0x3c5305c14160cc89, 0x3feff3c22b8f71f1,
#elif EXP64_N == 256
        0x0, 0x3ff0000000000000,
        0xbc84e82fc61851ac, 0x3feffb1afa5abcbf,
        0x3c9b3b4f1a88bf6e, 0x3feff63da9fb3335,
        0xbc82985dd8521d32, 0x3feff168143b0281,
        0xbc7160139cd8dc5d, 0x3fefec9a3e778061,
        0x3c651e617061bfbd, 0x3fefe7d42e11bbcc,
        0xbc905e7a108766d1, 0x3fefe315e86e7f85,
        0x3c845fad437fa426, 0x3fefde5f72f654b1,
        0x3c8cd2523567f613, 0x3fefd9b0d3158574,
        0xbc954529642b232f, 0x3fefd50a0e3c1f89,
        0xbc8bce8023f98efa, 0x3fefd06b29ddf6de,
        0x3c8293708ef5c32e, 0x3fefcbd42b72a836,
        0x3c60f74e61e6c861, 0x3fefc74518759bc8,
        0xbc95b9280905b2a4, 0x3fefc2bdf66607e0,
        0x3c90a3e45b33d399, 0x3fefbe3ecac6f383,
        0x3c84f31f32c4b7e7, 0x3fefb9c79b1f3919,
        0x3c979aa65d837b6d, 0x3fefb5586cf9890f,
        0x3c9407fb30d06420, 0x3fefb0f145e46c85,
        0x3c8eb51a92fdeffc, 0x3fefac922b7247f7,
        0xbc9a5d04b3b9911b, 0x3fefa83b23395dec,
        0x3c3ebe3d702f9cd1, 0x3fefa3ec32d3d1a2,
        0xbc937a01f0739546, 0x3fef9fa55fdfa9c5,
        0xbc6a033489906e0b, 0x3fef9b66affed31b,
        0x3c8b8268b04ef0a5, 0x3fef973028d7233e,
        0xbc9556522a2fbd0e, 0x3fef9301d0125b51,
        0xbc9ac46e44a2ebcc, 0x3fef8edbab5e2ab6,
        0xbc5080ef8c4eea55, 0x3fef8abdc06c31cc,
        0xbc65704e90c9f860, 0x3fef86a814f204ab,
        0xbc91c923b9d5f416, 0x3fef829aaea92de0,
        0xbc897cea57e46280, 0x3fef7e95934f312e,
        0x3c80d3e3e95c55af, 0x3fef7a98c8a58e51,
        0x3c56f01429e2b9d2, 0x3fef76a45471c3c2,
        0xbc801b15eaa59348, 0x3fef72b83c7d517b,
        0x3c6e653b2459034b, 0x3fef6ed48695bbc0,
        0xbc8f1ff055de323d, 0x3fef6af9388c8dea,
        0x3c92cc7ea345b7dc, 0x3fef672658375d2f,
        0x3c8b898c3f1353bf, 0x3fef635beb6fcb75,
        0x3c957bfb2876ea9e, 0x3fef5f99f8138a1c,
        0xbc96d99c7611eb26, 0x3fef5be084045cd4,
        0x3c8cdc1873af2155, 0x3fef582f95281c6b,
        0x3c9aecf73e3a2f60, 0x3fef54873168b9aa,
        0xbc9493684653a131, 0x3fef50e75eb44027,
        0xbc8fe782cb86389d, 0x3fef4d5022fcd91d,
        0xbc98e2899077520a, 0x3fef49c18438ce4d,
        0x3c8a6f4144a6c38d, 0x3fef463b88628cd6,
        0x3c9120fcd4f59273, 0x3fef42be3578a819,
        0x3c807a05b0e4047d, 0x3fef3f49917ddc96,
        0x3c89b788c188c9b8, 0x3fef3bdda27912d1,
        0x3c968efde3a8a894, 0x3fef387a6e756238,
        0x3c877afbca90ef84, 0x3fef351ffb82140a,
        0x3c875e18f274487d, 0x3fef31ce4fb2a63f,
        0x3c91512f082876ee, 0x3fef2e85711ece75,
        0x3c80472b981fe7f2, 0x3fef2b4565e27cdd,
        0x3c9a02f0c7d75ec6, 0x3fef280e341ddf29,
        0xbc96b87b3f71085e, 0x3fef24dfe1f56381,
        0xbc803297e78260bf, 0x3fef21ba7591bb70,
        0x3c82f7e16d09ab31, 0x3fef1e9df51fdee1,
        0xbc95b77e5ccd9fbf, 0x3fef1b8a66d10f13,
        0xbc3d219b1a6fbffa, 0x3fef187fd0dad990,
        0xbc91e75c40b4251e, 0x3fef157e39771b2f,
        0x3c8b3782720c0ab4, 0x3fef1285a6e4030b,
        0x3c98a911f1f7785a, 0x3fef0f961f641589,
        0x3c6e149289cecb8f, 0x3fef0cafa93e2f56,
        0xbc61e7c998db7dbb, 0x3fef09d24abd886b,
        0x3c834d754db0abb6, 0x3fef06fe0a31b715,
        0x3c85425c11faadf4, 0x3fef0432edeeb2fd,
        0x3c864201e2ac744c, 0x3fef0170fc4cd831,
        0xbc979517a03e2847, 0x3feefeb83ba8ea32,
        0x3c8fdd395dd3f84a, 0x3feefc08b26416ff,
        0xbc800e2a46da4bee, 0x3feef96266e3fa2d,
        0xbc86a3803b8e5b04, 0x3feef6c55f929ff1,
        0xbc87430803972b34, 0x3feef431a2de883b,
        0xbc924aedcc4b5068, 0x3feef1a7373aa9cb,
        0xbc954de30ae02d94, 0x3feeef26231e754a,
        0xbc9907f81b512d8e, 0x3feeecae6d05d866,
        0xbc94f2487e1c03ec, 0x3feeea401b7140ef,
        0xbc71d1e83e9436d2, 0x3feee7db34e59ff7,
        0x3c914a5432fcb2f4, 0x3feee57fbfec6cf4,
        0xbc991919b3ce1b15, 0x3feee32dc313a8e5,
        0x3c79c3bba5562a2f, 0x3feee0e544ede173,
        0x3c859f48a72a4c6d, 0x3feedea64c123422,
        0xbc85a71612e21658, 0x3feedc70df1c5175,
        0xbc9312607a28698a, 0x3feeda4504ac801c,
        0x3c86421f6f1d24d6, 0x3feed822c367a024,
        0xbc58a78f4817895b, 0x3feed60a21f72e2a,
        0xbc9348a6815fce65, 0x3feed3fb2709468a,
        0xbc7c2c9b67499a1b, 0x3feed1f5d950a897,
        0x3c835c43984d9871, 0x3feecffa3f84b9d4,
        0x3c4363ed60c2ac11, 0x3feece086061892d,
        0xbc632afc8d9473a0, 0x3feecc2042a7d232,
        0x3c9666093b0664ef, 0x3feeca41ed1d0057,
        0xbc95fc5e44de020e, 0x3feec86d668b3237,
        0x3c6ecce1daa10379, 0x3feec6a2b5c13cd0,
        0xbc7ea0148327c42f, 0x3feec4e1e192aed2,
        0x3c93ff8e3f0f1230, 0x3feec32af0d7d3de,
        0xbc7a843ad1a88022, 0x3feec17dea6db7d7,
        0x3c7690cebb7aafb0, 0x3feebfdad5362a27,
        0x3c892ca3bf144e63, 0x3feebe41b817c114,
        0x3c931dbdeb54e077, 0x3feebcb299fddd0d,
        0xbc902c99b04aa8b0, 0x3feebb2d81d8abff,
        0xbc8f94340071a38e, 0x3feeb9b2769d2ca7,
        0x3c73e34f67e67118, 0x3feeb8417f4531ee,
        0xbc87deccdc93a349, 0x3feeb6daa2cf6642,
        0xbc75a3b1197ba0f0, 0x3feeb57de83f4eef,
        0xbc78dec6bd0f385f, 0x3feeb42b569d4f82,
        0x3c81bd2888075068, 0x3feeb2e2f4f6ad27,
        0xbc861246ec7b5cf6, 0x3feeb1a4ca5d920f,
        0xbc896be8ae89ef8f, 0x3feeb070dde910d2,
        0x3c93350518fdd78e, 0x3feeaf4736b527da,
        0xbc88e6ac90348602, 0x3feeae27dbe2c4cf,
        0x3c7b98b72f8a9b05, 0x3feead12d497c7fd,
        0xbc91af7f1365c3ac, 0x3feeac0827ff07cc,
        0x3c9063e1e21c5409, 0x3feeab07dd485429,
        0xbc943a3540d1898a, 0x3feeaa11fba87a03,
        0x3c34c7855019c6ea, 0x3feea9268a5946b7,
        0xbc951f58ddaa8090, 0x3feea84590998b93,
        0x3c9432e62b64c035, 0x3feea76f15ad2148,
        0xbc82e1648e50a17c, 0x3feea6a320dceb71,
        0xbc8ce44a6199769f, 0x3feea5e1b976dc09,
        0x3c95f30eda98a575, 0x3feea52ae6cdf6f4,
        0xbc8c33c53bef4da8, 0x3feea47eb03a5585,
        0x3c917ecda8a72159, 0x3feea3dd1d1929fd,
        0xbc845378892be9ae, 0x3feea34634ccc320,
        0xbc9345f3cee1ae6e, 0x3feea2b9febc8fb7,
        0xbc93cedd78565858, 0x3feea23882552225,
        0xbc85c33fdf910406, 0x3feea1c1c70833f6,
        0x3c5710aa807e1964, 0x3feea155d44ca973,
        0x3c81079ab5789604, 0x3feea0f4b19e9538,
        0xbc93b3efbf5e2228, 0x3feea09e667f3bcd,
        0x3c727df161cd7778, 0x3feea052fa75173e,
        0xbc6a12ad8734b982, 0x3feea012750bdabf,
        0x3c93f9924a05b767, 0x3fee9fdcddd47645,
        0xbc6367efb86da9ee, 0x3fee9fb23c651a2f,
        0xbc87557939a8b5ef, 0x3fee9f9298593ae5,
        0xbc80dc3d54e08851, 0x3fee9f7df9519484,
        0x3c51ed2f56fa9d1a, 0x3fee9f7466f42e87,
        0xbc781f647e5a3ecf, 0x3fee9f75e8ec5f74,
        0xbc88e67a9006c909, 0x3fee9f8286ead08a,
        0xbc86ee4ac08b7db0, 0x3fee9f9a48a58174,
        0x3c86597566977ac8, 0x3fee9fbd35d7cbfd,
        0xbc8619321e55e68a, 0x3fee9feb564267c9,
        0x3c92c0b7028a5c3a, 0x3feea024b1ab6e09,
        0x3c909ccb5e09d4d3, 0x3feea0694fde5d3f,
        0x3c8a30faf49cc78c, 0x3feea0b938ac1cf6,
        0xbc7b32dcb94da51d, 0x3feea11473eb0187,
        0xbc92dad3519d7b5b, 0x3feea17b0976cfdb,
        0x3c94ecfd5467c06b, 0x3feea1ed0130c132,
        0x3c87d51410fd15c2, 0x3feea26a62ff86f0,
        0x3c65ebe1abd66c55, 0x3feea2f336cf4e62,
        0xbc760a3629969871, 0x3feea3878491c491,
        0xbc88a1c52fb3cf42, 0x3feea427543e1a12,
        0x3c8b18c6e3fdef5d, 0x3feea4d2add106d9,
        0xbc9369b6f13b3734, 0x3feea589994cce13,
        0x3c90ec1ddcb1390a, 0x3feea64c1eb941f7,
        0xbc805e843a19ff1e, 0x3feea71a4623c7ad,
        0xbc522cea4f3afa1e, 0x3feea7f4179f5b21,
        0xbc94d450d872576e, 0x3feea8d99b4492ed,
        0x3c7c88549b958471, 0x3feea9cad931a436,
        0x3c90ad675b0e8a00, 0x3feeaac7d98a6699,
        0x3c931143962f7877, 0x3feeabd0a478580f,
        0x3c8db72fc1f0eab4, 0x3feeace5422aa0db,
        0x3c93e9e96f112479, 0x3feeae05bad61778,
        0xbc65b6609cc5e7ff, 0x3feeaf3216b5448c,
        0xbc8dac42a4a38df0, 0x3feeb06a5e0866d9,
        0x3c7bf68359f35f44, 0x3feeb1ae99157736,
        0x3c8b99dd98b1ed84, 0x3feeb2fed0282c8a,
        0xbc93091fa71e3d83, 0x3feeb45b0b91ffc6,
        0xbc7885ad50cbb750, 0x3feeb5c353aa2fe2,
        0xbc5da9b88b6c1e29, 0x3feeb737b0cdc5e5,
        0xbc82d5e85f3e0301, 0x3feeb8b82b5f98e5,
        0xbc6c23f97c90b959, 0x3feeba44cbc8520f,
        0xbc51669428996971, 0x3feebbdd9a7670b3,
        0xbc92434322f4f9aa, 0x3feebd829fde4e50,
        0x3c71f2b2c1c4c014, 0x3feebf33e47a22a2,
        0xbc85ca6cd7668e4b, 0x3feec0f170ca07ba,
        0xbc9294f304f166b6, 0x3feec2bb4d53fe0d,
        0x3c71affc2b91ce27, 0x3feec49182a3f090,
        0xbc8a1e58414c07d3, 0x3feec674194bb8d5,
        0x3c6dd235e10a73bb, 0x3feec86319e32323,
        0xbc79740b58a20091, 0x3feeca5e8d07f29e,
        0xbc87c50422622263, 0x3feecc667b5de565,
        0x3c9165830a2b96c2, 0x3feece7aed8eb8bb,
        0x3c8b1c86e3e231d5, 0x3feed09bec4a2d33,
        0xbc903d5cbe27874b, 0x3feed2c980460ad8,
        0xbc91bbd1d3bcbb15, 0x3feed503b23e255d,
        0x3c5986178980fce0, 0x3feed74a8af46052,
        0x3c90cc319cee31d2, 0x3feed99e1330b358,
        0xbc89472975b1f2a5, 0x3feedbfe53c12e59,
        0x3c8469846e735ab3, 0x3feede6b5579fdbf,
        0x3c7d8157a34b7e7f, 0x3feee0e521356eba,
        0xbc82dfcd978e9db4, 0x3feee36bbfd3f37a,
        0x3c8c8a4e231ebb7d, 0x3feee5ff3a3c2774,
        0x3c8c1a7792cb3387, 0x3feee89f995ad3ad,
        0xbc888c8d11a142e5, 0x3feeeb4ce622f2ff,
        0xbc907b8f4ad1d9fa, 0x3feeee07298db666,
        0x3c889c2ea41433c7, 0x3feef0ce6c9a8952,
        0xbc55c3d956dcaeba, 0x3feef3a2b84f15fb,
        0xbc7274aedac8ff80, 0x3feef68415b749b1,
        0xbc90a40e3da6f640, 0x3feef9728de5593a,
        0x3c85c620ce76df06, 0x3feefc6e29f1c52a,
        0xbc68d6f438ad9334, 0x3feeff76f2fb5e47,
        0xbc8fda52e1b51e41, 0x3fef028cf22749e4,
        0xbc91eee26b588a35, 0x3fef05b030a1064a,
        0xbc32141a7b3e2cd8, 0x3fef08e0b79a6f1f,
        0x3c74ffd70a5fddcd, 0x3fef0c1e904bc1d2,
        0xbc302899507554e5, 0x3fef0f69c3f3a207,
        0xbc91bdfbfa9298ac, 0x3fef12c25bd71e09,
        0xbc80dda2d4c0010c, 0x3fef16286141b33d,
        0x3c736eae30af0cb3, 0x3fef199bdd85529c,
        0xbc8a007daadf8d68, 0x3fef1d1cd9fa652c,
        0x3c8ee3325c9ffd94, 0x3fef20ab5fffd07a,
        0x3c836909391181d3, 0x3fef244778fafb22,
        0x3c84e08fd10959ac, 0x3fef27f12e57d14b,
        0xbc811cd7dbdf9547, 0x3fef2ba88988c933,
        0x3c63cdaf384e1a67, 0x3fef2f6d9406e7b5,
        0xbc7ac28b7bef6621, 0x3fef33405751c4db,
        0x3c676b2c6c921968, 0x3fef3720dcef9069,
        0xbc7030587207b9e1, 0x3fef3b0f2e6d1675,
        0xbc808a1883ccb5d2, 0x3fef3f0b555dc3fa,
        0xbc8cc734592af7fc, 0x3fef43155b5bab74,
        0xbc8fad5d3ffffa6f, 0x3fef472d4a07897c,
        0x3c87752a44f587e8, 0x3fef4b532b08c968,
        0xbc900dae3875a949, 0x3fef4f87080d89f2,
        0x3c85b66fefeef52e, 0x3fef53c8eacaa1d6,
        0x3c74a385a63d07a7, 0x3fef5818dcfba487,
        0x3c5159d9d908a96e, 0x3fef5c76e862e6d3,
        0xbc82919e2040220f, 0x3fef60e316c98398,
        0x3c8c254d16117a68, 0x3fef655d71ff6075,
        0x3c8e5a50d5c192ac, 0x3fef69e603db3285,
        0xbc8d8c329fbd0e03, 0x3fef6e7cd63a8315,
        0x3c843a59ac016b4b, 0x3fef7321f301b460,
        0xbc8ea6e6fbd5f2a6, 0x3fef77d5641c0658,
        0xbc82d52107b43e1f, 0x3fef7c97337b9b5f,
        0xbc63e8e3eab2cbb4, 0x3fef81676b197d17,
        0xbc892ab93b470dc9, 0x3fef864614f5a129,
        0xbc8b7966cd0d2cd9, 0x3fef8b333b16ee12,
        0x3c74b604603a88d3, 0x3fef902ee78b3ff6,
        0xbc776caa4c2ff1cf, 0x3fef953924676d76,
        0x3c83c5ec519d7271, 0x3fef9a51fbc74c83,
        0xbc81d5fc525d9940, 0x3fef9f7977cdb740,
        0xbc8ff7128fd391f0, 0x3fefa4afa2a490da,
        0x3c855cd8aaea3d21, 0x3fefa9f4867cca6e,
        0xbc8dae98e223747d, 0x3fefaf482d8e67f1,
        0x3c8269947c2bed4a, 0x3fefb4aaa2188510,
        0x3c8ec3bc41aa2008, 0x3fefba1bee615a27,
        0xbc83b6137e9afe9e, 0x3fefbf9c1cb6412a,
        0x3c842b94c3a9eb32, 0x3fefc52b376bba97,
        0xbc69fa74878ba7c7, 0x3fefcac948dd7274,
        0x3c8a64a931d185ee, 0x3fefd0765b6e4540,
        0x3c901f3a75ee0efe, 0x3fefd632798844f8,
        0xbc8e37bae43be3ed, 0x3fefdbfdad9cbe14,
        0xbc516a9ce6ed84fa, 0x3fefe1d802243c89,
        0x3c77893b4d91cd9d, 0x3fefe7c1819e90d8,
        0xbc699c7db2effc76, 0x3fefedba3692d514,
        0x3c5305c14160cc89, 0x3feff3c22b8f71f1,
        0x3c64b458677f9840, 0x3feff9d96b2a23d9,
#endif
    },
};
#undef EXP64_N


#define LOG64_N (1 << LOG64_TABLE_BITS)
global const LogData gLogData64 = {
    .ln2hi = 0x1.62e42fefa3800p-1,
    .ln2lo = 0x1.ef35793c76730p-45,
    .poly = {
#if LOG64_N == 64 && LOG64_POLY_ORDER == 7
        // relative error: 0x1.906eb8ap-58
        // abs error: 0x1.d2cad5a8p-67
        // in -0x1.fp-8 0x1.fp-8
        -0x1.0000000000027p-1,
        0x1.555555555556ap-2,
        -0x1.fffffff0440bap-3,
        0x1.99999991906c3p-3,
        -0x1.555c8d7e8201ep-3,
        0x1.24978c59151fap-3,
#elif LOG64_N == 128 && LOG64_POLY_ORDER == 6
        // relative error: 0x1.926199e8p-56
        // abs error: 0x1.882ff33p-65
        // in -0x1.fp-9 0x1.fp-9
        -0x1.0000000000001p-1,
        0x1.555555551305bp-2,
        -0x1.fffffffeb459p-3,
        0x1.999b324f10111p-3,
        -0x1.55575e506c89fp-3,
#elif LOG64_N == 128 && LOG64_POLY_ORDER == 7
        // relative error: 0x1.649fc4bp-64
        // abs error: 0x1.c3b5769p-74
        // in -0x1.fp-9 0x1.fp-9
        -0x1.0000000000001p-1,
        0x1.5555555555556p-2,
        -0x1.fffffffea1a8p-3,
        0x1.99999998e9139p-3,
        -0x1.555776801b968p-3,
        0x1.2493c29331a5cp-3,
#endif
    },
    .poly1 = {
#if LOG64_POLY1_ORDER == 10
        // relative error: 0x1.32eccc6p-62
        // in -0x1p-5 0x1.1p-5 (|log(1+x)| > 0x1p-5 outside this interval)
        -0x1p-1,
        0x1.55555555554e5p-2,
        -0x1.0000000000af2p-2,
        0x1.9999999bbe436p-3,
        -0x1.55555537f9cdep-3,
        0x1.24922fc8127cfp-3,
        -0x1.0000b7d6bb612p-3,
        0x1.c806ee1ddbcafp-4,
        -0x1.972335a9c2d6ep-4,
#elif LOG64_POLY1_ORDER == 11
        // relative error: 0x1.52c8b708p-68
        // in -0x1p-5 0x1.1p-5 (|log(1+x)| > 0x1p-5 outside this interval)
        -0x1p-1,
        0x1.5555555555555p-2,
        -0x1.ffffffffffea9p-3,
        0x1.999999999c4d4p-3,
        -0x1.55555557f5541p-3,
        0x1.249248fbe33e4p-3,
        -0x1.ffffc9a3c825bp-4,
        0x1.c71e1f204435dp-4,
        -0x1.9a7f26377d06ep-4,
        0x1.71c30cf8f7364p-4,
#elif LOG64_POLY1_ORDER == 12
        // relative error: 0x1.c04d76cp-63
        // in -0x1p-4 0x1.09p-4 (|log(1+x)| > 0x1p-4 outside the interval)
        -0x1p-1,
        0x1.5555555555577p-2,
        -0x1.ffffffffffdcbp-3,
        0x1.999999995dd0cp-3,
        -0x1.55555556745a7p-3,
        0x1.24924a344de3p-3,
        -0x1.fffffa4423d65p-4,
        0x1.c7184282ad6cap-4,
        -0x1.999eb43b068ffp-4,
        0x1.78182f7afd085p-4,
        -0x1.5521375d145cdp-4,
#endif
    },

    /* Algorithm:

        x = 2^k z
        log(x) = k ln2 + log(c) + log(z/c)
        log(z/c) = poly(z/c - 1)

    where z is in [1.6p-1; 1.6p0] which is split into N subintervals and z falls
    into the ith one, then table entries are computed as

        tab[i].invc = 1/c
        tab[i].logc = (double)log(c)
        tab2[i].chi = (double)c
        tab2[i].clo = (double)(c - (double)c)

    where c is near the center of the subinterval and is chosen by trying +-2^29
    floating point invc candidates around 1/center and selecting one for which

        1) the rounding error in 0x1.8p9 + logc is 0,
        2) the rounding error in z - chi - clo is < 0x1p-66 and
        3) the rounding error in (double)log(c) is minimized (< 0x1p-66).

    Note: 1) ensures that k*ln2hi + logc can be computed without rounding error,
    2) ensures that z/c - 1 can be computed as (z - chi - clo)*invc with close to
    a single rounding error when there is no fast fma for z*invc - 1, 3) ensures
    that logc + poly(z/c - 1) has small error, however near x == 1 when
    |log(x)| < 0x1p-4, this is not enough so that is special cased.  */
    .tab = {
#if LOG64_N == 64
        {0x1.7242886495cd8p+0, -0x1.79e267bdfe000p-2},
        {0x1.6e1f769340dc9p+0, -0x1.6e60ee0ecb000p-2},
        {0x1.6a13ccc8f195cp+0, -0x1.63002fdbf6000p-2},
        {0x1.661ec72e86f3ap+0, -0x1.57bf76c597000p-2},
        {0x1.623fa6c447b16p+0, -0x1.4c9e07f0d2000p-2},
        {0x1.5e75bbca31702p+0, -0x1.419b42f027000p-2},
        {0x1.5ac05655adb10p+0, -0x1.36b67660e6000p-2},
        {0x1.571ed3e940191p+0, -0x1.2bef0839e4800p-2},
        {0x1.539094ac0fbbfp+0, -0x1.21445727cb000p-2},
        {0x1.5015007e7fc42p+0, -0x1.16b5ca3c3d000p-2},
        {0x1.4cab877c31cf9p+0, -0x1.0c42d3805f800p-2},
        {0x1.49539e76a88d3p+0, -0x1.01eae61b60800p-2},
        {0x1.460cbc12211dap+0, -0x1.ef5adb9fb0000p-3},
        {0x1.42d6624debe3ap+0, -0x1.db13daab99000p-3},
        {0x1.3fb0144f0d462p+0, -0x1.c6ffbe896e000p-3},
        {0x1.3c995a1f9a9b4p+0, -0x1.b31d84722d000p-3},
        {0x1.3991c23952500p+0, -0x1.9f6c3cf6eb000p-3},
        {0x1.3698df35eaa14p+0, -0x1.8beafe7f13000p-3},
        {0x1.33ae463091760p+0, -0x1.7898db878d000p-3},
        {0x1.30d190aae3d72p+0, -0x1.6574efe4ec000p-3},
        {0x1.2e025c9203c89p+0, -0x1.527e620845000p-3},
        {0x1.2b404a7244988p+0, -0x1.3fb457d798000p-3},
        {0x1.288b01dc19544p+0, -0x1.2d1615a077000p-3},
        {0x1.25e2268085f69p+0, -0x1.1aa2b431e5000p-3},
        {0x1.23456812abb74p+0, -0x1.08598f1d2b000p-3},
        {0x1.20b4703174157p+0, -0x1.ec738fee40000p-4},
        {0x1.1e2ef308b4e9bp+0, -0x1.c885768862000p-4},
        {0x1.1bb4a36b70a3fp+0, -0x1.a4e75b6a46000p-4},
        {0x1.194538e960658p+0, -0x1.8197efba9a000p-4},
        {0x1.16e0692a10ac8p+0, -0x1.5e95ad734e000p-4},
        {0x1.1485f1ba1568bp+0, -0x1.3bdf67117c000p-4},
        {0x1.12358e123ed6fp+0, -0x1.1973b744f0000p-4},
        {0x1.0fef01de37c8dp+0, -0x1.eea33446bc000p-5},
        {0x1.0db20b82be414p+0, -0x1.aaef4ab304000p-5},
        {0x1.0b7e6f67f69b3p+0, -0x1.67c962fd2c000p-5},
        {0x1.0953f342fc108p+0, -0x1.252f29acf8000p-5},
        {0x1.0732604ec956bp+0, -0x1.c63d19e9c0000p-6},
        {0x1.051980117f9b0p+0, -0x1.432ab6a388000p-6},
        {0x1.03091aa6810f1p+0, -0x1.8244357f50000p-7},
        {0x1.01010152cf066p+0, -0x1.0080a711c0000p-8},
        {0x1.fc07ef6b6e30bp-1, 0x1.fe03018e80000p-8},
        {0x1.f4465aa1024afp-1, 0x1.7b91986450000p-6},
        {0x1.ecc07a8fd3f5ep-1, 0x1.39e88608c8000p-5},
        {0x1.e573ad856b537p-1, 0x1.b42dc6e624000p-5},
        {0x1.de5d6dc7b8057p-1, 0x1.165372ec20000p-4},
        {0x1.d77b6498bddf7p-1, 0x1.51b07a0170000p-4},
        {0x1.d0cb580315c0fp-1, 0x1.8c3465c7ea000p-4},
        {0x1.ca4b30d1cf449p-1, 0x1.c5e544a290000p-4},
        {0x1.c3f8ef4810d8ep-1, 0x1.fec91aa0a6000p-4},
        {0x1.bdd2b8b311f44p-1, 0x1.1b72acdc5c000p-3},
        {0x1.b7d6c2eeac054p-1, 0x1.371fc65a98000p-3},
        {0x1.b20363474c8f5p-1, 0x1.526e61c1aa000p-3},
        {0x1.ac570165eeab1p-1, 0x1.6d60ffc240000p-3},
        {0x1.a6d019f331df4p-1, 0x1.87fa08a013000p-3},
        {0x1.a16d3ebc9e3c3p-1, 0x1.a23bc630c3000p-3},
        {0x1.9c2d14567ef45p-1, 0x1.bc286a3512000p-3},
        {0x1.970e4efae9169p-1, 0x1.d5c2195697000p-3},
        {0x1.920fb3bd0b802p-1, 0x1.ef0ae132d3000p-3},
        {0x1.8d3018b58699ap-1, 0x1.040259974e000p-2},
        {0x1.886e5ff170ee6p-1, 0x1.1058bd40e2000p-2},
        {0x1.83c977ad35d27p-1, 0x1.1c898c1137800p-2},
        {0x1.7f405ed16c520p-1, 0x1.2895a3e65b000p-2},
        {0x1.7ad220d0335c4p-1, 0x1.347dd8f6bd000p-2},
        {0x1.767dce53474fdp-1, 0x1.4043083cb3800p-2},
#elif LOG64_N == 128
        {0x1.734f0c3e0de9fp+0, -0x1.7cc7f79e69000p-2},
        {0x1.713786a2ce91fp+0, -0x1.76feec20d0000p-2},
        {0x1.6f26008fab5a0p+0, -0x1.713e31351e000p-2},
        {0x1.6d1a61f138c7dp+0, -0x1.6b85b38287800p-2},
        {0x1.6b1490bc5b4d1p+0, -0x1.65d5590807800p-2},
        {0x1.69147332f0cbap+0, -0x1.602d076180000p-2},
        {0x1.6719f18224223p+0, -0x1.5a8ca86909000p-2},
        {0x1.6524f99a51ed9p+0, -0x1.54f4356035000p-2},
        {0x1.63356aa8f24c4p+0, -0x1.4f637c36b4000p-2},
        {0x1.614b36b9ddc14p+0, -0x1.49da7fda85000p-2},
        {0x1.5f66452c65c4cp+0, -0x1.445923989a800p-2},
        {0x1.5d867b5912c4fp+0, -0x1.3edf439b0b800p-2},
        {0x1.5babccb5b90dep+0, -0x1.396ce448f7000p-2},
        {0x1.59d61f2d91a78p+0, -0x1.3401e17bda000p-2},
        {0x1.5805612465687p+0, -0x1.2e9e2ef468000p-2},
        {0x1.56397cee76bd3p+0, -0x1.2941b3830e000p-2},
        {0x1.54725e2a77f93p+0, -0x1.23ec58cda8800p-2},
        {0x1.52aff42064583p+0, -0x1.1e9e129279000p-2},
        {0x1.50f22dbb2bddfp+0, -0x1.1956d2b48f800p-2},
        {0x1.4f38f4734ded7p+0, -0x1.141679ab9f800p-2},
        {0x1.4d843cfde2840p+0, -0x1.0edd094ef9800p-2},
        {0x1.4bd3ec078a3c8p+0, -0x1.09aa518db1000p-2},
        {0x1.4a27fc3e0258ap+0, -0x1.047e65263b800p-2},
        {0x1.4880524d48434p+0, -0x1.feb224586f000p-3},
        {0x1.46dce1b192d0bp+0, -0x1.f474a7517b000p-3},
        {0x1.453d9d3391854p+0, -0x1.ea4443d103000p-3},
        {0x1.43a2744b4845ap+0, -0x1.e020d44e9b000p-3},
        {0x1.420b54115f8fbp+0, -0x1.d60a22977f000p-3},
        {0x1.40782da3ef4b1p+0, -0x1.cc00104959000p-3},
        {0x1.3ee8f5d57fe8fp+0, -0x1.c202956891000p-3},
        {0x1.3d5d9a00b4ce9p+0, -0x1.b81178d811000p-3},
        {0x1.3bd60c010c12bp+0, -0x1.ae2c9ccd3d000p-3},
        {0x1.3a5242b75dab8p+0, -0x1.a45402e129000p-3},
        {0x1.38d22cd9fd002p+0, -0x1.9a877681df000p-3},
        {0x1.3755bc5847a1cp+0, -0x1.90c6d69483000p-3},
        {0x1.35dce49ad36e2p+0, -0x1.87120a645c000p-3},
        {0x1.34679984dd440p+0, -0x1.7d68fb4143000p-3},
        {0x1.32f5cceffcb24p+0, -0x1.73cb83c627000p-3},
        {0x1.3187775a10d49p+0, -0x1.6a39a9b376000p-3},
        {0x1.301c8373e3990p+0, -0x1.60b3154b7a000p-3},
        {0x1.2eb4ebb95f841p+0, -0x1.5737d76243000p-3},
        {0x1.2d50a0219a9d1p+0, -0x1.4dc7b8fc23000p-3},
        {0x1.2bef9a8b7fd2ap+0, -0x1.4462c51d20000p-3},
        {0x1.2a91c7a0c1babp+0, -0x1.3b08abc830000p-3},
        {0x1.293726014b530p+0, -0x1.31b996b490000p-3},
        {0x1.27dfa5757a1f5p+0, -0x1.2875490a44000p-3},
        {0x1.268b39b1d3bbfp+0, -0x1.1f3b9f879a000p-3},
        {0x1.2539d838ff5bdp+0, -0x1.160c8252ca000p-3},
        {0x1.23eb7aac9083bp+0, -0x1.0ce7f57f72000p-3},
        {0x1.22a012ba940b6p+0, -0x1.03cdc49fea000p-3},
        {0x1.2157996cc4132p+0, -0x1.f57bdbc4b8000p-4},
        {0x1.201201dd2fc9bp+0, -0x1.e370896404000p-4},
        {0x1.1ecf4494d480bp+0, -0x1.d17983ef94000p-4},
        {0x1.1d8f5528f6569p+0, -0x1.bf9674ed8a000p-4},
        {0x1.1c52311577e7cp+0, -0x1.adc79202f6000p-4},
        {0x1.1b17c74cb26e9p+0, -0x1.9c0c3e7288000p-4},
        {0x1.19e010c2c1ab6p+0, -0x1.8a646b372c000p-4},
        {0x1.18ab07bb670bdp+0, -0x1.78d01b3ac0000p-4},
        {0x1.1778a25efbcb6p+0, -0x1.674f145380000p-4},
        {0x1.1648d354c31dap+0, -0x1.55e0e6d878000p-4},
        {0x1.151b990275fddp+0, -0x1.4485cdea1e000p-4},
        {0x1.13f0ea432d24cp+0, -0x1.333d94d6aa000p-4},
        {0x1.12c8b7210f9dap+0, -0x1.22079f8c56000p-4},
        {0x1.11a3028ecb531p+0, -0x1.10e4698622000p-4},
        {0x1.107fbda8434afp+0, -0x1.ffa6c6ad20000p-5},
        {0x1.0f5ee0f4e6bb3p+0, -0x1.dda8d4a774000p-5},
        {0x1.0e4065d2a9fcep+0, -0x1.bbcece4850000p-5},
        {0x1.0d244632ca521p+0, -0x1.9a1894012c000p-5},
        {0x1.0c0a77ce2981ap+0, -0x1.788583302c000p-5},
        {0x1.0af2f83c636d1p+0, -0x1.5715e67d68000p-5},
        {0x1.09ddb98a01339p+0, -0x1.35c8a49658000p-5},
        {0x1.08cabaf52e7dfp+0, -0x1.149e364154000p-5},
        {0x1.07b9f2f4e28fbp+0, -0x1.e72c082eb8000p-6},
        {0x1.06ab58c358f19p+0, -0x1.a55f152528000p-6},
        {0x1.059eea5ecf92cp+0, -0x1.63d62cf818000p-6},
        {0x1.04949cdd12c90p+0, -0x1.228fb8caa0000p-6},
        {0x1.038c6c6f0ada9p+0, -0x1.c317b20f90000p-7},
        {0x1.02865137932a9p+0, -0x1.419355daa0000p-7},
        {0x1.0182427ea7348p+0, -0x1.81203c2ec0000p-8},
        {0x1.008040614b195p+0, -0x1.0040979240000p-9},
        {0x1.fe01ff726fa1ap-1, 0x1.feff384900000p-9},
        {0x1.fa11cc261ea74p-1, 0x1.7dc41353d0000p-7},
        {0x1.f6310b081992ep-1, 0x1.3cea3c4c28000p-6},
        {0x1.f25f63ceeadcdp-1, 0x1.b9fc114890000p-6},
        {0x1.ee9c8039113e7p-1, 0x1.1b0d8ce110000p-5},
        {0x1.eae8078cbb1abp-1, 0x1.58a5bd001c000p-5},
        {0x1.e741aa29d0c9bp-1, 0x1.95c8340d88000p-5},
        {0x1.e3a91830a99b5p-1, 0x1.d276aef578000p-5},
        {0x1.e01e009609a56p-1, 0x1.07598e598c000p-4},
        {0x1.dca01e577bb98p-1, 0x1.253f5e30d2000p-4},
        {0x1.d92f20b7c9103p-1, 0x1.42edd8b380000p-4},
        {0x1.d5cac66fb5ccep-1, 0x1.606598757c000p-4},
        {0x1.d272caa5ede9dp-1, 0x1.7da76356a0000p-4},
        {0x1.cf26e3e6b2ccdp-1, 0x1.9ab434e1c6000p-4},
        {0x1.cbe6da2a77902p-1, 0x1.b78c7bb0d6000p-4},
        {0x1.c8b266d37086dp-1, 0x1.d431332e72000p-4},
        {0x1.c5894bd5d5804p-1, 0x1.f0a3171de6000p-4},
        {0x1.c26b533bb9f8cp-1, 0x1.067152b914000p-3},
        {0x1.bf583eeece73fp-1, 0x1.147858292b000p-3},
        {0x1.bc4fd75db96c1p-1, 0x1.2266ecdca3000p-3},
        {0x1.b951e0c864a28p-1, 0x1.303d7a6c55000p-3},
        {0x1.b65e2c5ef3e2cp-1, 0x1.3dfc33c331000p-3},
        {0x1.b374867c9888bp-1, 0x1.4ba366b7a8000p-3},
        {0x1.b094b211d304ap-1, 0x1.5933928d1f000p-3},
        {0x1.adbe885f2ef7ep-1, 0x1.66acd2418f000p-3},
        {0x1.aaf1d31603da2p-1, 0x1.740f8ec669000p-3},
        {0x1.a82e63fd358a7p-1, 0x1.815c0f51af000p-3},
        {0x1.a5740ef09738bp-1, 0x1.8e92954f68000p-3},
        {0x1.a2c2a90ab4b27p-1, 0x1.9bb3602f84000p-3},
        {0x1.a01a01393f2d1p-1, 0x1.a8bed1c2c0000p-3},
        {0x1.9d79f24db3c1bp-1, 0x1.b5b515c01d000p-3},
        {0x1.9ae2505c7b190p-1, 0x1.c2967ccbcc000p-3},
        {0x1.9852ef297ce2fp-1, 0x1.cf635d5486000p-3},
        {0x1.95cbaeea44b75p-1, 0x1.dc1bd3446c000p-3},
        {0x1.934c69de74838p-1, 0x1.e8c01b8cfe000p-3},
        {0x1.90d4f2f6752e6p-1, 0x1.f5509c0179000p-3},
        {0x1.8e6528effd79dp-1, 0x1.00e6c121fb800p-2},
        {0x1.8bfce9fcc007cp-1, 0x1.071b80e93d000p-2},
        {0x1.899c0dabec30ep-1, 0x1.0d46b9e867000p-2},
        {0x1.87427aa2317fbp-1, 0x1.13687334bd000p-2},
        {0x1.84f00acb39a08p-1, 0x1.1980d67234800p-2},
        {0x1.82a49e8653e55p-1, 0x1.1f8ffe0cc8000p-2},
        {0x1.8060195f40260p-1, 0x1.2595fd7636800p-2},
        {0x1.7e22563e0a329p-1, 0x1.2b9300914a800p-2},
        {0x1.7beb377dcb5adp-1, 0x1.3187210436000p-2},
        {0x1.79baa679725c2p-1, 0x1.377266dec1800p-2},
        {0x1.77907f2170657p-1, 0x1.3d54ffbaf3000p-2},
        {0x1.756cadbd6130cp-1, 0x1.432eee32fe000p-2},
#endif
    },
#if !MATS_HAVE_FAST_FMA
    .tab2 = {
#if LOG64_N == 64
        {0x1.61ffff94c4fecp-1, -0x1.9fe4fc998f325p-56},
        {0x1.66000020377ddp-1, 0x1.e804c7a9519f2p-55},
        {0x1.6a00004c41678p-1, 0x1.902c675d9ecfep-55},
        {0x1.6dffff7384f87p-1, -0x1.2fd6b95e55043p-56},
        {0x1.720000b37216ep-1, 0x1.802bc8d437043p-55},
        {0x1.75ffffbeb3c9dp-1, 0x1.6047ad0a0d4e4p-57},
        {0x1.7a0000628daep-1, -0x1.e00434b49313dp-56},
        {0x1.7dffffd7abd1ap-1, -0x1.6015f8a083576p-56},
        {0x1.81ffffdf40c54p-1, 0x1.7f54bf76a42c9p-57},
        {0x1.860000f334e11p-1, 0x1.60054cb5344d7p-56},
        {0x1.8a0001238aca7p-1, 0x1.c03c9bd132f55p-57},
        {0x1.8dffffb81d212p-1, -0x1.001e519f2764fp-55},
        {0x1.92000086adc7cp-1, 0x1.1fe40f88f49c6p-55},
        {0x1.960000135d8eap-1, -0x1.f832268dc3095p-55},
        {0x1.99ffff9435acp-1, 0x1.7031d8b835edcp-56},
        {0x1.9e00003478565p-1, -0x1.0030b221ce3eep-58},
        {0x1.a20000b592948p-1, 0x1.8fd2f1dbd4639p-55},
        {0x1.a600000ad0bcfp-1, 0x1.901d6a974e6bep-55},
        {0x1.a9ffff55953a5p-1, 0x1.a07556192db98p-57},
        {0x1.adffff29ce03dp-1, -0x1.fff0717ec71c2p-56},
        {0x1.b1ffff34f3ac8p-1, 0x1.8005573de89d1p-57},
        {0x1.b60000894c55bp-1, -0x1.ff2fb51b044c7p-57},
        {0x1.b9fffef45ec7dp-1, -0x1.9ff7c4e8730fp-56},
        {0x1.be0000cda7b2ap-1, 0x1.57d058dbf3c1dp-55},
        {0x1.c1ffff2c57917p-1, 0x1.7e66d7e48dbc9p-58},
        {0x1.c60000ea5b82ap-1, -0x1.47f5e132ed4bep-55},
        {0x1.ca0001121ae98p-1, -0x1.40958c8d5e00ap-58},
        {0x1.ce0000f9241cbp-1, -0x1.7da063caa81c8p-59},
        {0x1.d1fffe8be95a4p-1, -0x1.82e3a411afcd9p-59},
        {0x1.d5ffff035932bp-1, -0x1.00f901b3fe87dp-58},
        {0x1.d9fffe8b54ba7p-1, 0x1.ffef55d6e3a4p-55},
        {0x1.de0000ad95d19p-1, 0x1.5feb2efd4c7c7p-55},
        {0x1.e1fffe925ce47p-1, 0x1.c8085484eaf08p-55},
        {0x1.e5fffe3ddf853p-1, -0x1.fd5ed02c5cadp-60},
        {0x1.e9fffed0a0e5fp-1, -0x1.a80aaef411586p-55},
        {0x1.ee00008f82eep-1, -0x1.b000aeaf97276p-55},
        {0x1.f20000a22d2f4p-1, -0x1.8f8906e13eba3p-56},
        {0x1.f5fffee35b57dp-1, 0x1.1fdd33b2d3714p-57},
        {0x1.fa00014eec3a6p-1, -0x1.3ee0b7a18c1a5p-58},
        {0x1.fdffff5daa89fp-1, -0x1.c1e24c8e3b503p-58},
        {0x1.0200005b93349p+0, -0x1.50197fe6bedcap-54},
        {0x1.05ffff9d597acp+0, 0x1.20160d062d0dcp-55},
        {0x1.0a00005687a63p+0, -0x1.27f3f9307696ep-54},
        {0x1.0dffff779164ep+0, 0x1.b7eb40bb9c4f4p-54},
        {0x1.12000044a0aa8p+0, 0x1.efbc914d512c4p-55},
        {0x1.16000069685bcp+0, -0x1.c0bea3eb2d82cp-57},
        {0x1.1a000093f0d78p+0, 0x1.1fecbf1e8c52p-54},
        {0x1.1dffffb2b1457p+0, -0x1.3fc91365637d6p-55},
        {0x1.2200008824a1p+0, -0x1.dff7e9feb578ap-54},
        {0x1.25ffffeef953p+0, -0x1.b00a61ec912f7p-55},
        {0x1.2a0000a1e7783p+0, 0x1.60048318b0483p-56},
        {0x1.2e0000853d4c7p+0, -0x1.77fbedf2c8cf3p-54},
        {0x1.320000324c55bp+0, 0x1.f81983997354fp-54},
        {0x1.360000594f796p+0, -0x1.cfe4beff900a9p-54},
        {0x1.3a0000a4c1c0fp+0, 0x1.07dbb2e268d0ep-54},
        {0x1.3e0000751c61bp+0, 0x1.80583ed1c566ep-56},
        {0x1.42000069e8a9fp+0, 0x1.f01f1edf82045p-54},
        {0x1.460000b5a1e34p+0, -0x1.dfdf0cf45c14ap-55},
        {0x1.4a0000187e513p+0, 0x1.401306b83a98dp-55},
        {0x1.4dffff3ba420bp+0, 0x1.9fc6539a6454ep-56},
        {0x1.51fffffe391c9p+0, -0x1.601ef3353ac83p-54},
        {0x1.560000e342455p+0, 0x1.3fb7fac8ac151p-55},
        {0x1.59ffffc39676fp+0, 0x1.4fe7dd6659cc2p-55},
        {0x1.5dfffff10ef42p+0, -0x1.48154cb592bcbp-54},
#elif LOG64_N == 128
        {0x1.61000014fb66bp-1, 0x1.e026c91425b3cp-56},
        {0x1.63000034db495p-1, 0x1.dbfea48005d41p-55},
        {0x1.650000d94d478p-1, 0x1.e7fa786d6a5b7p-55},
        {0x1.67000074e6fadp-1, 0x1.1fcea6b54254cp-57},
        {0x1.68ffffedf0faep-1, -0x1.c7e274c590efdp-56},
        {0x1.6b0000763c5bcp-1, -0x1.ac16848dcda01p-55},
        {0x1.6d0001e5cc1f6p-1, 0x1.33f1c9d499311p-55},
        {0x1.6efffeb05f63ep-1, -0x1.e80041ae22d53p-56},
        {0x1.710000e86978p-1, 0x1.bff6671097952p-56},
        {0x1.72ffffc67e912p-1, 0x1.c00e226bd8724p-55},
        {0x1.74fffdf81116ap-1, -0x1.e02916ef101d2p-57},
        {0x1.770000f679c9p-1, -0x1.7fc71cd549c74p-57},
        {0x1.78ffffa7ec835p-1, 0x1.1bec19ef50483p-55},
        {0x1.7affffe20c2e6p-1, -0x1.07e1729cc6465p-56},
        {0x1.7cfffed3fc9p-1, -0x1.08072087b8b1cp-55},
        {0x1.7efffe9261a76p-1, 0x1.dc0286d9df9aep-55},
        {0x1.81000049ca3e8p-1, 0x1.97fd251e54c33p-55},
        {0x1.8300017932c8fp-1, -0x1.afee9b630f381p-55},
        {0x1.850000633739cp-1, 0x1.9bfbf6b6535bcp-55},
        {0x1.87000204289c6p-1, -0x1.bbf65f3117b75p-55},
        {0x1.88fffebf57904p-1, -0x1.9006ea23dcb57p-55},
        {0x1.8b00022bc04dfp-1, -0x1.d00df38e04b0ap-56},
        {0x1.8cfffe50c1b8ap-1, -0x1.8007146ff9f05p-55},
        {0x1.8effffc918e43p-1, 0x1.3817bd07a7038p-55},
        {0x1.910001efa5fc7p-1, 0x1.93e9176dfb403p-55},
        {0x1.9300013467bb9p-1, 0x1.f804e4b980276p-56},
        {0x1.94fffe6ee076fp-1, -0x1.f7ef0d9ff622ep-55},
        {0x1.96fffde3c12d1p-1, -0x1.082aa962638bap-56},
        {0x1.98ffff4458a0dp-1, -0x1.7801b9164a8efp-55},
        {0x1.9afffdd982e3ep-1, -0x1.740e08a5a9337p-55},
        {0x1.9cfffed49fb66p-1, 0x1.fce08c19bep-60},
        {0x1.9f00020f19c51p-1, -0x1.a3faa27885b0ap-55},
        {0x1.a10001145b006p-1, 0x1.4ff489958da56p-56},
        {0x1.a300007bbf6fap-1, 0x1.cbeab8a2b6d18p-55},
        {0x1.a500010971d79p-1, 0x1.8fecadd78793p-55},
        {0x1.a70001df52e48p-1, -0x1.f41763dd8abdbp-55},
        {0x1.a90001c593352p-1, -0x1.ebf0284c27612p-55},
        {0x1.ab0002a4f3e4bp-1, -0x1.9fd043cff3f5fp-57},
        {0x1.acfffd7ae1ed1p-1, -0x1.23ee7129070b4p-55},
        {0x1.aefffee510478p-1, 0x1.a063ee00edea3p-57},
        {0x1.b0fffdb650d5bp-1, 0x1.a06c8381f0ab9p-58},
        {0x1.b2ffffeaaca57p-1, -0x1.9011e74233c1dp-56},
        {0x1.b4fffd995badcp-1, -0x1.9ff1068862a9fp-56},
        {0x1.b7000249e659cp-1, 0x1.aff45d0864f3ep-55},
        {0x1.b8ffff987164p-1, 0x1.cfe7796c2c3f9p-56},
        {0x1.bafffd204cb4fp-1, -0x1.3ff27eef22bc4p-57},
        {0x1.bcfffd2415c45p-1, -0x1.cffb7ee3bea21p-57},
        {0x1.beffff86309dfp-1, -0x1.14103972e0b5cp-55},
        {0x1.c0fffe1b57653p-1, 0x1.bc16494b76a19p-55},
        {0x1.c2ffff1fa57e3p-1, -0x1.4feef8d30c6edp-57},
        {0x1.c4fffdcbfe424p-1, -0x1.43f68bcec4775p-55},
        {0x1.c6fffed54b9f7p-1, 0x1.47ea3f053e0ecp-55},
        {0x1.c8fffeb998fd5p-1, 0x1.383068df992f1p-56},
        {0x1.cb0002125219ap-1, -0x1.8fd8e64180e04p-57},
        {0x1.ccfffdd94469cp-1, 0x1.e7ebe1cc7ea72p-55},
        {0x1.cefffeafdc476p-1, 0x1.ebe39ad9f88fep-55},
        {0x1.d1000169af82bp-1, 0x1.57d91a8b95a71p-56},
        {0x1.d30000d0ff71dp-1, 0x1.9c1906970c7dap-55},
        {0x1.d4fffea790fc4p-1, -0x1.80e37c558fe0cp-58},
        {0x1.d70002edc87e5p-1, -0x1.f80d64dc10f44p-56},
        {0x1.d900021dc82aap-1, -0x1.47c8f94fd5c5cp-56},
        {0x1.dafffd86b0283p-1, 0x1.c7f1dc521617ep-55},
        {0x1.dd000296c4739p-1, 0x1.8019eb2ffb153p-55},
        {0x1.defffe54490f5p-1, 0x1.e00d2c652cc89p-57},
        {0x1.e0fffcdabf694p-1, -0x1.f8340202d69d2p-56},
        {0x1.e2fffdb52c8ddp-1, 0x1.b00c1ca1b0864p-56},
        {0x1.e4ffff24216efp-1, 0x1.2ffa8b094ab51p-56},
        {0x1.e6fffe88a5e11p-1, -0x1.7f673b1efbe59p-58},
        {0x1.e9000119eff0dp-1, -0x1.4808d5e0bc801p-55},
        {0x1.eafffdfa51744p-1, 0x1.80006d54320b5p-56},
        {0x1.ed0001a127fa1p-1, -0x1.002f860565c92p-58},
        {0x1.ef00007babcc4p-1, -0x1.540445d35e611p-55},
        {0x1.f0ffff57a8d02p-1, -0x1.ffb3139ef9105p-59},
        {0x1.f30001ee58ac7p-1, 0x1.a81acf2731155p-55},
        {0x1.f4ffff5823494p-1, 0x1.a3f41d4d7c743p-55},
        {0x1.f6ffffca94c6bp-1, -0x1.202f41c987875p-57},
        {0x1.f8fffe1f9c441p-1, 0x1.77dd1f477e74bp-56},
        {0x1.fafffd2e0e37ep-1, -0x1.f01199a7ca331p-57},
        {0x1.fd0001c77e49ep-1, 0x1.181ee4bceacb1p-56},
        {0x1.feffff7e0c331p-1, -0x1.e05370170875ap-57},
        {0x1.00ffff465606ep+0, -0x1.a7ead491c0adap-55},
        {0x1.02ffff3867a58p+0, -0x1.77f69c3fcb2ep-54},
        {0x1.04ffffdfc0d17p+0, 0x1.7bffe34cb945bp-54},
        {0x1.0700003cd4d82p+0, 0x1.20083c0e456cbp-55},
        {0x1.08ffff9f2cbe8p+0, -0x1.dffdfbe37751ap-57},
        {0x1.0b000010cda65p+0, -0x1.13f7faee626ebp-54},
        {0x1.0d00001a4d338p+0, 0x1.07dfa79489ff7p-55},
        {0x1.0effffadafdfdp+0, -0x1.7040570d66bcp-56},
        {0x1.110000bbafd96p+0, 0x1.e80d4846d0b62p-55},
        {0x1.12ffffae5f45dp+0, 0x1.dbffa64fd36efp-54},
        {0x1.150000dd59ad9p+0, 0x1.a0077701250aep-54},
        {0x1.170000f21559ap+0, 0x1.dfdf9e2e3deeep-55},
        {0x1.18ffffc275426p+0, 0x1.10030dc3b7273p-54},
        {0x1.1b000123d3c59p+0, 0x1.97f7980030188p-54},
        {0x1.1cffff8299eb7p+0, -0x1.5f932ab9f8c67p-57},
        {0x1.1effff48ad4p+0, 0x1.37fbf9da75bebp-54},
        {0x1.210000c8b86a4p+0, 0x1.f806b91fd5b22p-54},
        {0x1.2300003854303p+0, 0x1.3ffc2eb9fbf33p-54},
        {0x1.24fffffbcf684p+0, 0x1.601e77e2e2e72p-56},
        {0x1.26ffff52921d9p+0, 0x1.ffcbb767f0c61p-56},
        {0x1.2900014933a3cp+0, -0x1.202ca3c02412bp-56},
        {0x1.2b00014556313p+0, -0x1.2808233f21f02p-54},
        {0x1.2cfffebfe523bp+0, -0x1.8ff7e384fdcf2p-55},
        {0x1.2f0000bb8ad96p+0, -0x1.5ff51503041c5p-55},
        {0x1.30ffffb7ae2afp+0, -0x1.10071885e289dp-55},
        {0x1.32ffffeac5f7fp+0, -0x1.1ff5d3fb7b715p-54},
        {0x1.350000ca66756p+0, 0x1.57f82228b82bdp-54},
        {0x1.3700011fbf721p+0, 0x1.000bac40dd5ccp-55},
        {0x1.38ffff9592fb9p+0, -0x1.43f9d2db2a751p-54},
        {0x1.3b00004ddd242p+0, 0x1.57f6b707638e1p-55},
        {0x1.3cffff5b2c957p+0, 0x1.a023a10bf1231p-56},
        {0x1.3efffeab0b418p+0, 0x1.87f6d66b152bp-54},
        {0x1.410001532aff4p+0, 0x1.7f8375f198524p-57},
        {0x1.4300017478b29p+0, 0x1.301e672dc5143p-55},
        {0x1.44fffe795b463p+0, 0x1.9ff69b8b2895ap-55},
        {0x1.46fffe80475ep+0, -0x1.5c0b19bc2f254p-54},
        {0x1.48fffef6fc1e7p+0, 0x1.b4009f23a2a72p-54},
        {0x1.4afffe5bea704p+0, -0x1.4ffb7bf0d7d45p-54},
        {0x1.4d000171027dep+0, -0x1.9c06471dc6a3dp-54},
        {0x1.4f0000ff03ee2p+0, 0x1.77f890b85531cp-54},
        {0x1.5100012dc4bd1p+0, 0x1.004657166a436p-57},
        {0x1.530001605277ap+0, -0x1.6bfcece233209p-54},
        {0x1.54fffecdb704cp+0, -0x1.902720505a1d7p-55},
        {0x1.56fffef5f54a9p+0, 0x1.bbfe60ec96412p-54},
        {0x1.5900017e61012p+0, 0x1.87ec581afef9p-55},
        {0x1.5b00003c93e92p+0, -0x1.f41080abf0ccp-54},
        {0x1.5d0001d4919bcp+0, -0x1.8812afb254729p-54},
        {0x1.5efffe7b87a89p+0, -0x1.47eb780ed6904p-54},
#endif
    },
#endif /* !MATS_HAVE_FAST_FMA */
};
#undef LOG64_N

// TODO(michiel): Compile assert LOG2_64_N == 64
global const Log2Data gLog2Data64 =
{
    .invln2hi = 0x1.7154765200000p+0,
    .invln2lo = 0x1.705fc2eefa200p-33,
    .poly = {
        // relative error: 0x1.a72c2bf8p-58
        // abs error: 0x1.67a552c8p-66
        // in -0x1.f45p-8 0x1.f45p-8
        -0x1.71547652b8339p-1,
        0x1.ec709dc3a04bep-2,
        -0x1.7154764702ffbp-2,
        0x1.2776c50034c48p-2,
        -0x1.ec7b328ea92bcp-3,
        0x1.a6225e117f92ep-3,
    },
    .poly1 = {
        // relative error: 0x1.2fad8188p-63
        // in -0x1.5b51p-5 0x1.6ab2p-5
        -0x1.71547652b82fep-1,
        0x1.ec709dc3a03f7p-2,
        -0x1.71547652b7c3fp-2,
        0x1.2776c50f05be4p-2,
        -0x1.ec709dd768fe5p-3,
        0x1.a61761ec4e736p-3,
        -0x1.7153fbc64a79bp-3,
        0x1.484d154f01b4ap-3,
        -0x1.289e4a72c383cp-3,
        0x1.0b32f285aee66p-3,
    },
    /* Algorithm:

	x = 2^k z
	log2(x) = k + log2(c) + log2(z/c)
	log2(z/c) = poly(z/c - 1)

where z is in [1.6p-1; 1.6p0] which is split into N subintervals and z falls
into the ith one, then table entries are computed as

	tab[i].invc = 1/c
	tab[i].logc = (double)log2(c)
	tab2[i].chi = (double)c
	tab2[i].clo = (double)(c - (double)c)

where c is near the center of the subinterval and is chosen by trying +-2^29
floating point invc candidates around 1/center and selecting one for which

	1) the rounding error in 0x1.8p10 + logc is 0,
	2) the rounding error in z - chi - clo is < 0x1p-64 and
	3) the rounding error in (double)log2(c) is minimized (< 0x1p-68).

Note: 1) ensures that k + logc can be computed without rounding error, 2)
ensures that z/c - 1 can be computed as (z - chi - clo)*invc with close to a
single rounding error when there is no fast fma for z*invc - 1, 3) ensures
that logc + poly(z/c - 1) has small error, however near x == 1 when
|log2(x)| < 0x1p-4, this is not enough so that is special cased.  */
    .tab = {
        {0x1.724286bb1acf8p+0, -0x1.1095feecdb000p-1},
        {0x1.6e1f766d2cca1p+0, -0x1.08494bd76d000p-1},
        {0x1.6a13d0e30d48ap+0, -0x1.00143aee8f800p-1},
        {0x1.661ec32d06c85p+0, -0x1.efec5360b4000p-2},
        {0x1.623fa951198f8p+0, -0x1.dfdd91ab7e000p-2},
        {0x1.5e75ba4cf026cp+0, -0x1.cffae0cc79000p-2},
        {0x1.5ac055a214fb8p+0, -0x1.c043811fda000p-2},
        {0x1.571ed0f166e1ep+0, -0x1.b0b67323ae000p-2},
        {0x1.53909590bf835p+0, -0x1.a152f5a2db000p-2},
        {0x1.5014fed61adddp+0, -0x1.9217f5af86000p-2},
        {0x1.4cab88e487bd0p+0, -0x1.8304db0719000p-2},
        {0x1.49539b4334feep+0, -0x1.74189f9a9e000p-2},
        {0x1.460cbdfafd569p+0, -0x1.6552bb5199000p-2},
        {0x1.42d664ee4b953p+0, -0x1.56b23a29b1000p-2},
        {0x1.3fb01111dd8a6p+0, -0x1.483650f5fa000p-2},
        {0x1.3c995b70c5836p+0, -0x1.39de937f6a000p-2},
        {0x1.3991c4ab6fd4ap+0, -0x1.2baa1538d6000p-2},
        {0x1.3698e0ce099b5p+0, -0x1.1d98340ca4000p-2},
        {0x1.33ae48213e7b2p+0, -0x1.0fa853a40e000p-2},
        {0x1.30d191985bdb1p+0, -0x1.01d9c32e73000p-2},
        {0x1.2e025cab271d7p+0, -0x1.e857da2fa6000p-3},
        {0x1.2b404cf13cd82p+0, -0x1.cd3c8633d8000p-3},
        {0x1.288b02c7ccb50p+0, -0x1.b26034c14a000p-3},
        {0x1.25e2263944de5p+0, -0x1.97c1c2f4fe000p-3},
        {0x1.234563d8615b1p+0, -0x1.7d6023f800000p-3},
        {0x1.20b46e33eaf38p+0, -0x1.633a71a05e000p-3},
        {0x1.1e2eefdcda3ddp+0, -0x1.494f5e9570000p-3},
        {0x1.1bb4a580b3930p+0, -0x1.2f9e424e0a000p-3},
        {0x1.19453847f2200p+0, -0x1.162595afdc000p-3},
        {0x1.16e06c0d5d73cp+0, -0x1.f9c9a75bd8000p-4},
        {0x1.1485f47b7e4c2p+0, -0x1.c7b575bf9c000p-4},
        {0x1.12358ad0085d1p+0, -0x1.960c60ff48000p-4},
        {0x1.0fef00f532227p+0, -0x1.64ce247b60000p-4},
        {0x1.0db2077d03a8fp+0, -0x1.33f78b2014000p-4},
        {0x1.0b7e6d65980d9p+0, -0x1.0387d1a42c000p-4},
        {0x1.0953efe7b408dp+0, -0x1.a6f9208b50000p-5},
        {0x1.07325cac53b83p+0, -0x1.47a954f770000p-5},
        {0x1.05197e40d1b5cp+0, -0x1.d23a8c50c0000p-6},
        {0x1.03091c1208ea2p+0, -0x1.16a2629780000p-6},
        {0x1.0101025b37e21p+0, -0x1.720f8d8e80000p-8},
        {0x1.fc07ef9caa76bp-1, 0x1.6fe53b1500000p-7},
        {0x1.f4465d3f6f184p-1, 0x1.11ccce10f8000p-5},
        {0x1.ecc079f84107fp-1, 0x1.c4dfc8c8b8000p-5},
        {0x1.e573a99975ae8p-1, 0x1.3aa321e574000p-4},
        {0x1.de5d6f0bd3de6p-1, 0x1.918a0d08b8000p-4},
        {0x1.d77b681ff38b3p-1, 0x1.e72e9da044000p-4},
        {0x1.d0cb5724de943p-1, 0x1.1dcd2507f6000p-3},
        {0x1.ca4b2dc0e7563p-1, 0x1.476ab03dea000p-3},
        {0x1.c3f8ee8d6cb51p-1, 0x1.7074377e22000p-3},
        {0x1.bdd2b4f020c4cp-1, 0x1.98ede8ba94000p-3},
        {0x1.b7d6c006015cap-1, 0x1.c0db86ad2e000p-3},
        {0x1.b20366e2e338fp-1, 0x1.e840aafcee000p-3},
        {0x1.ac57026295039p-1, 0x1.0790ab4678000p-2},
        {0x1.a6d01bc2731ddp-1, 0x1.1ac056801c000p-2},
        {0x1.a16d3bc3ff18bp-1, 0x1.2db11d4fee000p-2},
        {0x1.9c2d14967feadp-1, 0x1.406464ec58000p-2},
        {0x1.970e4f47c9902p-1, 0x1.52dbe093af000p-2},
        {0x1.920fb3982bcf2p-1, 0x1.651902050d000p-2},
        {0x1.8d30187f759f1p-1, 0x1.771d2cdeaf000p-2},
        {0x1.886e5ebb9f66dp-1, 0x1.88e9c857d9000p-2},
        {0x1.83c97b658b994p-1, 0x1.9a80155e16000p-2},
        {0x1.7f405ffc61022p-1, 0x1.abe186ed3d000p-2},
        {0x1.7ad22181415cap-1, 0x1.bd0f2aea0e000p-2},
        {0x1.767dcf99eff8cp-1, 0x1.ce0a43dbf4000p-2},
    },
#if !MATS_HAVE_FAST_FMA
    .tab2 = {
        {0x1.6200012b90a8ep-1, 0x1.904ab0644b605p-55},
        {0x1.66000045734a6p-1, 0x1.1ff9bea62f7a9p-57},
        {0x1.69fffc325f2c5p-1, 0x1.27ecfcb3c90bap-55},
        {0x1.6e00038b95a04p-1, 0x1.8ff8856739326p-55},
        {0x1.71fffe09994e3p-1, 0x1.afd40275f82b1p-55},
        {0x1.7600015590e1p-1, -0x1.2fd75b4238341p-56},
        {0x1.7a00012655bd5p-1, 0x1.808e67c242b76p-56},
        {0x1.7e0003259e9a6p-1, -0x1.208e426f622b7p-57},
        {0x1.81fffedb4b2d2p-1, -0x1.402461ea5c92fp-55},
        {0x1.860002dfafcc3p-1, 0x1.df7f4a2f29a1fp-57},
        {0x1.89ffff78c6b5p-1, -0x1.e0453094995fdp-55},
        {0x1.8e00039671566p-1, -0x1.a04f3bec77b45p-55},
        {0x1.91fffe2bf1745p-1, -0x1.7fa34400e203cp-56},
        {0x1.95fffcc5c9fd1p-1, -0x1.6ff8005a0695dp-56},
        {0x1.9a0003bba4767p-1, 0x1.0f8c4c4ec7e03p-56},
        {0x1.9dfffe7b92da5p-1, 0x1.e7fd9478c4602p-55},
        {0x1.a1fffd72efdafp-1, -0x1.a0c554dcdae7ep-57},
        {0x1.a5fffde04ff95p-1, 0x1.67da98ce9b26bp-55},
        {0x1.a9fffca5e8d2bp-1, -0x1.284c9b54c13dep-55},
        {0x1.adfffddad03eap-1, 0x1.812c8ea602e3cp-58},
        {0x1.b1ffff10d3d4dp-1, -0x1.efaddad27789cp-55},
        {0x1.b5fffce21165ap-1, 0x1.3cb1719c61237p-58},
        {0x1.b9fffd950e674p-1, 0x1.3f7d94194cep-56},
        {0x1.be000139ca8afp-1, 0x1.50ac4215d9bcp-56},
        {0x1.c20005b46df99p-1, 0x1.beea653e9c1c9p-57},
        {0x1.c600040b9f7aep-1, -0x1.c079f274a70d6p-56},
        {0x1.ca0006255fd8ap-1, -0x1.a0b4076e84c1fp-56},
        {0x1.cdfffd94c095dp-1, 0x1.8f933f99ab5d7p-55},
        {0x1.d1ffff975d6cfp-1, -0x1.82c08665fe1bep-58},
        {0x1.d5fffa2561c93p-1, -0x1.b04289bd295f3p-56},
        {0x1.d9fff9d228b0cp-1, 0x1.70251340fa236p-55},
        {0x1.de00065bc7e16p-1, -0x1.5011e16a4d80cp-56},
        {0x1.e200002f64791p-1, 0x1.9802f09ef62ep-55},
        {0x1.e600057d7a6d8p-1, -0x1.e0b75580cf7fap-56},
        {0x1.ea00027edc00cp-1, -0x1.c848309459811p-55},
        {0x1.ee0006cf5cb7cp-1, -0x1.f8027951576f4p-55},
        {0x1.f2000782b7dccp-1, -0x1.f81d97274538fp-55},
        {0x1.f6000260c450ap-1, -0x1.071002727ffdcp-59},
        {0x1.f9fffe88cd533p-1, -0x1.81bdce1fda8bp-58},
        {0x1.fdfffd50f8689p-1, 0x1.7f91acb918e6ep-55},
        {0x1.0200004292367p+0, 0x1.b7ff365324681p-54},
        {0x1.05fffe3e3d668p+0, 0x1.6fa08ddae957bp-55},
        {0x1.0a0000a85a757p+0, -0x1.7e2de80d3fb91p-58},
        {0x1.0e0001a5f3fccp+0, -0x1.1823305c5f014p-54},
        {0x1.11ffff8afbaf5p+0, -0x1.bfabb6680bac2p-55},
        {0x1.15fffe54d91adp+0, -0x1.d7f121737e7efp-54},
        {0x1.1a00011ac36e1p+0, 0x1.c000a0516f5ffp-54},
        {0x1.1e00019c84248p+0, -0x1.082fbe4da5dap-54},
        {0x1.220000ffe5e6ep+0, -0x1.8fdd04c9cfb43p-55},
        {0x1.26000269fd891p+0, 0x1.cfe2a7994d182p-55},
        {0x1.2a00029a6e6dap+0, -0x1.00273715e8bc5p-56},
        {0x1.2dfffe0293e39p+0, 0x1.b7c39dab2a6f9p-54},
        {0x1.31ffff7dcf082p+0, 0x1.df1336edc5254p-56},
        {0x1.35ffff05a8b6p+0, -0x1.e03564ccd31ebp-54},
        {0x1.3a0002e0eaeccp+0, 0x1.5f0e74bd3a477p-56},
        {0x1.3e000043bb236p+0, 0x1.c7dcb149d8833p-54},
        {0x1.4200002d187ffp+0, 0x1.e08afcf2d3d28p-56},
        {0x1.460000d387cb1p+0, 0x1.20837856599a6p-55},
        {0x1.4a00004569f89p+0, -0x1.9fa5c904fbcd2p-55},
        {0x1.4e000043543f3p+0, -0x1.81125ed175329p-56},
        {0x1.51fffcc027f0fp+0, 0x1.883d8847754dcp-54},
        {0x1.55ffffd87b36fp+0, -0x1.709e731d02807p-55},
        {0x1.59ffff21df7bap+0, 0x1.7f79f68727b02p-55},
        {0x1.5dfffebfc3481p+0, -0x1.180902e30e93ep-54},
    },
#endif
};

// TODO(michiel): Compile assert (1 << POW64_LOG_TABLE_BITS) == 128
// TODO(michiel): Compile assert POW64_LOG_POLY_ORDER == 8
global const PowLogData gPowLogData64 =
{
    .ln2hi = 0x1.62e42fefa3800p-1,
    .ln2lo = 0x1.ef35793c76730p-45,
    .poly = {
        // relative error: 0x1.11922ap-70
        // in -0x1.6bp-8 0x1.6bp-8
        // Coefficients are scaled to match the scaling during evaluation.
        -0x1p-1,
        0x1.555555555556p-2 * -2,
        -0x1.0000000000006p-2 * -2,
        0x1.999999959554ep-3 * 4,
        -0x1.555555529a47ap-3 * 4,
        0x1.2495b9b4845e9p-3 * -8,
        -0x1.0002b8b263fc3p-3 * -8,
    },
    /* Algorithm:

	x = 2^k z
	log(x) = k ln2 + log(c) + log(z/c)
	log(z/c) = poly(z/c - 1)

where z is in [0x1.69555p-1; 0x1.69555p0] which is split into N subintervals
and z falls into the ith one, then table entries are computed as

	tab[i].invc = 1/c
	tab[i].logc = round(0x1p43*log(c))/0x1p43
	tab[i].logctail = (double)(log(c) - logc)

where c is chosen near the center of the subinterval such that 1/c has only a
few precision bits so z/c - 1 is exactly representible as double:

	1/c = center < 1 ? round(N/center)/N : round(2*N/center)/N/2

Note: |z/c - 1| < 1/N for the chosen c, |log(c) - logc - logctail| < 0x1p-97,
the last few bits of logc are rounded away so k*ln2hi + logc has no rounding
error and the interval for z is selected such that near x == 1, where log(x)
is tiny, large cancellation error is avoided in logc + poly(z/c - 1).  */
    .tab = {
        {0x1.6a00000000000p+0, 0, -0x1.62c82f2b9c800p-2, 0x1.ab42428375680p-48},
        {0x1.6800000000000p+0, 0, -0x1.5d1bdbf580800p-2, -0x1.ca508d8e0f720p-46},
        {0x1.6600000000000p+0, 0, -0x1.5767717455800p-2, -0x1.362a4d5b6506dp-45},
        {0x1.6400000000000p+0, 0, -0x1.51aad872df800p-2, -0x1.684e49eb067d5p-49},
        {0x1.6200000000000p+0, 0, -0x1.4be5f95777800p-2, -0x1.41b6993293ee0p-47},
        {0x1.6000000000000p+0, 0, -0x1.4618bc21c6000p-2, 0x1.3d82f484c84ccp-46},
        {0x1.5e00000000000p+0, 0, -0x1.404308686a800p-2, 0x1.c42f3ed820b3ap-50},
        {0x1.5c00000000000p+0, 0, -0x1.3a64c55694800p-2, 0x1.0b1c686519460p-45},
        {0x1.5a00000000000p+0, 0, -0x1.347dd9a988000p-2, 0x1.5594dd4c58092p-45},
        {0x1.5800000000000p+0, 0, -0x1.2e8e2bae12000p-2, 0x1.67b1e99b72bd8p-45},
        {0x1.5600000000000p+0, 0, -0x1.2895a13de8800p-2, 0x1.5ca14b6cfb03fp-46},
        {0x1.5600000000000p+0, 0, -0x1.2895a13de8800p-2, 0x1.5ca14b6cfb03fp-46},
        {0x1.5400000000000p+0, 0, -0x1.22941fbcf7800p-2, -0x1.65a242853da76p-46},
        {0x1.5200000000000p+0, 0, -0x1.1c898c1699800p-2, -0x1.fafbc68e75404p-46},
        {0x1.5000000000000p+0, 0, -0x1.1675cababa800p-2, 0x1.f1fc63382a8f0p-46},
        {0x1.4e00000000000p+0, 0, -0x1.1058bf9ae4800p-2, -0x1.6a8c4fd055a66p-45},
        {0x1.4c00000000000p+0, 0, -0x1.0a324e2739000p-2, -0x1.c6bee7ef4030ep-47},
        {0x1.4a00000000000p+0, 0, -0x1.0402594b4d000p-2, -0x1.036b89ef42d7fp-48},
        {0x1.4a00000000000p+0, 0, -0x1.0402594b4d000p-2, -0x1.036b89ef42d7fp-48},
        {0x1.4800000000000p+0, 0, -0x1.fb9186d5e4000p-3, 0x1.d572aab993c87p-47},
        {0x1.4600000000000p+0, 0, -0x1.ef0adcbdc6000p-3, 0x1.b26b79c86af24p-45},
        {0x1.4400000000000p+0, 0, -0x1.e27076e2af000p-3, -0x1.72f4f543fff10p-46},
        {0x1.4200000000000p+0, 0, -0x1.d5c216b4fc000p-3, 0x1.1ba91bbca681bp-45},
        {0x1.4000000000000p+0, 0, -0x1.c8ff7c79aa000p-3, 0x1.7794f689f8434p-45},
        {0x1.4000000000000p+0, 0, -0x1.c8ff7c79aa000p-3, 0x1.7794f689f8434p-45},
        {0x1.3e00000000000p+0, 0, -0x1.bc286742d9000p-3, 0x1.94eb0318bb78fp-46},
        {0x1.3c00000000000p+0, 0, -0x1.af3c94e80c000p-3, 0x1.a4e633fcd9066p-52},
        {0x1.3a00000000000p+0, 0, -0x1.a23bc1fe2b000p-3, -0x1.58c64dc46c1eap-45},
        {0x1.3a00000000000p+0, 0, -0x1.a23bc1fe2b000p-3, -0x1.58c64dc46c1eap-45},
        {0x1.3800000000000p+0, 0, -0x1.9525a9cf45000p-3, -0x1.ad1d904c1d4e3p-45},
        {0x1.3600000000000p+0, 0, -0x1.87fa06520d000p-3, 0x1.bbdbf7fdbfa09p-45},
        {0x1.3400000000000p+0, 0, -0x1.7ab890210e000p-3, 0x1.bdb9072534a58p-45},
        {0x1.3400000000000p+0, 0, -0x1.7ab890210e000p-3, 0x1.bdb9072534a58p-45},
        {0x1.3200000000000p+0, 0, -0x1.6d60fe719d000p-3, -0x1.0e46aa3b2e266p-46},
        {0x1.3000000000000p+0, 0, -0x1.5ff3070a79000p-3, -0x1.e9e439f105039p-46},
        {0x1.3000000000000p+0, 0, -0x1.5ff3070a79000p-3, -0x1.e9e439f105039p-46},
        {0x1.2e00000000000p+0, 0, -0x1.526e5e3a1b000p-3, -0x1.0de8b90075b8fp-45},
        {0x1.2c00000000000p+0, 0, -0x1.44d2b6ccb8000p-3, 0x1.70cc16135783cp-46},
        {0x1.2c00000000000p+0, 0, -0x1.44d2b6ccb8000p-3, 0x1.70cc16135783cp-46},
        {0x1.2a00000000000p+0, 0, -0x1.371fc201e9000p-3, 0x1.178864d27543ap-48},
        {0x1.2800000000000p+0, 0, -0x1.29552f81ff000p-3, -0x1.48d301771c408p-45},
        {0x1.2600000000000p+0, 0, -0x1.1b72ad52f6000p-3, -0x1.e80a41811a396p-45},
        {0x1.2600000000000p+0, 0, -0x1.1b72ad52f6000p-3, -0x1.e80a41811a396p-45},
        {0x1.2400000000000p+0, 0, -0x1.0d77e7cd09000p-3, 0x1.a699688e85bf4p-47},
        {0x1.2400000000000p+0, 0, -0x1.0d77e7cd09000p-3, 0x1.a699688e85bf4p-47},
        {0x1.2200000000000p+0, 0, -0x1.fec9131dbe000p-4, -0x1.575545ca333f2p-45},
        {0x1.2000000000000p+0, 0, -0x1.e27076e2b0000p-4, 0x1.a342c2af0003cp-45},
        {0x1.2000000000000p+0, 0, -0x1.e27076e2b0000p-4, 0x1.a342c2af0003cp-45},
        {0x1.1e00000000000p+0, 0, -0x1.c5e548f5bc000p-4, -0x1.d0c57585fbe06p-46},
        {0x1.1c00000000000p+0, 0, -0x1.a926d3a4ae000p-4, 0x1.53935e85baac8p-45},
        {0x1.1c00000000000p+0, 0, -0x1.a926d3a4ae000p-4, 0x1.53935e85baac8p-45},
        {0x1.1a00000000000p+0, 0, -0x1.8c345d631a000p-4, 0x1.37c294d2f5668p-46},
        {0x1.1a00000000000p+0, 0, -0x1.8c345d631a000p-4, 0x1.37c294d2f5668p-46},
        {0x1.1800000000000p+0, 0, -0x1.6f0d28ae56000p-4, -0x1.69737c93373dap-45},
        {0x1.1600000000000p+0, 0, -0x1.51b073f062000p-4, 0x1.f025b61c65e57p-46},
        {0x1.1600000000000p+0, 0, -0x1.51b073f062000p-4, 0x1.f025b61c65e57p-46},
        {0x1.1400000000000p+0, 0, -0x1.341d7961be000p-4, 0x1.c5edaccf913dfp-45},
        {0x1.1400000000000p+0, 0, -0x1.341d7961be000p-4, 0x1.c5edaccf913dfp-45},
        {0x1.1200000000000p+0, 0, -0x1.16536eea38000p-4, 0x1.47c5e768fa309p-46},
        {0x1.1000000000000p+0, 0, -0x1.f0a30c0118000p-5, 0x1.d599e83368e91p-45},
        {0x1.1000000000000p+0, 0, -0x1.f0a30c0118000p-5, 0x1.d599e83368e91p-45},
        {0x1.0e00000000000p+0, 0, -0x1.b42dd71198000p-5, 0x1.c827ae5d6704cp-46},
        {0x1.0e00000000000p+0, 0, -0x1.b42dd71198000p-5, 0x1.c827ae5d6704cp-46},
        {0x1.0c00000000000p+0, 0, -0x1.77458f632c000p-5, -0x1.cfc4634f2a1eep-45},
        {0x1.0c00000000000p+0, 0, -0x1.77458f632c000p-5, -0x1.cfc4634f2a1eep-45},
        {0x1.0a00000000000p+0, 0, -0x1.39e87b9fec000p-5, 0x1.502b7f526feaap-48},
        {0x1.0a00000000000p+0, 0, -0x1.39e87b9fec000p-5, 0x1.502b7f526feaap-48},
        {0x1.0800000000000p+0, 0, -0x1.f829b0e780000p-6, -0x1.980267c7e09e4p-45},
        {0x1.0800000000000p+0, 0, -0x1.f829b0e780000p-6, -0x1.980267c7e09e4p-45},
        {0x1.0600000000000p+0, 0, -0x1.7b91b07d58000p-6, -0x1.88d5493faa639p-45},
        {0x1.0400000000000p+0, 0, -0x1.fc0a8b0fc0000p-7, -0x1.f1e7cf6d3a69cp-50},
        {0x1.0400000000000p+0, 0, -0x1.fc0a8b0fc0000p-7, -0x1.f1e7cf6d3a69cp-50},
        {0x1.0200000000000p+0, 0, -0x1.fe02a6b100000p-8, -0x1.9e23f0dda40e4p-46},
        {0x1.0200000000000p+0, 0, -0x1.fe02a6b100000p-8, -0x1.9e23f0dda40e4p-46},
        {0x1.0000000000000p+0, 0, 0x0.0000000000000p+0, 0x0.0000000000000p+0},
        {0x1.0000000000000p+0, 0, 0x0.0000000000000p+0, 0x0.0000000000000p+0},
        {0x1.fc00000000000p-1, 0, 0x1.0101575890000p-7, -0x1.0c76b999d2be8p-46},
        {0x1.f800000000000p-1, 0, 0x1.0205658938000p-6, -0x1.3dc5b06e2f7d2p-45},
        {0x1.f400000000000p-1, 0, 0x1.8492528c90000p-6, -0x1.aa0ba325a0c34p-45},
        {0x1.f000000000000p-1, 0, 0x1.0415d89e74000p-5, 0x1.111c05cf1d753p-47},
        {0x1.ec00000000000p-1, 0, 0x1.466aed42e0000p-5, -0x1.c167375bdfd28p-45},
        {0x1.e800000000000p-1, 0, 0x1.894aa149fc000p-5, -0x1.97995d05a267dp-46},
        {0x1.e400000000000p-1, 0, 0x1.ccb73cdddc000p-5, -0x1.a68f247d82807p-46},
        {0x1.e200000000000p-1, 0, 0x1.eea31c006c000p-5, -0x1.e113e4fc93b7bp-47},
        {0x1.de00000000000p-1, 0, 0x1.1973bd1466000p-4, -0x1.5325d560d9e9bp-45},
        {0x1.da00000000000p-1, 0, 0x1.3bdf5a7d1e000p-4, 0x1.cc85ea5db4ed7p-45},
        {0x1.d600000000000p-1, 0, 0x1.5e95a4d97a000p-4, -0x1.c69063c5d1d1ep-45},
        {0x1.d400000000000p-1, 0, 0x1.700d30aeac000p-4, 0x1.c1e8da99ded32p-49},
        {0x1.d000000000000p-1, 0, 0x1.9335e5d594000p-4, 0x1.3115c3abd47dap-45},
        {0x1.cc00000000000p-1, 0, 0x1.b6ac88dad6000p-4, -0x1.390802bf768e5p-46},
        {0x1.ca00000000000p-1, 0, 0x1.c885801bc4000p-4, 0x1.646d1c65aacd3p-45},
        {0x1.c600000000000p-1, 0, 0x1.ec739830a2000p-4, -0x1.dc068afe645e0p-45},
        {0x1.c400000000000p-1, 0, 0x1.fe89139dbe000p-4, -0x1.534d64fa10afdp-45},
        {0x1.c000000000000p-1, 0, 0x1.1178e8227e000p-3, 0x1.1ef78ce2d07f2p-45},
        {0x1.be00000000000p-1, 0, 0x1.1aa2b7e23f000p-3, 0x1.ca78e44389934p-45},
        {0x1.ba00000000000p-1, 0, 0x1.2d1610c868000p-3, 0x1.39d6ccb81b4a1p-47},
        {0x1.b800000000000p-1, 0, 0x1.365fcb0159000p-3, 0x1.62fa8234b7289p-51},
        {0x1.b400000000000p-1, 0, 0x1.4913d8333b000p-3, 0x1.5837954fdb678p-45},
        {0x1.b200000000000p-1, 0, 0x1.527e5e4a1b000p-3, 0x1.633e8e5697dc7p-45},
        {0x1.ae00000000000p-1, 0, 0x1.6574ebe8c1000p-3, 0x1.9cf8b2c3c2e78p-46},
        {0x1.ac00000000000p-1, 0, 0x1.6f0128b757000p-3, -0x1.5118de59c21e1p-45},
        {0x1.aa00000000000p-1, 0, 0x1.7898d85445000p-3, -0x1.c661070914305p-46},
        {0x1.a600000000000p-1, 0, 0x1.8beafeb390000p-3, -0x1.73d54aae92cd1p-47},
        {0x1.a400000000000p-1, 0, 0x1.95a5adcf70000p-3, 0x1.7f22858a0ff6fp-47},
        {0x1.a000000000000p-1, 0, 0x1.a93ed3c8ae000p-3, -0x1.8724350562169p-45},
        {0x1.9e00000000000p-1, 0, 0x1.b31d8575bd000p-3, -0x1.c358d4eace1aap-47},
        {0x1.9c00000000000p-1, 0, 0x1.bd087383be000p-3, -0x1.d4bc4595412b6p-45},
        {0x1.9a00000000000p-1, 0, 0x1.c6ffbc6f01000p-3, -0x1.1ec72c5962bd2p-48},
        {0x1.9600000000000p-1, 0, 0x1.db13db0d49000p-3, -0x1.aff2af715b035p-45},
        {0x1.9400000000000p-1, 0, 0x1.e530effe71000p-3, 0x1.212276041f430p-51},
        {0x1.9200000000000p-1, 0, 0x1.ef5ade4dd0000p-3, -0x1.a211565bb8e11p-51},
        {0x1.9000000000000p-1, 0, 0x1.f991c6cb3b000p-3, 0x1.bcbecca0cdf30p-46},
        {0x1.8c00000000000p-1, 0, 0x1.07138604d5800p-2, 0x1.89cdb16ed4e91p-48},
        {0x1.8a00000000000p-1, 0, 0x1.0c42d67616000p-2, 0x1.7188b163ceae9p-45},
        {0x1.8800000000000p-1, 0, 0x1.1178e8227e800p-2, -0x1.c210e63a5f01cp-45},
        {0x1.8600000000000p-1, 0, 0x1.16b5ccbacf800p-2, 0x1.b9acdf7a51681p-45},
        {0x1.8400000000000p-1, 0, 0x1.1bf99635a6800p-2, 0x1.ca6ed5147bdb7p-45},
        {0x1.8200000000000p-1, 0, 0x1.214456d0eb800p-2, 0x1.a87deba46baeap-47},
        {0x1.7e00000000000p-1, 0, 0x1.2bef07cdc9000p-2, 0x1.a9cfa4a5004f4p-45},
        {0x1.7c00000000000p-1, 0, 0x1.314f1e1d36000p-2, -0x1.8e27ad3213cb8p-45},
        {0x1.7a00000000000p-1, 0, 0x1.36b6776be1000p-2, 0x1.16ecdb0f177c8p-46},
        {0x1.7800000000000p-1, 0, 0x1.3c25277333000p-2, 0x1.83b54b606bd5cp-46},
        {0x1.7600000000000p-1, 0, 0x1.419b423d5e800p-2, 0x1.8e436ec90e09dp-47},
        {0x1.7400000000000p-1, 0, 0x1.4718dc271c800p-2, -0x1.f27ce0967d675p-45},
        {0x1.7200000000000p-1, 0, 0x1.4c9e09e173000p-2, -0x1.e20891b0ad8a4p-45},
        {0x1.7000000000000p-1, 0, 0x1.522ae0738a000p-2, 0x1.ebe708164c759p-45},
        {0x1.6e00000000000p-1, 0, 0x1.57bf753c8d000p-2, 0x1.fadedee5d40efp-46},
        {0x1.6c00000000000p-1, 0, 0x1.5d5bddf596000p-2, -0x1.a0b2a08a465dcp-47},
    },
};
