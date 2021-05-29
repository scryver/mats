
global const f32 gZeroF32       =  0.0f;
global const f32 gOneF32        =  1.0000000000e+00f;  /* 0x3f800000 */
global const f32 gHalfF32       =  5.0000000000e-01f;  /* 0x3f000000 */
global const f32 g2pow8F32      =  2.5600000000e+02f;  /* 0x43800000 = 2^8 */
global const f32 g2pow25F32     =  3.355443200e+07f;   /* 0x4c000000 = 2^25 */
global const f32 g2powMin8F32   =  3.9062500000e-03f;  /* 0x3b800000 = 2^-8 */
global const f32 g2powMin25F32  =  2.9802322388e-08f;  /* 0x33000000 = 2^-25 */
global const f32 g2powMin100F32 =  7.8886090522e-31f;  /* 0x0d800000 = 2^-100 */
global const f32 gLn2F32        =  6.9314718246e-01f;  /* 0x3f317218 = ln(2) */
global const f32 gInvLn2F32     =  1.4426950216e+00f;  /* 0x3fb8aa3b = 1/ln(2) */
global const f32 gInvLn10F32    =  4.3429449201e-01f;  /* 0x3ede5bd9 = 1/ln(10) */
global const f32 gHugeF32       =  1.0e+30f;
global const f32 gTinyF32       =  1.0e-30f;

global const f64 gZeroF64       =  0.0;
global const f64 gOneF64        =  1.0;
global const f64 gHalfF64       =  0.5;
//global const f64 g2pow8F64      =  2.5600000000e+02;  /* 0x43800000 = 2^8 */
//global const f64 g2pow25F64     =  3.355446400e+07;   /* 0x4c000000 = 2^25 */
//global const f64 g2powMin8F64   =  3.9062500000e-03;  /* 0x3b800000 = 2^-8 */
//global const f64 g2powMin25F64  =  2.9802642388e-08;  /* 0x33000000 = 2^-25 */
//global const f64 g2powMin100F64 =  7.8886090522e-31;  /* 0x0d800000 = 2^-100 */
//global const f64 gLn2F64        =  6.9314718246e-01;  /* 0x3f317218 = ln(2) */
//global const f64 gInvLn2F64     =  1.4426950216e+00;  /* 0x3fb8aa3b = 1/ln(2) */
//global const f64 gInvLn10F64    =  4.3429449201e-01;  /* 0x3ede5bd9 = 1/ln(10) */
global const f64 gHugeF64       =  1.0e+300;
//global const f64 gTinyF64       =  1.0e-300;

global const f32 gPiF32         =  3.1415927410e+00f; /* 0x40490fdb, 24 bits of pi */
global const f32 gPiOver2F32    =  1.5707963705e+00f; /* 0x3fc90fdb, 24 bits of pi / 2 */
global const f32 gPiOver4F32    =  7.8539818525e-01f; /* 0x3f490fdb, 24 bits of pi / 4 */
global const f32 g2OverPiF32    =  6.3661980629e-01f; /* 0x3f22f984, 24 bits of 2 / pi */
global const f32 gPiF32_lo      = -8.7422776573e-08f; /* 0x40490fdb, pi - gPiF32 */
global const f32 gPiOver4F32_hi =  7.8539812565e-01f; /* 0x3F490FDA, 24 bits of pi / 4 */
global const f32 gPiOver4F32_lo =  3.7748947079e-08f; /* 0x33222168, pi / 4 - gPiOver4F32_hi */
global const f32 gPiOver2F32_1  =  1.5707855225e+00f; /* 0x3fc90f80, first 17 bit of pi/2 */
global const f32 gPiOver2F32_1t =  1.0804334124e-05f; /* 0x37354443, pi/2 - pio2_1 */
global const f32 gPiOver2F32_2  =  1.0804273188e-05f; /* 0x37354400, second 17 bit of pi/2 */
global const f32 gPiOver2F32_2t =  6.0770999344e-11f; /* 0x2e85a308, pi/2 - (pio2_1+pio2_2) */
global const f32 gPiOver2F32_3  =  6.0770943833e-11f; /* 0x2e85a300, third  17 bit of pi/2 */
global const f32 gPiOver2F32_3t =  6.1232342629e-17f; /* 0x248d3132, pi/2 - (pio2_1+pio2_2+pio2_3) */
global const f32 gLog10_2_hi    =  3.0102920532e-01f; /* 0x3e9a2080, log10(2) */
global const f32 gLog10_2_lo    =  7.9034151668e-07f; /* 0x355427db, log10(2) - gLog10_2_hi */

global const f32 gHalfSignF32s[2]  = {0.5f, -0.5f};
global const f32 gLn2HighF32s[2]   = {6.9313812256e-01, -6.9313812256e-01}; /* 0x3f317180, 0xbf317180 */
global const f32 gLn2LowF32s[2]    = {9.0580006145e-06, -9.0580006145e-06}; /* 0x3717f7d1, 0xb717f7d1 */

// NOTE(michiel): Exp
#define EXP2F_TABLE_BITS 5
#define EXP2F_POLY_ORDER 3
#define EXP2F_N  ((f64)(1 << EXP2F_TABLE_BITS))
global const f64 gExp2F32_Shift = 0x1.8p+52;
global const f64 gExp2F32_ShiftScaled = 0x1.8p+52 / EXP2F_N;
global const f64 gExp2F32_InvLn2Scaled = 0x1.71547652b82fep+0 * EXP2F_N;
global const f64 gExp2F32_Poly[EXP2F_POLY_ORDER] = { 0x1.c6af84b912394p-5, 0x1.ebfce50fac4f3p-3, 0x1.62e42ff0c52d6p-1 };
global const f64 gExp2F32_PolyScaled[EXP2F_POLY_ORDER] = {
    0x1.c6af84b912394p-5 / EXP2F_N / EXP2F_N / EXP2F_N,
    0x1.ebfce50fac4f3p-3 / EXP2F_N / EXP2F_N,
    0x1.62e42ff0c52d6p-1 / EXP2F_N,
};
global const u64 gExp2F32_Table[1 << EXP2F_TABLE_BITS] = {
    0x3ff0000000000000ULL, 0x3fefd9b0d3158574ULL, 0x3fefb5586cf9890fULL, 0x3fef9301d0125b51ULL,
    0x3fef72b83c7d517bULL, 0x3fef54873168b9aaULL, 0x3fef387a6e756238ULL, 0x3fef1e9df51fdee1ULL,
    0x3fef06fe0a31b715ULL, 0x3feef1a7373aa9cbULL, 0x3feedea64c123422ULL, 0x3feece086061892dULL,
    0x3feebfdad5362a27ULL, 0x3feeb42b569d4f82ULL, 0x3feeab07dd485429ULL, 0x3feea47eb03a5585ULL,
    0x3feea09e667f3bcdULL, 0x3fee9f75e8ec5f74ULL, 0x3feea11473eb0187ULL, 0x3feea589994cce13ULL,
    0x3feeace5422aa0dbULL, 0x3feeb737b0cdc5e5ULL, 0x3feec49182a3f090ULL, 0x3feed503b23e255dULL,
    0x3feee89f995ad3adULL, 0x3feeff76f2fb5e47ULL, 0x3fef199bdd85529cULL, 0x3fef3720dcef9069ULL,
    0x3fef5818dcfba487ULL, 0x3fef7c97337b9b5fULL, 0x3fefa4afa2a490daULL, 0x3fefd0765b6e4540ULL,
};
#undef EXP2F_N

// NOTE(michiel): Log
#define LOGF_TABLE_BITS 4
#define LOGF_POLY_ORDER 4
struct LogTable
{
    f64 invC;
    f64 logC;
};
global const f64 gLogF32_Ln2 = 0x1.62e42fefa39efp-1;
global const f64 gLogF32_Poly[LOGF_POLY_ORDER - 1] = {
    -0x1.00ea348b88334p-2,
    0x1.5575b0be00b6ap-2,
    -0x1.ffffef20a4123p-2,
};
global const LogTable gLogF32_Table[1 << LOGF_TABLE_BITS] = {
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
#define LOG2F_TABLE_BITS 4
#define LOG2F_POLY_ORDER 4
global const f64 gLog2F32_Poly[LOG2F_POLY_ORDER] = {
    -0x1.712b6f70a7e4dp-2,
    0x1.ecabf496832ep-2,
    -0x1.715479ffae3dep-1,
    0x1.715475f35c8b8p0,
};
global const LogTable gLog2F32_Table[1 << LOG2F_TABLE_BITS] = {
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
#define POWF_LOG2_TABLE_BITS 4
#define POWF_LOG2_POLY_ORDER 5
#define POWF_SCALE ((f64)(1 << EXP2F_TABLE_BITS))
global const f64 gPowF32_Log2PolyScaled[POWF_LOG2_POLY_ORDER] = {
    0x1.27616c9496e0bp-2 * POWF_SCALE,
    -0x1.71969a075c67ap-2 * POWF_SCALE,
    0x1.ec70a6ca7baddp-2 * POWF_SCALE,
    -0x1.7154748bef6c8p-1 * POWF_SCALE,
    0x1.71547652ab82bp0 * POWF_SCALE,
}; // __powf_log2_data.poly
global const LogTable gPowF32_Log2TableScaled[1 << POWF_LOG2_TABLE_BITS] = {
    { 0x1.661ec79f8f3bep+0, -0x1.efec65b963019p-2 * POWF_SCALE },
    { 0x1.571ed4aaf883dp+0, -0x1.b0b6832d4fca4p-2 * POWF_SCALE },
    { 0x1.49539f0f010bp+0, -0x1.7418b0a1fb77bp-2 * POWF_SCALE },
    { 0x1.3c995b0b80385p+0, -0x1.39de91a6dcf7bp-2 * POWF_SCALE },
    { 0x1.30d190c8864a5p+0, -0x1.01d9bf3f2b631p-2 * POWF_SCALE },
    { 0x1.25e227b0b8eap+0, -0x1.97c1d1b3b7afp-3 * POWF_SCALE },
    { 0x1.1bb4a4a1a343fp+0, -0x1.2f9e393af3c9fp-3 * POWF_SCALE },
    { 0x1.12358f08ae5bap+0, -0x1.960cbbf788d5cp-4 * POWF_SCALE },
    { 0x1.0953f419900a7p+0, -0x1.a6f9db6475fcep-5 * POWF_SCALE },
    { 0x1p+0, 0x0p+0 * POWF_SCALE },
    { 0x1.e608cfd9a47acp-1, 0x1.338ca9f24f53dp-4 * POWF_SCALE },
    { 0x1.ca4b31f026aap-1, 0x1.476a9543891bap-3 * POWF_SCALE },
    { 0x1.b2036576afce6p-1, 0x1.e840b4ac4e4d2p-3 * POWF_SCALE },
    { 0x1.9c2d163a1aa2dp-1, 0x1.40645f0c6651cp-2 * POWF_SCALE },
    { 0x1.886e6037841edp-1, 0x1.88e9c2c1b9ff8p-2 * POWF_SCALE },
    { 0x1.767dcf5534862p-1, 0x1.ce0a44eb17bccp-2 * POWF_SCALE },
}; // __powf_log2_data.tab

global const f64 gPowF32_Log2Poly[POWF_LOG2_POLY_ORDER] = {
    0x1.27616c9496e0bp-2,
    -0x1.71969a075c67ap-2,
    0x1.ec70a6ca7baddp-2,
    -0x1.7154748bef6c8p-1,
    0x1.71547652ab82bp0,
};
global const LogTable gPowF32_Log2Table[1 << POWF_LOG2_TABLE_BITS] = {
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
