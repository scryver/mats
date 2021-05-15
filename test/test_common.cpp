#define MATH_FUNC_F32_FROM_F32(name)                f32 name(f32 x)
typedef MATH_FUNC_F32_FROM_F32(MathFuncF32FromF32);
#define MATH_FUNC_F32_FROM_F32_F32(name)            f32 name(f32 a, f32 b)
typedef MATH_FUNC_F32_FROM_F32_F32(MathFuncF32FromF32F32);
#define MATH_FUNC_F64_FROM_F64(name)                f64 name(f64 x)
typedef MATH_FUNC_F64_FROM_F64(MathFuncF64FromF64);
#define MATH_FUNC_F64_FROM_F64_F64(name)            f64 name(f64 a, f64 b)
typedef MATH_FUNC_F64_FROM_F64_F64(MathFuncF64FromF64F64);
#define MATH_FUNC_F32_4x_FROM_F32_4x(name)          f32_4x name(f32_4x x)
typedef MATH_FUNC_F32_4x_FROM_F32_4x(MathFuncF32_4xFromF32_4x);
#define MATH_FUNC_F32_4x_FROM_F32_4x_F32_4x(name)   f32_4x name(f32_4x a, f32_4x b)
typedef MATH_FUNC_F32_4x_FROM_F32_4x_F32_4x(MathFuncF32_4xFromF32_4xF32_4x);

#define call_comp_x(lib, name, suf, sec)  \
run_comp_f32_x(string(#lib " " #name), 20, #name, tests, minVal, maxVal, name, name ## suf, sec)
#define call_spd(lib, name, suf, sec) \
run_speed_f32(string(#lib " " #name), 20, #name, tests, minVal, maxVal, name ## suf, sec)

#define call_comp_x2(lib, name, suf, sec)  \
run_comp_f32_f32_x(string(#lib " " #name), 20, #name, testsA, minValA, maxValA, testsB, minValB, maxValB, name, name ## suf, sec)
#define call_spd2(lib, name, suf, sec) \
run_speed_f32_f32(string(#lib " " #name), 20, #name, testsA, minValA, maxValA, testsB, minValB, maxValB, name ## suf, sec)

#define call_comp_x_4x(lib, name, suf, sec)  \
run_comp_f32_4x_x(string(#lib " " #name), 20, #name "_4x", tests, minVal, maxVal, name, name ## suf, sec)
#define call_spd_4x(lib, name, suf, sec) \
run_speed_f32_4x(string(#lib " " #name), 20, #name "_4x", tests, minVal, maxVal, name ## suf, sec)

#define call_comp_x2_4x(lib, name, suf, sec)  \
run_comp_f32_f32_4x_x(string(#lib " " #name), 20, #name "_4x", testsA, minValA, maxValA, testsB, minValB, maxValB, name, name ## suf, sec)
#define call_spd2_4x(lib, name, suf, sec) \
run_speed_f32_f32_4x(string(#lib " " #name), 20, #name "_4x", testsA, minValA, maxValA, testsB, minValB, maxValB, name ## suf, sec)

internal void
set_default_fp_behavior(void)
{
#define FLUSH_TO_ZERO_BIT (1 << 15)
#define ROUNDING_CONTROL_BITS (3 << 13)
#define PRECISION_MASK (1 << 12)
#define UNDERFLOW_MASK (1 << 11)
#define OVERFLOW_MASK (1 << 10)
#define DBZ_MASK (1 << 9)
#define DENORMAL_OP_MASK (1 << 8)
#define INVALID_OP_MASK (1 << 7)
#define DENORMALS_ARE_ZERO (1 << 6)
    u32 fpControlMask = (FLUSH_TO_ZERO_BIT |
                         // ROUNDING_CONTROL_BITS |
                         PRECISION_MASK |
                         UNDERFLOW_MASK |
                         OVERFLOW_MASK |
                         DBZ_MASK |
                         DENORMAL_OP_MASK |
                         INVALID_OP_MASK |
                         DENORMALS_ARE_ZERO);
    u32 desiredBits = fpControlMask;
    u32 oldControlBits = _mm_getcsr();
    u32 newControlBits = (oldControlBits & ~fpControlMask) | desiredBits;
    _mm_setcsr(newControlBits);
}

internal f32
absolute(f32 x)
{
    U32F32 xf = u32f32(x);
    xf.u &= ~F32_SIGN_MASK;
    return xf.f;
}

internal f64
absolute(f64 x)
{
    U64F64 xf = u64f64(x);
    xf.u &= ~F64_SIGN_MASK;
    return xf.f;
}

internal f32
square(f32 x)
{
    return x * x;
}

internal u32
string_length(const char *cString)
{
    const char *start = cString;
    while (cString && *cString++);
    return safe_truncate_to_u32(cString - start);
}

internal String
string(char *cStr)
{
    return {string_length(cStr), (u8 *)cStr};
}

internal b32
strings_are_equal(const char *a, const char *b)
{
    u32 lenA = string_length(a);
    u32 lenB = string_length(b);
    b32 result = lenA == lenB;
    if (result && (a != b))
    {
        for (u32 index = 0; index < lenA; ++index) {
            if (*a++ != *b++) {
                result = false;
                break;
            }
        }

    }
    return result;
}

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

struct CompInfo
{
    f64 totalAbsErr;
    f64 totalRelErr;
    f64 minAbsErr;
    f64 maxAbsErr;
    f64 minRelErr;
    f64 maxRelErr;
    u32 numInputs;
    f32 minAbsInputA;
    f32 minAbsInputB;
    f32 maxAbsInputA;
    f32 maxAbsInputB;
    f32 minRelInputA;
    f32 minRelInputB;
    f32 maxRelInputA;
    f32 maxRelInputB;
};

internal void
print_comp_info(String name, u32 maxNameSize, char *func, u32 tests, f32 seconds, CompInfo *info, f32 secondsRatio)
{
    f64 testCalls = (f64)tests;
    f64 oneOverTestCalls = 1.0 / testCalls;
    persist char *spaces = "                                                                ";
    if (info->numInputs == 1)
    {
        fprintf(stdout, "%.*s%.*s: %f, %f sec, %10.0f %s/sec, avg %f usec\n    abs %e (min(%a): %e, max(%a): %e)\n    rel %e (min(%a): %e, max(%a): %e)\n",
                STR_FMT(name), (u32)(maxNameSize - name.size), spaces, secondsRatio,
                seconds, testCalls / seconds, func, seconds / (testCalls / 1000000.0),
                info->totalAbsErr * oneOverTestCalls, info->minAbsInputA, info->minAbsErr, info->maxAbsInputA, info->maxAbsErr,
                info->totalRelErr * oneOverTestCalls, info->minRelInputA, info->minRelErr, info->maxRelInputA, info->maxRelErr);
    }
    else
    {
        i_expect(info->numInputs == 2);

        fprintf(stdout, "%.*s%.*s: %f, %f sec, %10.0f %s/sec, avg %f usec\n    abs %e (min(%a, %a): %e, max(%a, %a): %e)\n    rel %e (min(%a, %a): %e, max(%a, %a): %e)\n",
                STR_FMT(name), (u32)(maxNameSize - name.size), spaces, secondsRatio,
                seconds, testCalls / seconds, func, seconds / (testCalls / 1000000.0),
                info->totalAbsErr * oneOverTestCalls,
                info->minAbsInputA, info->minAbsInputB, info->minAbsErr, info->maxAbsInputA, info->maxAbsInputB, info->maxAbsErr,
                info->totalRelErr * oneOverTestCalls,
                info->minRelInputA, info->minRelInputB, info->minRelErr, info->maxRelInputA, info->maxRelInputB, info->maxRelErr);
    }
}

internal void
print_speed_info(String name, u32 maxNameSize, char *func, u32 tests, f32 seconds, f64 sum, f32 secondsRatio, f64 extraCall = 0.0)
{
    volatile f64 x = sum;
    unused(x);
    f64 testCalls = (f64)tests;
    persist char *spaces = "                                                                ";
    if (extraCall)
    {
        fprintf(stdout, "%.*s%.*s: %f, %f sec, %10.0f %s/sec, avg %f usec (%f usec/call)\n",
                STR_FMT(name), (u32)(maxNameSize - name.size), spaces, secondsRatio,
                seconds, testCalls / seconds, func, seconds / (testCalls / 1000000.0), seconds / (extraCall * testCalls / 1000000.0));
    }
    else
    {
        fprintf(stdout, "%.*s%.*s: %f, %f sec, %10.0f %s/sec, avg %f usec\n",
                STR_FMT(name), (u32)(maxNameSize - name.size), spaces, secondsRatio,
                seconds, testCalls / seconds, func, seconds / (testCalls / 1000000.0));
    }
}

internal void
update_comp(CompInfo *info, f64 origValue, f64 resultValue, f32 inputA, f32 inputB = 0.0f)
{
    f64 absErr = origValue - resultValue;
    if (absolute(info->minAbsErr) > absolute(absErr))
    {
        info->minAbsErr = absErr;
        info->minAbsInputA = inputA;
        info->minAbsInputB = inputB;
    }
    if (absolute(info->maxAbsErr) < absolute(absErr))
    {
        info->maxAbsErr = absErr;
        info->maxAbsInputA = inputA;
        info->maxAbsInputB = inputB;
    }
    info->totalAbsErr += absolute(absErr);
    if (origValue)
    {
        //f64 relErr = (origValue - resultValue) / origValue;
        f64 relErr = 1.0 - resultValue / origValue;
        if (absolute(info->minRelErr) > absolute(relErr))
        {
            info->minRelErr = relErr;
            info->minRelInputA = inputA;
            info->minRelInputB = inputB;
        }
        if (absolute(info->maxRelErr) < absolute(relErr))
        {
            info->maxRelErr = relErr;
            info->maxRelInputA = inputA;
            info->maxRelInputB = inputB;
        }
        info->totalRelErr += absolute(relErr);
    }
}

internal f32
run_comp_f32(String name, u32 maxNameSize, char *func, u32 tests, f32 minVal, f32 maxVal,
             MathFuncF32FromF32 *origFunc, MathFuncF32FromF32 *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 1;
    compInfo.minAbsErr = F32_MAX;
    compInfo.minRelErr = F32_MAX;

    f32 oneOverTests = 1.0f / (f32)tests;
    f32 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f32 inputVal = (f32)index * scale + minVal;
        f32 origRes = origFunc(inputVal);
        f32 testRes = testFunc(inputVal);
        update_comp(&compInfo, origRes, testRes, inputVal);
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    if (secondsBase == 0.0f) {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, 1.0f);
    } else {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, seconds / secondsBase);
    }
    return seconds;
}

internal f32
run_comp_f32_x(String name, u32 maxNameSize, char *func, u32 tests, f32 minVal, f32 maxVal,
               MathFuncF64FromF64 *origFunc, MathFuncF32FromF32 *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 1;
    compInfo.minAbsErr = F32_MAX;
    compInfo.minRelErr = F32_MAX;

    f32 oneOverTests = 1.0f / (f32)tests;
    f32 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f32 inputVal = (f32)index * scale + minVal;
        f64 origRes = origFunc(inputVal);
        f32 testRes = testFunc(inputVal);
        update_comp(&compInfo, origRes, testRes, inputVal);
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    if (secondsBase == 0.0f) {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, 1.0f);
    } else {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, seconds / secondsBase);
    }
    return seconds;
}

global volatile f64 gRunSpeedSum;
internal f32
run_speed_f32(String name, u32 maxNameSize, char *func, u32 tests, f32 minVal, f32 maxVal,
              MathFuncF32FromF32 *testFunc, f32 secondsBase)
{
    gRunSpeedSum = 0.0;
    f32 oneOverTests = 1.0f / (f32)tests;
    f32 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f32 inputVal = (f32)index * scale + minVal;
        f32 testRes0 = testFunc(inputVal);
        inputVal = clamp(minVal, inputVal * testRes0 + inputVal, maxVal);
        f32 testRes1 = testFunc(inputVal);
        inputVal = clamp(minVal, inputVal * testRes1 + inputVal, maxVal);
        f32 testRes2 = testFunc(inputVal);
        inputVal = clamp(minVal, inputVal * testRes2 + inputVal, maxVal);
        f32 testRes3 = testFunc(inputVal);
        gRunSpeedSum += inputVal + testRes3;
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    if (secondsBase == 0.0f) {
        print_speed_info(name, maxNameSize, func, 4*tests, seconds, gRunSpeedSum, 1.0f);
    } else {
        print_speed_info(name, maxNameSize, func, 4*tests, seconds, gRunSpeedSum, seconds / secondsBase);
    }
    return seconds;
}

internal f32
run_comp_f32_4x(String name, u32 maxNameSize, char *func, u32 tests, f32 minVal, f32 maxVal,
                MathFuncF32FromF32 *origFunc, MathFuncF32_4xFromF32_4x *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 1;
    compInfo.minAbsErr = F32_MAX;
    compInfo.minRelErr = F32_MAX;

    f32 oneOverTests = 1.0f / (f32)tests;
    f32 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f32 baseInput = (f32)index * scale + minVal;
        f32_4x inputVal = F32_4x(baseInput, 1.0e-9f * baseInput, 1.0e9f * baseInput, -1.0e9f * baseInput);
        f32 origRes = origFunc(baseInput);
        f32_4x testRes = testFunc(inputVal);
        update_comp(&compInfo, origRes, testRes.e[0], baseInput);

    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    if (secondsBase == 0.0f) {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, 1.0f);
    } else {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, seconds / secondsBase);
    }
    return seconds;
}

internal f32
run_comp_f32_4x_x(String name, u32 maxNameSize, char *func, u32 tests, f32 minVal, f32 maxVal,
                  MathFuncF64FromF64 *origFunc, MathFuncF32_4xFromF32_4x  *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 1;
    compInfo.minAbsErr = F32_MAX;
    compInfo.minRelErr = F32_MAX;

    f32 oneOverTests = 1.0f / (f32)tests;
    f32 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f32 baseInput = (f32)index * scale + minVal;
        f32_4x inputVal = F32_4x(baseInput, 1.0e-9f * baseInput, 1.0e9f * baseInput, -1.0e9f * baseInput);
        f64 origRes = origFunc(baseInput);
        f32_4x testRes = testFunc(inputVal);
        update_comp(&compInfo, origRes, testRes.e[0], baseInput);
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    if (secondsBase == 0.0f) {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, 1.0f);
    } else {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, seconds / secondsBase);
    }
    return seconds;
}

global volatile f32_4x gRunSpeedSum4x;
internal f32
run_speed_f32_4x(String name, u32 maxNameSize, char *func, u32 tests, f32 minVal, f32 maxVal,
                 MathFuncF32_4xFromF32_4x *testFunc, f32 secondsBase)
{
    gRunSpeedSum4x.m = _mm_setzero_ps();
    f32 oneOverTests = 1.0f / (f32)tests;
    f32 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f32 baseInput = (f32)index * scale + minVal;
        f32_4x inputVal = F32_4x(baseInput, 1.0e-9f * baseInput, 1.0e9f * baseInput, -1.0e9f * baseInput);
        f32_4x testRes0 = testFunc(inputVal);
        inputVal = inputVal * testRes0 + inputVal;
        f32_4x testRes1 = testFunc(inputVal);
        inputVal = inputVal * testRes1 + testRes0;
        f32_4x testRes2 = testFunc(inputVal);
        inputVal = inputVal * testRes2 + testRes1;
        f32_4x testRes3 = testFunc(inputVal);
        gRunSpeedSum4x.m = _mm_add_ps(gRunSpeedSum4x.m, (inputVal + testRes3).m);
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    if (secondsBase == 0.0f) {
        print_speed_info(name, maxNameSize, func, 4*tests, seconds, gRunSpeedSum4x.e[0], 1.0f, 4.0);
    } else {
        print_speed_info(name, maxNameSize, func, 4*tests, seconds, gRunSpeedSum4x.e[0], seconds / secondsBase, 4.0);
    }
    return seconds;
}

internal f32
run_comp_f32_f32(String name, u32 maxNameSize, char *func, u32 testsA, f32 minValA, f32 maxValA, u32 testsB, f32 minValB, f32 maxValB,
                 MathFuncF32FromF32F32 *origFunc, MathFuncF32FromF32F32 *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 2;
    compInfo.minAbsErr = F32_MAX;
    compInfo.minRelErr = F32_MAX;

    f32 scaleA = (maxValA - minValA) / (f32)testsA;
    f32 scaleB = (maxValB - minValB) / (f32)testsB;
    struct timespec start = linux_get_wall_clock();
    for (u32 aIdx = 0; aIdx < testsA; ++aIdx)
    {
        f32 inputValA = (f32)aIdx * scaleA + minValA;
        for (u32 bIdx = 0; bIdx < testsB; ++bIdx)
        {
            f32 inputValB = (f32)bIdx * scaleB + minValB;
            f32 origRes = origFunc(inputValA, inputValB);
            f32 testRes = testFunc(inputValA, inputValB);
            update_comp(&compInfo, origRes, testRes, inputValA, inputValB);
        }
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 tests = testsA * testsB;
    if (secondsBase == 0.0f) {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, 1.0f);
    } else {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, seconds / secondsBase);
    }
    return seconds;
}

internal f32
run_comp_f32_f32_x(String name, u32 maxNameSize, char *func, u32 testsA, f32 minValA, f32 maxValA, u32 testsB, f32 minValB, f32 maxValB,
                   MathFuncF64FromF64F64 *origFunc, MathFuncF32FromF32F32 *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 2;
    compInfo.minAbsErr = F32_MAX;
    compInfo.minRelErr = F32_MAX;

    f32 scaleA = (maxValA - minValA) / (f32)testsA;
    f32 scaleB = (maxValB - minValB) / (f32)testsB;
    struct timespec start = linux_get_wall_clock();
    for (u32 aIdx = 0; aIdx < testsA; ++aIdx)
    {
        f32 inputValA = (f32)aIdx * scaleA + minValA;
        for (u32 bIdx = 0; bIdx < testsB; ++bIdx)
        {
            f32 inputValB = (f32)bIdx * scaleB + minValB;
            f64 origRes = origFunc(inputValA, inputValB);
            f32 testRes = testFunc(inputValA, inputValB);
            update_comp(&compInfo, origRes, testRes, inputValA, inputValB);
        }
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 tests = testsA * testsB;
    if (secondsBase == 0.0f) {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, 1.0f);
    } else {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, seconds / secondsBase);
    }
    return seconds;
}

global volatile f64 gRunSpeedSumF32F32;
internal f32
run_speed_f32_f32(String name, u32 maxNameSize, char *func, u32 testsA, f32 minValA, f32 maxValA, u32 testsB, f32 minValB, f32 maxValB,
                  MathFuncF32FromF32F32 *testFunc, f32 secondsBase)
{
    gRunSpeedSumF32F32 = 0.0;
    f32 scaleA = (maxValA - minValA) / (f32)testsA;
    f32 scaleB = (maxValB - minValB) / (f32)testsB;
    struct timespec start = linux_get_wall_clock();
    for (u32 aIdx = 0; aIdx < testsA; ++aIdx)
    {
        f32 inputValA = (f32)aIdx * scaleA + minValA;
        for (u32 bIdx = 0; bIdx < testsB; ++bIdx)
        {
            f32 inputValB = (f32)bIdx * scaleB + minValB;
            f32 inputVal = inputValA;
            f32 testRes0 = testFunc(inputVal, inputValB);
            inputVal = clamp(minValA, inputVal * testRes0 + inputVal, maxValA);
            f32 testRes1 = testFunc(inputVal, inputValB);
            inputVal = clamp(minValA, inputVal * testRes1 + testRes0, maxValA);
            f32 testRes2 = testFunc(inputVal, inputValB);
            inputVal = clamp(minValA, inputVal * testRes2 + testRes1, maxValA);
            f32 testRes3 = testFunc(inputVal, inputValB);
            gRunSpeedSumF32F32 += inputVal + testRes3;
        }
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 tests = testsA * testsB;
    if (secondsBase == 0.0f) {
        print_speed_info(name, maxNameSize, func, 4*tests, seconds, gRunSpeedSumF32F32, 1.0f);
    } else {
        print_speed_info(name, maxNameSize, func, 4*tests, seconds, gRunSpeedSumF32F32, seconds / secondsBase);
    }
    return seconds;
}

internal f32
run_comp_f32_f32_4x(String name, u32 maxNameSize, char *func,
                    u32 testsA, f32 minValA, f32 maxValA, u32 testsB, f32 minValB, f32 maxValB,
                    MathFuncF32FromF32F32 *origFunc, MathFuncF32_4xFromF32_4xF32_4x *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 2;
    compInfo.minAbsErr = F32_MAX;
    compInfo.minRelErr = F32_MAX;

    f32 scaleA = (maxValA - minValA) / (f32)testsA;
    f32 scaleB = (maxValB - minValB) / (f32)testsB;
    struct timespec start = linux_get_wall_clock();
    for (u32 aIdx = 0; aIdx < testsA; ++aIdx)
    {
        f32 inputValA = (f32)aIdx * scaleA + minValA;
        f32_4x inputValA4 = F32_4x(inputValA, 1.0e-9f * inputValA, 1.0e9f * inputValA, -1.0e9f * inputValA);
        for (u32 bIdx = 0; bIdx < testsB; ++bIdx)
        {
            f32 inputValB = (f32)bIdx * scaleB + minValB;
            f32_4x inputValB4 = F32_4x(inputValB, 1.0e-9f * inputValB, 1.0e9f * inputValB, -1.0e9f * inputValB);
            f32 origRes = origFunc(inputValA, inputValB);
            f32_4x testRes = testFunc(inputValA4, inputValB4);
            update_comp(&compInfo, origRes, testRes.e[0], inputValA, inputValB);
        }
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 tests = testsA * testsB;
    if (secondsBase == 0.0f) {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, 1.0f);
    } else {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, seconds / secondsBase);
    }
    return seconds;
}

internal f32
run_comp_f32_f32_4x_x(String name, u32 maxNameSize, char *func,
                      u32 testsA, f32 minValA, f32 maxValA, u32 testsB, f32 minValB, f32 maxValB,
                      MathFuncF64FromF64F64 *origFunc, MathFuncF32_4xFromF32_4xF32_4x *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 2;
    compInfo.minAbsErr = F32_MAX;
    compInfo.minRelErr = F32_MAX;

    f32 scaleA = (maxValA - minValA) / (f32)testsA;
    f32 scaleB = (maxValB - minValB) / (f32)testsB;
    struct timespec start = linux_get_wall_clock();
    for (u32 aIdx = 0; aIdx < testsA; ++aIdx)
    {
        f32 inputValA = (f32)aIdx * scaleA + minValA;
        f32_4x inputValA4 = F32_4x(inputValA, 1.0e-9f * inputValA, 1.0e9f * inputValA, -1.0e9f * inputValA);
        for (u32 bIdx = 0; bIdx < testsB; ++bIdx)
        {
            f32 inputValB = (f32)bIdx * scaleB + minValB;
            f32_4x inputValB4 = F32_4x(inputValB, 1.0e-9f * inputValB, 1.0e9f * inputValB, -1.0e9f * inputValB);
            f64 origRes = origFunc(inputValA, inputValB);
            f32_4x testRes = testFunc(inputValA4, inputValB4);
            update_comp(&compInfo, origRes, testRes.e[0], inputValA, inputValB);
        }
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 tests = testsA * testsB;
    if (secondsBase == 0.0f) {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, 1.0f);
    } else {
        print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, seconds / secondsBase);
    }
    return seconds;
}

global volatile f32_4x gRunSpeedSumF32F324x;
internal f32
run_speed_f32_f32_4x(String name, u32 maxNameSize, char *func, u32 testsA, f32 minValA, f32 maxValA, u32 testsB, f32 minValB, f32 maxValB,
                     MathFuncF32_4xFromF32_4xF32_4x *testFunc, f32 secondsBase)
{
    gRunSpeedSumF32F324x.m = _mm_setzero_ps();
    f32 scaleA = (maxValA - minValA) / (f32)testsA;
    f32 scaleB = (maxValB - minValB) / (f32)testsB;
    struct timespec start = linux_get_wall_clock();
    for (u32 aIdx = 0; aIdx < testsA; ++aIdx)
    {
        f32 inputValA = (f32)aIdx * scaleA + minValA;
        f32_4x inputValA4 = F32_4x(inputValA, 1.0e-9f * inputValA, 1.0e9f * inputValA, -1.0e9f * inputValA);
        for (u32 bIdx = 0; bIdx < testsB; ++bIdx)
        {
            f32 inputValB = (f32)bIdx * scaleB + minValB;
            f32_4x inputValB4 = F32_4x(inputValB, 1.0e-9f * inputValB, 1.0e9f * inputValB, -1.0e9f * inputValB);
            f32_4x inputVal = inputValA4;
            f32_4x testRes0 = testFunc(inputVal, inputValB4);
            f32_4x testRes1 = testFunc(inputVal, inputValB4);
            inputVal = inputVal * testRes1 + testRes0;
            f32_4x testRes2 = testFunc(inputVal, inputValB4);
            inputVal = inputVal * testRes2 + testRes1;
            f32_4x testRes3 = testFunc(inputVal, inputValB4);
            gRunSpeedSumF32F324x.m = _mm_add_ps(gRunSpeedSumF32F324x.m, (inputVal + testRes3).m);
        }
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 tests = testsA * testsB;
    if (secondsBase == 0.0f) {
        print_speed_info(name, maxNameSize, func, 4*tests, seconds, gRunSpeedSumF32F324x.e[0], 1.0f, 4.0);
    } else {
        print_speed_info(name, maxNameSize, func, 4*tests, seconds, gRunSpeedSumF32F324x.e[0], seconds / secondsBase, 4.0);
    }
    return seconds;
}
