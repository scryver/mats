#define MATH_FUNC_F32_FROM_F32(name)                f32 name(f32 x)
typedef MATH_FUNC_F32_FROM_F32(MathFuncF32FromF32);
#define MATH_FUNC_F32_FROM_F32_F32(name)            f32 name(f32 a, f32 b)
typedef MATH_FUNC_F32_FROM_F32_F32(MathFuncF32FromF32F32);
#define MATH_FUNC_F64_FROM_F64(name)                f64 name(f64 x)
typedef MATH_FUNC_F64_FROM_F64(MathFuncF64FromF64);

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
    f32 minAbsInput;
    f32 maxAbsInput;
    f32 minRelInput;
    f32 maxRelInput;
};

internal void
print_comp_info(String name, u32 maxNameSize, char *func, u32 tests, f32 seconds, CompInfo *info, f32 secondsRatio)
{
    f64 testCalls = (f64)tests;
    f64 oneOverTestCalls = 1.0 / testCalls;
    persist char *spaces = "                                                                ";
    fprintf(stdout, "%.*s%.*s: %f, %f sec, %10.0f %s/sec, avg %f usec\n    abs %e (min(%a): %e, max(%a): %e)\n    rel %e (min(%a): %e, max(%a): %e)\n",
            STR_FMT(name), (u32)(maxNameSize - name.size), spaces, secondsRatio,
            seconds, testCalls / seconds, func, seconds / (testCalls / 1000000.0),
            info->totalAbsErr * oneOverTestCalls, info->minAbsInput, info->minAbsErr, info->maxAbsInput, info->maxAbsErr,
            info->totalRelErr * oneOverTestCalls, info->minRelInput, info->minRelErr, info->maxRelInput, info->maxRelErr);
}

internal void
print_speed_info(String name, u32 maxNameSize, char *func, u32 tests, f32 seconds, f64 sum, f32 secondsRatio)
{
    volatile f64 x = sum;
    unused(x);
    f64 testCalls = (f64)tests;
    persist char *spaces = "                                                                ";
    fprintf(stdout, "%.*s%.*s: %f, %f sec, %10.0f %s/sec, avg %f usec\n",
            STR_FMT(name), (u32)(maxNameSize - name.size), spaces, secondsRatio,
            seconds, testCalls / seconds, func, seconds / (testCalls / 1000000.0));
}

internal f32
run_comp_f32(String name, u32 maxNameSize, char *func, u32 tests, f32 minVal, f32 maxVal,
             MathFuncF32FromF32 *origFunc, MathFuncF32FromF32 *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.minAbsErr = F32_MAX;
    compInfo.maxAbsErr = -F32_MAX;
    compInfo.minRelErr = F32_MAX;
    compInfo.maxRelErr = -F32_MAX;

    f32 oneOverTests = 1.0f / (f32)tests;
    f32 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f32 inputVal = (f32)index * scale + minVal;
        f32 origRes = origFunc(inputVal);
        f32 testRes = testFunc(inputVal);

        f64 absErr = origRes - testRes;
        if (compInfo.minAbsErr > absErr) {
            compInfo.minAbsErr = absErr;
            compInfo.minAbsInput = inputVal;
        }
        if (compInfo.maxAbsErr < absErr) {
            compInfo.maxAbsErr = absErr;
            compInfo.maxAbsInput = inputVal;
        }
        compInfo.totalAbsErr += absErr;

        if (inputVal)
        {
            f64 relErr = (origRes - testRes) / inputVal;
            if (compInfo.minRelErr > relErr) {
                compInfo.minRelErr = relErr;
                compInfo.minRelInput = inputVal;
            }
            if (compInfo.maxRelErr < relErr) {
                compInfo.maxRelErr = relErr;
                compInfo.maxRelInput = inputVal;
            }
            compInfo.totalRelErr += relErr;
        }
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
    compInfo.minAbsErr = F32_MAX;
    compInfo.maxAbsErr = -F32_MAX;
    compInfo.minRelErr = F32_MAX;
    compInfo.maxRelErr = -F32_MAX;

    f32 oneOverTests = 1.0f / (f32)tests;
    f32 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f32 inputVal = (f32)index * scale + minVal;
        f64 origRes = origFunc(inputVal);
        f32 testRes = testFunc(inputVal);

        f64 absErr = origRes - testRes;
        if (compInfo.minAbsErr > absErr) {
            compInfo.minAbsErr = absErr;
            compInfo.minAbsInput = inputVal;
        }
        if (compInfo.maxAbsErr < absErr) {
            compInfo.maxAbsErr = absErr;
            compInfo.maxAbsInput = inputVal;
        }
        compInfo.totalAbsErr += absErr;

        if (inputVal)
        {
            f64 relErr = (origRes - testRes) / inputVal;
            if (compInfo.minRelErr > relErr) {
                compInfo.minRelErr = relErr;
                compInfo.minRelInput = inputVal;
            }
            if (compInfo.maxRelErr < relErr) {
                compInfo.maxRelErr = relErr;
                compInfo.maxRelInput = inputVal;
            }
            compInfo.totalRelErr += relErr;
        }
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
        f32 testRes1 = testFunc(testRes0);
        f32 testRes2 = testFunc(testRes1);
        f32 testRes3 = testFunc(testRes2);
        gRunSpeedSum += testRes3;
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    if (secondsBase == 0.0f) {
        print_speed_info(name, maxNameSize, func, 4*tests, seconds, gRunSpeedSum, 1.0f);
    } else {
        print_speed_info(name, maxNameSize, func, 4*tests, seconds, gRunSpeedSum, seconds / secondsBase);
    }
    return seconds;
}
