//
// NOTE(michiel): 32-bit
//
#define MATH_FUNC_F32_FROM_F32(name)                f32 name(f32 x)
typedef MATH_FUNC_F32_FROM_F32(MathFuncF32FromF32);
#define MATH_FUNC_F32_FROM_F32_F32(name)            f32 name(f32 a, f32 b)
typedef MATH_FUNC_F32_FROM_F32_F32(MathFuncF32FromF32F32);
#define MATH_FUNC_F32_F32_FROM_F32(name)            void name(f32 a, f32 *x, f32 *y)
typedef MATH_FUNC_F32_F32_FROM_F32(MathFuncF32F32FromF32);
#define MATH_FUNC_F32_4x_FROM_F32_4x(name)          f32_4x name(f32_4x x)
typedef MATH_FUNC_F32_4x_FROM_F32_4x(MathFuncF32_4xFromF32_4x);
#define MATH_FUNC_F32_4x_FROM_F32_4x_F32_4x(name)   f32_4x name(f32_4x a, f32_4x b)
typedef MATH_FUNC_F32_4x_FROM_F32_4x_F32_4x(MathFuncF32_4xFromF32_4xF32_4x);

#define call_comp(lib, name, suf, sec)  \
run_comp_f32(string(#lib " " #name "32"), 30, #name, tests, minVal, maxVal, name, name ## suf, sec)
#define call_spd(lib, name, suf, sec) \
run_speed_f32(string(#lib " " #name "32"), 30, #name, tests, minVal, maxVal, name ## suf, sec)

#define call_comp_r2(lib, name, suf, sec)  \
run_comp_f32_f32(string(#lib " " #name "32"), 30, #name, tests, minVal, maxVal, name, name ## suf, sec)

#define call_comp2(lib, name, suf, sec)  \
run_comp_f32_f32(string(#lib " " #name "32"), 30, #name, testsA, minValA, maxValA, testsB, minValB, maxValB, name, name ## suf, sec)
#define call_spd2(lib, name, suf, sec) \
run_speed_f32_f32(string(#lib " " #name "32"), 30, #name, testsA, minValA, maxValA, testsB, minValB, maxValB, name ## suf, sec)

#define call_comp_4x(lib, name, suf, sec)  \
run_comp_f32_4x(string(#lib " " #name "32"), 30, #name "_4x", tests, minVal, maxVal, name, name ## suf, sec)
#define call_spd_4x(lib, name, suf, sec) \
run_speed_f32_4x(string(#lib " " #name "32"), 30, #name "_4x", tests, minVal, maxVal, name ## suf, sec)

#define call_comp2_4x(lib, name, suf, sec)  \
run_comp_f32_f32_4x(string(#lib " " #name "32"), 30, #name "_4x", testsA, minValA, maxValA, testsB, minValB, maxValB, name, name ## suf, sec)
#define call_spd2_4x(lib, name, suf, sec) \
run_speed_f32_f32_4x(string(#lib " " #name "32"), 30, #name "_4x", testsA, minValA, maxValA, testsB, minValB, maxValB, name ## suf, sec)

#define WIDE_FUNC_FROM_F32(name) \
func f32_4x \
name##_4x(f32_4x x) \
{ \
f32_4x result; \
result.e[0] = name(x.e[0]); \
result.e[1] = name(x.e[1]); \
result.e[2] = name(x.e[2]); \
result.e[3] = name(x.e[3]); \
return result; \
}

#define WIDE_FUNC_FROM_F32_F32(name) \
func f32_4x \
name##_4x(f32_4x x, f32_4x y) \
{ \
f32_4x result; \
result.e[0] = name(x.e[0], y.e[0]); \
result.e[1] = name(x.e[1], y.e[1]); \
result.e[2] = name(x.e[2], y.e[2]); \
result.e[3] = name(x.e[3], y.e[3]); \
return result; \
}

#define BEGIN_TEST(d, name, func) \
if (d & DoTest_##name) { \
fprintf(stdout, #name "32\n"); \
f32 stdSec = func(stdlib, name, f, 0.0f)

#define BEGIN_TEST_WIDE(d, name, func) \
if (d & DoTest_##name) { \
fprintf(stdout, #name "32 wide\n"); \
f32 stdSec = func(stdlib, name, f_4x, 0.0f)

#define END_TEST() \
fprintf(stdout, "\n"); \
}

//
// NOTE(michiel): 64-bit
//
#define MATH_FUNC_F64_FROM_F64(name)                f64 name(f64 x)
typedef MATH_FUNC_F64_FROM_F64(MathFuncF64FromF64);
#define MATH_FUNC_F64_FROM_F64_F64(name)            f64 name(f64 a, f64 b)
typedef MATH_FUNC_F64_FROM_F64_F64(MathFuncF64FromF64F64);
#define MATH_FUNC_F64_F64_FROM_F64(name)            void name(f64 a, f64 *x, f64 *y)
typedef MATH_FUNC_F64_F64_FROM_F64(MathFuncF64F64FromF64);
#define MATH_FUNC_F64_2x_FROM_F64_2x(name)          f64_2x name(f64_2x x)
typedef MATH_FUNC_F64_2x_FROM_F64_2x(MathFuncF64_2xFromF64_2x);
#define MATH_FUNC_F64_2x_FROM_F64_2x_F64_2x(name)   f64_2x name(f64_2x a, f64_2x b)
typedef MATH_FUNC_F64_2x_FROM_F64_2x_F64_2x(MathFuncF64_2xFromF64_2xF64_2x);

#define call_comp64(lib, name, suf, sec)  \
run_comp_f64(string(#lib " " #name "64"), 30, #name, tests, minVal, maxVal, name, name ## suf, sec)
#define call_spd64(lib, name, suf, sec) \
run_speed_f64(string(#lib " " #name "64"), 30, #name, tests, minVal, maxVal, name ## suf, sec)

#define call_comp64_r2(lib, name, suf, sec)  \
run_comp_f64_f64(string(#lib " " #name "64"), 30, #name, tests, minVal, maxVal, name, name ## suf, sec)

#define call_comp64_2(lib, name, suf, sec)  \
run_comp_f64_f64(string(#lib " " #name "64"), 30, #name, testsA, minValA, maxValA, testsB, minValB, maxValB, name, name ## suf, sec)
#define call_spd64_2(lib, name, suf, sec) \
run_speed_f64_f64(string(#lib " " #name "64"), 30, #name, testsA, minValA, maxValA, testsB, minValB, maxValB, name ## suf, sec)

#define call_comp64_2x(lib, name, suf, sec)  \
run_comp_f64_2x(string(#lib " " #name "64"), 30, #name "_4x", tests, minVal, maxVal, name, name ## suf, sec)
#define call_spd64_2x(lib, name, suf, sec) \
run_speed_f64_2x(string(#lib " " #name "64"), 30, #name "_4x", tests, minVal, maxVal, name ## suf, sec)

#define call_comp64_2_2x(lib, name, suf, sec)  \
run_comp_f64_f64_2x(string(#lib " " #name "64"), 30, #name "_4x", testsA, minValA, maxValA, testsB, minValB, maxValB, name, name ## suf, sec)
#define call_spd64_2_2x(lib, name, suf, sec) \
run_speed_f64_f64_2x(string(#lib " " #name "64"), 30, #name "_4x", testsA, minValA, maxValA, testsB, minValB, maxValB, name ## suf, sec)

#define WIDE_FUNC_FROM_F64(name) \
func f64_2x \
name##_2x(f64_2x x) \
{ \
f64_2x result; \
result.e[0] = name(x.e[0]); \
result.e[1] = name(x.e[1]); \
return result; \
}

#define WIDE_FUNC_FROM_F64_F64(name) \
func f64_2x \
name##_2x(f64_2x x, f64_2x y) \
{ \
f64_2x result; \
result.e[0] = name(x.e[0], y.e[0]); \
result.e[1] = name(x.e[1], y.e[1]); \
return result; \
}

#define BEGIN_TEST64(d, name, func) \
if (d & DoTest_##name) { \
fprintf(stdout, #name "64\n"); \
f32 stdSec = func(stdlib, name, , 0.0f)

#define BEGIN_TEST_WIDE64(d, name, func) \
if (d & DoTest_##name) { \
fprintf(stdout, #name "64 wide\n"); \
f32 stdSec = func(stdlib, name, _2x, 0.0f)

func void
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

func f32
square(f32 x)
{
    return x * x;
}

func u32
string_length(const char *cString)
{
    const char *start = cString;
    while (cString && *cString++);
    return safe_truncate_to_u32(cString - start);
}

func String
string(char *cStr)
{
    return {string_length(cStr), (u8 *)cStr};
}

func b32
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

func struct timespec
linux_get_wall_clock()
{
    struct timespec clock;
    clock_gettime(CLOCK_MONOTONIC, &clock);
    return clock;
}

func f32
linux_get_seconds_elapsed(struct timespec start, struct timespec end)
{
    return ((f32)(end.tv_sec - start.tv_sec)
            + ((f32)(end.tv_nsec - start.tv_nsec) * 1e-9f));
}

struct CompInfo
{
    u32 numInputs;
    f64 totalAbsErr;
    f64 totalRelErr;
    f64 minAbsErr;
    f64 maxAbsErr;
    f64 minRelErr;
    f64 maxRelErr;
    f64 minAbsInputA;
    f64 minAbsInputB;
    f64 maxAbsInputA;
    f64 maxAbsInputB;
    f64 minRelInputA;
    f64 minRelInputB;
    f64 maxRelInputA;
    f64 maxRelInputB;
};

func void
print_comp_info(String name, u32 maxNameSize, char *func, u32 tests, f32 seconds, CompInfo *info, f32 secondsBase)
{
    f64 secondsRatio =  secondsBase == 0.0f ? 1.0 : (f64)secondsBase / (f64)seconds;
    f64 testCalls = (f64)tests;
    f64 oneOverTestCalls = 1.0 / testCalls;
    persist char *spaces = "                                                                ";
    if (info->numInputs == 1)
    {
        fprintf(stdout, "%.*s%.*s: %9.6f, %9.6f sec, %10.0f %s/sec, avg %f usec\n    abs %g (min(%g): %a, max(%g): %a)\n    rel %g (min(%g): %a, max(%g): %a)\n",
                STR_FMT(name), (u32)(maxNameSize - name.size), spaces, secondsRatio,
                seconds, testCalls / seconds, func, seconds / (testCalls / 1000000.0),
                info->totalAbsErr * oneOverTestCalls, info->minAbsInputA, info->minAbsErr, info->maxAbsInputA, info->maxAbsErr,
                info->totalRelErr * oneOverTestCalls, info->minRelInputA, info->minRelErr, info->maxRelInputA, info->maxRelErr);
    }
    else
    {
        i_expect(info->numInputs == 2);

        fprintf(stdout, "%.*s%.*s: %9.6f, %9.6f sec, %10.0f %s/sec, avg %f usec\n    abs %g (min(%g, %g): %a, max(%g, %g): %a)\n    rel %g (min(%g, %g): %a, max(%g, %g): %a)\n",
                STR_FMT(name), (u32)(maxNameSize - name.size), spaces, secondsRatio,
                seconds, testCalls / seconds, func, seconds / (testCalls / 1000000.0),
                info->totalAbsErr * oneOverTestCalls,
                info->minAbsInputA, info->minAbsInputB, info->minAbsErr, info->maxAbsInputA, info->maxAbsInputB, info->maxAbsErr,
                info->totalRelErr * oneOverTestCalls,
                info->minRelInputA, info->minRelInputB, info->minRelErr, info->maxRelInputA, info->maxRelInputB, info->maxRelErr);
    }
}

func void
print_speed_info(String name, u32 maxNameSize, char *func, u32 tests, f32 seconds, f64 sum, f32 secondsBase, f64 extraCall = 0.0)
{
    f64 secondsRatio = secondsBase == 0.0f ? 1.0 : (f64)secondsBase / (f64)seconds;
    volatile f64 x = sum;
    unused(x);
    f64 testCalls = (f64)tests;
    persist char *spaces = "                                                                ";
    if (extraCall)
    {
        fprintf(stdout, "%.*s%.*s: %9.6f, %9.6f sec, %10.0f %s/sec, avg %f usec (%f usec/call)\n",
                STR_FMT(name), (u32)(maxNameSize - name.size), spaces, secondsRatio,
                seconds, testCalls / seconds, func, seconds / (testCalls / 1000000.0), seconds / (extraCall * testCalls / 1000000.0));
    }
    else
    {
        fprintf(stdout, "%.*s%.*s: %9.6f, %9.6f sec, %10.0f %s/sec, avg %f usec\n",
                STR_FMT(name), (u32)(maxNameSize - name.size), spaces, secondsRatio,
                seconds, testCalls / seconds, func, seconds / (testCalls / 1000000.0));
    }
}

func void
update_comp(CompInfo *info, f64 origValue, f64 resultValue, f64 inputA, f64 inputB = 0.0)
{
    f64 absErr = origValue - resultValue;
    if (absolute64(info->minAbsErr) > absolute64(absErr))
    {
        info->minAbsErr = absErr;
        info->minAbsInputA = inputA;
        info->minAbsInputB = inputB;
    }
    if (absolute64(info->maxAbsErr) < absolute64(absErr))
    {
        info->maxAbsErr = absErr;
        info->maxAbsInputA = inputA;
        info->maxAbsInputB = inputB;
    }
    info->totalAbsErr += absolute64(absErr);
    if (origValue)
    {
        //f64 relErr = (origValue - resultValue) / origValue;
        f64 relErr = 1.0 - resultValue / origValue;
        if (absolute64(info->minRelErr) > absolute64(relErr))
        {
            info->minRelErr = relErr;
            info->minRelInputA = inputA;
            info->minRelInputB = inputB;
        }
        if (absolute64(info->maxRelErr) < absolute64(relErr))
        {
            info->maxRelErr = relErr;
            info->maxRelInputA = inputA;
            info->maxRelInputB = inputB;
        }
        info->totalRelErr += absolute64(relErr);
    }
}

func f32
run_comp_f32(String name, u32 maxNameSize, char *func, u32 tests, f32 minVal, f32 maxVal,
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

    print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, secondsBase);
    return seconds;
}

func f32
run_comp_f64(String name, u32 maxNameSize, char *func, u32 tests, f64 minVal, f64 maxVal,
             MathFuncF64FromF64 *origFunc, MathFuncF64FromF64 *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 1;
    compInfo.minAbsErr = F64_MAX;
    compInfo.minRelErr = F64_MAX;

    f64 oneOverTests = 1.0 / (f64)tests;
    f64 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f64 inputVal = (f64)index * scale + minVal;
        f64 origRes = origFunc(inputVal);
        f64 testRes = testFunc(inputVal);
        update_comp(&compInfo, origRes, testRes, inputVal);
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, secondsBase);
    return seconds;
}

func f32
run_comp_f32_f32(String name, u32 maxNameSize, char *func, u32 tests, f32 minVal, f32 maxVal,
                 MathFuncF64F64FromF64 *origFunc, MathFuncF32F32FromF32 *testFunc, f32 secondsBase)
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
        f64 orig0, orig1;
        origFunc(inputVal, &orig0, &orig1);
        f32 test0, test1;
        testFunc(inputVal, &test0, &test1);
        update_comp(&compInfo, orig0, test0, inputVal);
        update_comp(&compInfo, orig1, test1, inputVal);
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, secondsBase);
    return seconds;
}

func f32
run_comp_f64_f64(String name, u32 maxNameSize, char *func, u32 tests, f64 minVal, f64 maxVal,
                 MathFuncF64F64FromF64 *origFunc, MathFuncF64F64FromF64 *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 1;
    compInfo.minAbsErr = F64_MAX;
    compInfo.minRelErr = F64_MAX;

    f64 oneOverTests = 1.0 / (f64)tests;
    f64 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f64 inputVal = (f64)index * scale + minVal;
        f64 orig0, orig1;
        origFunc(inputVal, &orig0, &orig1);
        f64 test0, test1;
        testFunc(inputVal, &test0, &test1);
        update_comp(&compInfo, orig0, test0, inputVal);
        update_comp(&compInfo, orig1, test1, inputVal);
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, secondsBase);
    return seconds;
}

global volatile f64 gRunSpeedSum;
func f32
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

    print_speed_info(name, maxNameSize, func, tests, seconds, gRunSpeedSum, secondsBase);
    return seconds;
}
;
func f32
run_speed_f64(String name, u32 maxNameSize, char *func, u32 tests, f64 minVal, f64 maxVal,
              MathFuncF64FromF64 *testFunc, f32 secondsBase)
{
    gRunSpeedSum = 0.0;
    f64 oneOverTests = 1.0 / (f64)tests;
    f64 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f64 inputVal = (f64)index * scale + minVal;
        f64 testRes0 = testFunc(inputVal);
        inputVal = clamp(minVal, inputVal * testRes0 + inputVal, maxVal);
        f64 testRes1 = testFunc(inputVal);
        inputVal = clamp(minVal, inputVal * testRes1 + inputVal, maxVal);
        f64 testRes2 = testFunc(inputVal);
        inputVal = clamp(minVal, inputVal * testRes2 + inputVal, maxVal);
        f64 testRes3 = testFunc(inputVal);
        gRunSpeedSum += inputVal + testRes3;
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    print_speed_info(name, maxNameSize, func, tests, seconds, gRunSpeedSum, secondsBase);
    return seconds;
}

func f32
run_comp_f32_4x(String name, u32 maxNameSize, char *func, u32 tests, f32 minVal, f32 maxVal,
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

    print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, secondsBase);
    return seconds;
}

func f32
run_comp_f64_2x(String name, u32 maxNameSize, char *func, u32 tests, f64 minVal, f64 maxVal,
                MathFuncF64FromF64 *origFunc, MathFuncF64_2xFromF64_2x  *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 1;
    compInfo.minAbsErr = F32_MAX;
    compInfo.minRelErr = F32_MAX;

    f64 oneOverTests = 1.0 / (f64)tests;
    f64 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f64 baseInput = (f64)index * scale + minVal;
        f64_2x inputVal = F64_2x(baseInput, 1.0e-9f * baseInput);
        f64 origRes = origFunc(baseInput);
        f64_2x testRes = testFunc(inputVal);
        update_comp(&compInfo, origRes, testRes.e[0], baseInput);
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, secondsBase);
    return seconds;
}

global volatile f32_4x gRunSpeedSum4x;
func f32
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
        f32_4x inputVal = F32_4x(baseInput, 1.0e-3f * baseInput, 1.0e3f * baseInput, -1.0e3f * baseInput);
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

    print_speed_info(name, maxNameSize, func, tests, seconds, gRunSpeedSum4x.e[0], secondsBase, 4.0);
    return seconds;
}

global volatile f64_2x gRunSpeedSum2x;
func f32
run_speed_f64_2x(String name, u32 maxNameSize, char *func, u32 tests, f64 minVal, f64 maxVal,
                 MathFuncF64_2xFromF64_2x *testFunc, f32 secondsBase)
{
    gRunSpeedSum2x.md = _mm_setzero_pd();
    f64 oneOverTests = 1.0 / (f64)tests;
    f64 scale = (maxVal - minVal) * oneOverTests;
    struct timespec start = linux_get_wall_clock();
    for (u32 index = 0; index < tests; ++index)
    {
        f64 baseInput = (f64)index * scale + minVal;
        f64_2x inputVal = F64_2x(baseInput, 1.0e-3f * baseInput);
        f64_2x testRes0 = testFunc(inputVal);
        inputVal = inputVal * testRes0 + inputVal;
        f64_2x testRes1 = testFunc(inputVal);
        inputVal = inputVal * testRes1 + testRes0;
        f64_2x testRes2 = testFunc(inputVal);
        inputVal = inputVal * testRes2 + testRes1;
        f64_2x testRes3 = testFunc(inputVal);
        gRunSpeedSum2x.md = _mm_add_pd(gRunSpeedSum2x.md, (inputVal + testRes3).md);
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    print_speed_info(name, maxNameSize, func, tests, seconds, gRunSpeedSum2x.e[0], secondsBase, 4.0);
    return seconds;
}

func f32
run_comp_f32_f32(String name, u32 maxNameSize, char *func, u32 testsA, f32 minValA, f32 maxValA, u32 testsB, f32 minValB, f32 maxValB,
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
    print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, secondsBase);
    return seconds;
}

func f32
run_comp_f64_f64(String name, u32 maxNameSize, char *func, u32 testsA, f64 minValA, f64 maxValA, u32 testsB, f64 minValB, f64 maxValB,
                 MathFuncF64FromF64F64 *origFunc, MathFuncF64FromF64F64 *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 2;
    compInfo.minAbsErr = F32_MAX;
    compInfo.minRelErr = F32_MAX;

    f64 scaleA = (maxValA - minValA) / (f64)testsA;
    f64 scaleB = (maxValB - minValB) / (f64)testsB;
    struct timespec start = linux_get_wall_clock();
    for (u32 aIdx = 0; aIdx < testsA; ++aIdx)
    {
        f64 inputValA = (f64)aIdx * scaleA + minValA;
        for (u32 bIdx = 0; bIdx < testsB; ++bIdx)
        {
            f64 inputValB = (f64)bIdx * scaleB + minValB;
            f64 origRes = origFunc(inputValA, inputValB);
            f64 testRes = testFunc(inputValA, inputValB);
            update_comp(&compInfo, origRes, testRes, inputValA, inputValB);
        }
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 tests = testsA * testsB;
    print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, secondsBase);
    return seconds;
}

global volatile f64 gRunSpeedSumF32F32;
func f32
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
    print_speed_info(name, maxNameSize, func, tests, seconds, gRunSpeedSumF32F32, secondsBase);
    return seconds;
}

global volatile f64 gRunSpeedSumF64F64;
func f32
run_speed_f64_f64(String name, u32 maxNameSize, char *func, u32 testsA, f64 minValA, f64 maxValA, u32 testsB, f64 minValB, f64 maxValB,
                  MathFuncF64FromF64F64 *testFunc, f32 secondsBase)
{
    gRunSpeedSumF32F32 = 0.0;
    f64 scaleA = (maxValA - minValA) / (f64)testsA;
    f64 scaleB = (maxValB - minValB) / (f64)testsB;
    struct timespec start = linux_get_wall_clock();
    for (u32 aIdx = 0; aIdx < testsA; ++aIdx)
    {
        f64 inputValA = (f64)aIdx * scaleA + minValA;
        for (u32 bIdx = 0; bIdx < testsB; ++bIdx)
        {
            f64 inputValB = (f64)bIdx * scaleB + minValB;
            f64 inputVal = inputValA;
            f64 testRes0 = testFunc(inputVal, inputValB);
            inputVal = clamp(minValA, inputVal * testRes0 + inputVal, maxValA);
            f64 testRes1 = testFunc(inputVal, inputValB);
            inputVal = clamp(minValA, inputVal * testRes1 + testRes0, maxValA);
            f64 testRes2 = testFunc(inputVal, inputValB);
            inputVal = clamp(minValA, inputVal * testRes2 + testRes1, maxValA);
            f64 testRes3 = testFunc(inputVal, inputValB);
            gRunSpeedSumF64F64 += inputVal + testRes3;
        }
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 tests = testsA * testsB;
    print_speed_info(name, maxNameSize, func, tests, seconds, gRunSpeedSumF64F64, secondsBase);
    return seconds;
}

func f32
run_comp_f32_f32_4x(String name, u32 maxNameSize, char *func,
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
    print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, secondsBase);
    return seconds;
}

func f32
run_comp_f64_f64_2x(String name, u32 maxNameSize, char *func,
                    u32 testsA, f64 minValA, f64 maxValA, u32 testsB, f64 minValB, f64 maxValB,
                    MathFuncF64FromF64F64 *origFunc, MathFuncF64_2xFromF64_2xF64_2x *testFunc, f32 secondsBase)
{
    CompInfo compInfo = {};
    compInfo.numInputs = 2;
    compInfo.minAbsErr = F32_MAX;
    compInfo.minRelErr = F32_MAX;

    f64 scaleA = (maxValA - minValA) / (f64)testsA;
    f64 scaleB = (maxValB - minValB) / (f64)testsB;
    struct timespec start = linux_get_wall_clock();
    for (u32 aIdx = 0; aIdx < testsA; ++aIdx)
    {
        f64 inputValA = (f64)aIdx * scaleA + minValA;
        f64_2x inputValA2 = F64_2x(inputValA, 1.0e-9f * inputValA);
        for (u32 bIdx = 0; bIdx < testsB; ++bIdx)
        {
            f64 inputValB = (f64)bIdx * scaleB + minValB;
            f64_2x inputValB2 = F64_2x(inputValB, 1.0e-9f * inputValB);
            f64 origRes = origFunc(inputValA, inputValB);
            f64_2x testRes = testFunc(inputValA2, inputValB2);
            update_comp(&compInfo, origRes, testRes.e[0], inputValA, inputValB);
        }
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 tests = testsA * testsB;
    print_comp_info(name, maxNameSize, func, tests, seconds, &compInfo, secondsBase);
    return seconds;
}

global volatile f32_4x gRunSpeedSumF32F324x;
func f32
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
    print_speed_info(name, maxNameSize, func, tests, seconds, gRunSpeedSumF32F324x.e[0], secondsBase, 4.0);
    return seconds;
}

global volatile f64_2x gRunSpeedSumF64F642x;
func f32
run_speed_f64_f64_2x(String name, u32 maxNameSize, char *func,
                     u32 testsA, f64 minValA, f64 maxValA, u32 testsB, f64 minValB, f64 maxValB,
                     MathFuncF64_2xFromF64_2xF64_2x *testFunc, f32 secondsBase)
{
    gRunSpeedSumF64F642x.md = _mm_setzero_pd();
    f64 scaleA = (maxValA - minValA) / (f64)testsA;
    f64 scaleB = (maxValB - minValB) / (f64)testsB;
    struct timespec start = linux_get_wall_clock();
    for (u32 aIdx = 0; aIdx < testsA; ++aIdx)
    {
        f64 inputValA = (f64)aIdx * scaleA + minValA;
        f64_2x inputValA2 = F64_2x(inputValA, 1.0e-9f * inputValA);
        for (u32 bIdx = 0; bIdx < testsB; ++bIdx)
        {
            f64 inputValB = (f64)bIdx * scaleB + minValB;
            f64_2x inputValB2 = F64_2x(inputValB, 1.0e-9f * inputValB);
            f64_2x inputVal = inputValA2;
            f64_2x testRes0 = testFunc(inputVal, inputValB2);
            f64_2x testRes1 = testFunc(inputVal, inputValB2);
            inputVal = inputVal * testRes1 + testRes0;
            f64_2x testRes2 = testFunc(inputVal, inputValB2);
            inputVal = inputVal * testRes2 + testRes1;
            f64_2x testRes3 = testFunc(inputVal, inputValB2);
            gRunSpeedSumF64F642x.md = _mm_add_pd(gRunSpeedSumF64F642x.md, (inputVal + testRes3).md);
        }
    }
    f32 seconds = linux_get_seconds_elapsed(start, linux_get_wall_clock());

    u32 tests = testsA * testsB;
    print_speed_info(name, maxNameSize, func, tests, seconds, gRunSpeedSumF64F642x.e[0], secondsBase, 4.0);
    return seconds;
}

