#!/bin/bash

set -e

curDir="$(pwd)"
codeDir="$curDir/src"
testDir="$curDir/test"
buildDir="$curDir/gebouw"

optsGcc="-O2 -std=gnu++17 -g -ggdb -march=x86-64 -msse4 -Wall -Werror -Wno-gnu-zero-variadic-macro-arguments -Wno-gnu-anonymous-struct -Wno-nested-anon-types -Wno-missing-braces -Wno-unused-function -Wno-writable-strings -Wno-c99-extensions -Wno-strict-aliasing -Wno-write-strings -Wno-unused-variable -Wno-unused-but-set-variable"

opts="-O2 -std=c++17 -g -ggdb -msse4 -Wall -Werror -Wno-gnu-zero-variadic-macro-arguments -Wno-gnu-anonymous-struct -Wno-nested-anon-types -Wno-missing-braces -Wno-unused-function -Wno-writable-strings -Wno-c99-extensions"

echo "Building std math"

mkdir -p "$buildDir"

pushd "$buildDir" > /dev/null
    clang++ $opts "$codeDir/simd_test.cpp" -o simd-test &

    g++ -DLIBBERDIP_EXPECT=0 $optsGcc "$testDir/test_math.cpp" -o test-math-gcc &
    clang++ -DLIBBERDIP_EXPECT=0 $opts "$testDir/test_math.cpp" -o test-math-clang &

    clang++ -DLIBBERDIP_EXPECT=0 $opts "$testDir/test_mats_rounding.cpp" -o test-mats-rounding &
    clang++ -DLIBBERDIP_EXPECT=0 $opts "$testDir/test_mats_elem.cpp"     -o test-mats-elem &
    clang++ -DLIBBERDIP_EXPECT=0 $opts "$testDir/test_mats_trig.cpp"     -o test-mats-trig &
    g++ -DLIBBERDIP_EXPECT=0 $optsGcc  "$testDir/test_mats_trig.cpp"     -o test-mats-trig-gcc &
    g++ -DLIBBERDIP_EXPECT=0 $optsGcc  "$testDir/test_mats_elem.cpp"     -o test-mats-elem-gcc &

    clang++ -DLIBBERDIP_EXPECT=0 $opts "$testDir/test_specials.cpp"      -o test-mats-special &

    clang++ -DLIBBERDIP_EXPECT=0 $opts "$testDir/test_fft.cpp" -o test-mats-fft &
    g++ -DLIBBERDIP_EXPECT=0 $optsGcc  "$testDir/test_fft.cpp" -o test-mats-fft-gcc &
popd > /dev/null

wait
