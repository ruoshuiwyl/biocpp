cmake_minimum_required(VERSION 3.6)
project(src)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_COMPILER /usr/bin/g++-5)
set(CMAKE_BUILD_TYPE Debug)

set(SOURCE_FILES
    pairhmm/common_pairHMM.cc
    pairhmm/common_pairHMM.h
    pairhmm/Log10_PairHMM.cc
    pairhmm/Log10_PairHMM.h
    pairhmm/pairHMM.cc
    pairhmm/pairHMM.h
    util/read.h
    util/util.cc
    util/util.h
    sw_align/sw_align.cc
    sw_align/sw_align.h)

set(TEST_FILES
    pairhmm/pairHMM_test.cc
    main.cc sw_align/sw_align_test.cc)


add_executable(pairhmm ${SOURCE_FILES} ${TEST_FILES})
target_link_libraries(pairhmm gtest pthread)
