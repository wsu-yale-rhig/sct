## get source files for libsct
file(GLOB tmp *.cc)
set(SCT_SRCS ${SCT_SRCS} ${tmp})

## get header files for libsct
file(GLOB tmp *.h)
set(SCT_HDRS ${SCT_HDRS} ${tmp})

## get test source files
file(GLOB tmp *_test.cc)
set(SCT_TEST_SRCS ${SCT_TEST_SRCS} ${tmp})

## and exclude the test files from libsct
if (SCT_TEST_SRCS)
list(REMOVE_ITEM SCT_SRCS ${SCT_TEST_SRCS})
endif(SCT_TEST_SRCS)

## get test benchmark files
file(GLOB tmp *_benchmark.cc)
set(SCT_BENCH_SRCS ${SCT_BENCH_SRCS} ${tmp})

## and exclude the benchmark files from libsct
if (SCT_BENCH_SRCS)
list(REMOVE_ITEM SCT_SRCS ${SCT_BENCH_SRCS})
endif(SCT_BENCH_SRCS)

## find gtest_main if it exists
file(GLOB tmp *gtest_main*)
set(SCT_TEST_MAIN ${SCT_TEST_MAIN} ${tmp})

## and exclude gtest_main from libsct
if (SCT_TEST_MAIN)
list(REMOVE_ITEM SCT_SRCS ${SCT_TEST_MAIN})
endif(SCT_TEST_MAIN)

## install headers
install(FILES ${SCT_HDRS} DESTINATION include/sct/stglauber)

## add subdirectories
## none

## export file lists to parent scope
set(SCT_SRCS ${SCT_SRCS} PARENT_SCOPE)
set(SCT_HDRS ${SCT_HDRS} PARENT_SCOPE)
set(SCT_TEST_SRCS ${SCT_TEST_SRCS} PARENT_SCOPE)
set(SCT_TEST_MAIN ${SCT_TEST_MAIN} PARENT_SCOPE)
set(SCT_BENCH_SRCS ${SCT_BENCH_SRCS} PARENT_SCOPE)
