# get binary source files
file(GLOB tmp *.cc)
set(SCT_BINARY_SRCS ${SCT_BINARY_SRCS} ${tmp})

# get test source files
file(GLOB tmp *_test.cc)
set(SCT_TEST_SRCS ${SCT_TEST_SRCS} ${tmp})

# and exclude the test files from analysis routines
if(SCT_TEST_SRCS)
  list(REMOVE_ITEM SCT_BINARY_SRCS ${SCT_TEST_SRCS})
endif(SCT_TEST_SRCS)

# get test benchmark files
file(GLOB tmp *_benchmark.cc)
set(SCT_BENCH_SRCS ${SCT_BENCH_SRCS} ${tmp})

# and exclude the benchmark files from analysis routines
if(SCT_BENCH_SRCS)
  list(REMOVE_ITEM SCT_BINARY_SRCS ${SCT_BENCH_SRCS})
endif(SCT_BENCH_SRCS)

# export file lists to parent scope
set(SCT_BINARY_SRCS ${SCT_BINARY_SRCS} PARENT_SCOPE)
set(SCT_TEST_SRCS ${SCT_TEST_SRCS} PARENT_SCOPE)
set(SCT_BENCH_SRCS ${SCT_BENCH_SRCS} PARENT_SCOPE)
