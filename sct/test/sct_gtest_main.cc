// Entrypoint to run all tests registered to gtest
// for more information on googletest, see 
// https://github.com/google/googletest

#include "gtest/gtest.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}