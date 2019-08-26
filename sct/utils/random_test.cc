#include "sct/utils/random.h"
#include "sct/lib/logging.h"

#include "gtest/gtest.h"

#include <future>
#include <stdlib.h>
#include <thread>
#include <vector>

// tests the boundaries and the averages of the
// distributions implemented in random.h

TEST(unit_random, bounds) {
  sct::Random &random = sct::Random::instance();
  unsigned n_throws = 1e6;
  for (int i = 0; i < n_throws; ++i) {
    double tmp = random.uniform();
    EXPECT_LE(0.0, tmp);
    EXPECT_GE(1.0, tmp);
  }
}

TEST(unit_random, average) {
  sct::Random &random = sct::Random::instance();
  unsigned n_throws = 1e6;
  double counter = 0.0;
  for (int i = 0; i < n_throws; ++i) {
    counter += random.uniform();
  }
  EXPECT_NEAR(0.5, counter / n_throws, 1e-3);
}

TEST(two_unit_random, bounds) {
  sct::Random &random = sct::Random::instance();
  unsigned n_throws = 1e6;
  for (int i = 0; i < n_throws; ++i) {
    double tmp = random.uniform2();
    EXPECT_LE(0.0, tmp);
    EXPECT_GE(2.0, tmp);
  }
}

TEST(two_unit_random, average) {
  sct::Random &random = sct::Random::instance();
  unsigned n_throws = 5e6;
  double counter = 0.0;
  for (int i = 0; i < n_throws; ++i) {
    counter += random.uniform2();
  }
  EXPECT_NEAR(1.0, counter / n_throws, 1e-3);
}

TEST(two_unit_center_random, bounds) {
  sct::Random &random = sct::Random::instance();
  unsigned n_throws = 1e6;
  for (int i = 0; i < n_throws; ++i) {
    double tmp = random.centeredUniform();
    EXPECT_LE(-1.0, tmp);
    EXPECT_GE(1.0, tmp);
  }
}

TEST(two_unit_center_random, average) {
  sct::Random &random = sct::Random::instance();
  unsigned n_throws = 5e6;
  double counter = 0.0;
  for (int i = 0; i < n_throws; ++i) {
    counter += random.centeredUniform();
  }
  EXPECT_NEAR(0.0, counter / n_throws, 1e-3);
}

TEST(linear_random, bounds) {
  sct::Random &random = sct::Random::instance();
  unsigned n_throws = 1e5;
  for (int i = 0; i < n_throws; ++i) {
    double tmp = random.linear();
    EXPECT_LE(0.0, tmp);
    EXPECT_GE(1.0, tmp);
  }
}

TEST(linear_random, average) {
  sct::Random &random = sct::Random::instance();
  unsigned n_throws = 1e6;
  double counter = 0.0;
  for (int i = 0; i < n_throws; ++i) {
    counter += random.linear();
  }
  EXPECT_NEAR(2.0 / 3.0, counter / n_throws, 1e-3);
}

double seed_thread_test(int seed) {

  sct::Random::instance().seed(seed);
  sleep(4);
  return sct::Random::instance().uniform();
}

TEST(random, threaded) {
  std::vector<double> results;
  auto future1 = std::async(seed_thread_test, 1);
  sleep(1);
  auto future2 = std::async(seed_thread_test, 1);
  sleep(1);
  auto future3 = std::async(seed_thread_test, 1);
  auto future4 = std::async(seed_thread_test, 1323);
  double result1 = future1.get();
  double result2 = future2.get();
  double result3 = future3.get();
  double result4 = future4.get();
  EXPECT_EQ(result1, result2);
  EXPECT_EQ(result1, result3);
  EXPECT_NE(result3, result4);
}