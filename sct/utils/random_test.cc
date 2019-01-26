#include "gtest/gtest.h"
#include "sct/utils/random.h"

// tests the boundaries and the averages of the
// distributions implemented in random.h

TEST(unit_random, bounds) {
  sct::Random& random = sct::Random::instance();
  unsigned n_throws = 1e5;
  for (int i = 0; i < n_throws; ++i) {
    double tmp = random.uniform();
    EXPECT_LE(0.0, tmp);
    EXPECT_GE(1.0, tmp);
  }
}

TEST(unit_random, average) {
  sct::Random& random = sct::Random::instance();
  unsigned n_throws = 1e6;
  double counter = 0.0;
  for (int i = 0; i < n_throws; ++i) {
    counter += random.uniform();
  }
  EXPECT_NEAR(0.5, counter / n_throws, 1e-3);
}

TEST(two_unit_random, bounds) {
  sct::Random& random = sct::Random::instance();
  unsigned n_throws = 1e5;
  for (int i = 0; i < n_throws; ++i) {
    double tmp = random.uniform2();
    EXPECT_LE(0.0, tmp);
    EXPECT_GE(2.0, tmp);
  }
}

TEST(two_unit_random, average) {
  sct::Random& random = sct::Random::instance();
  unsigned n_throws = 1e6;
  double counter = 0.0;
  for (int i = 0; i < n_throws; ++i) {
    counter += random.uniform2();
  }
  EXPECT_NEAR(1.0, counter / n_throws, 1e-3);
}

TEST(two_unit_center_random, bounds) {
  sct::Random& random = sct::Random::instance();
  unsigned n_throws = 1e5;
  for (int i = 0; i < n_throws; ++i) {
    double tmp = random.centeredUniform();
    EXPECT_LE(-1.0, tmp);
    EXPECT_GE(1.0, tmp);
  }
}

TEST(two_unit_center_random, average) {
  sct::Random& random = sct::Random::instance();
  unsigned n_throws = 1e6;
  double counter = 0.0;
  for (int i = 0; i < n_throws; ++i) {
    counter += random.centeredUniform();
  }
  EXPECT_NEAR(0.0, counter / n_throws, 1e-3);
}

TEST(linear_random, bounds) {
  sct::Random& random = sct::Random::instance();
  unsigned n_throws = 1e5;
  for (int i = 0; i < n_throws; ++i) {
    double tmp = random.linear();
    EXPECT_LE(0.0, tmp);
    EXPECT_GE(1.0, tmp);
  }
}

TEST(linear_random, average) {
  sct::Random& random = sct::Random::instance();
  unsigned n_throws = 1e6;
  double counter = 0.0;
  for (int i = 0; i < n_throws; ++i) {
    counter += random.linear();
  }
  EXPECT_NEAR(2.0/3.0, counter / n_throws, 1e-3);
}
