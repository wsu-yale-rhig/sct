#include "sct/utils/random.h"

#include "sct/lib/math.h"

namespace sct {

Random::Random()
    : unit_uniform_(0.0, 1.0), two_unit_uniform_(0.0, 2.0),
      two_unit_centered_uniform_(-1.0, 1.0), zero_to_pi_(0.0, pi), counter_(0) {
  std::vector<double> x = {0.0, 1.0};
  std::vector<double> w = {0.0, 1.0};
  unit_linear_ =
      std::piecewise_linear_distribution<>(x.begin(), x.end(), w.begin());

  generator_.seed(54854);
}

Random::~Random() {}

Random &Random::instance() {
  static Random instance;
  return instance;
}

double Random::uniform() {
  rng_mutex_.lock();
  double ret = unit_uniform_(generator_);
  rng_mutex_.unlock();
  return ret;
}

double Random::uniform2() {
  rng_mutex_.lock();
  double ret = two_unit_uniform_(generator_);
  rng_mutex_.unlock();
  return ret;
}

double Random::centeredUniform() {
  rng_mutex_.lock();
  double ret = two_unit_centered_uniform_(generator_);
  rng_mutex_.unlock();
  return ret;
}

double Random::zeroToPi() {
  rng_mutex_.lock();
  double ret = zero_to_pi_(generator_);
  rng_mutex_.unlock();
  return ret;
}

double Random::linear() {
  rng_mutex_.lock();
  double ret = unit_linear_(generator_);
  rng_mutex_.unlock();
  return ret;
}

void Random::seed(int seed) {
  rng_mutex_.lock();
  if (seed < 0) {
    std::random_device rd;
    generator_.seed(rd());
  } else {
    generator_.seed(seed);
  }
  rng_mutex_.unlock();
}
} // namespace sct
