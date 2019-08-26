#include "sct/utils/random.h"

#include "sct/lib/math.h"

namespace sct {

Counter::Counter() : counter_(0) {}

Counter::~Counter() {}

Counter &Counter::instance() {
  static Counter instance;
  return instance;
}

Random::Random()
    : unit_uniform_(0.0, 1.0), two_unit_uniform_(0.0, 2.0),
      two_unit_centered_uniform_(-1.0, 1.0), zero_to_pi_(0.0, pi) {
  std::vector<double> x = {0.0, 1.0};
  std::vector<double> w = {0.0, 1.0};
  unit_linear_ =
      std::piecewise_linear_distribution<>(x.begin(), x.end(), w.begin());

  generator_.seed(54854);
}

Random &Random::instance() {
  thread_local static Random instance;
  return instance;
}

Random::~Random() {}

double Random::uniform() {
  double ret = unit_uniform_(generator_);
  return ret;
}

double Random::uniform2() {
  double ret = two_unit_uniform_(generator_);
  return ret;
}

double Random::centeredUniform() {
  double ret = two_unit_centered_uniform_(generator_);
  return ret;
}

double Random::zeroToPi() {
  double ret = zero_to_pi_(generator_);
  return ret;
}

double Random::linear() {
  double ret = unit_linear_(generator_);
  return ret;
}

void Random::seed(int seed) {
  if (seed < 0) {
    std::random_device rd;
    generator_.seed(rd());
  } else {
    generator_.seed(seed);
  }
}
} // namespace sct
