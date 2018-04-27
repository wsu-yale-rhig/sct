// sct/utils/random.cc

#include "sct/utils/random.hh"

#include "sct/core/base.hh"

namespace sct {
    
  Random::Random()
  : unit_uniform_(0.0, 1.0),
  two_unit_uniform_(0.0, 2.0),
  two_unit_centered_uniform_(-1.0, 1.0),
  zero_to_pi_(0.0, pi),
  counter_(0) {
    vector<double> x = {0.0, 1.0};
    vector<double> w = {0.0, 1.0};
    unit_linear_ = std::piecewise_linear_distribution<>(x.begin(), x.end(), w.begin());
    
    std::random_device rd;
    generator_.seed(rd());
  }
    
  Random::~Random() {
    
  }
    
  Random& Random::instance() {
    static Random instance;
    return instance;
  }
    
  double Random::uniform() {
    return unit_uniform_(generator_);
  }
    
  double Random::uniform2() {
    return two_unit_uniform_(generator_);
  }
    
  double Random::centeredUniform() {
    return two_unit_centered_uniform_(generator_);
  }
    
  double Random::zeroToPi() {
    return zero_to_pi_(generator_);
  }
  
  double Random::linear() {
    return unit_linear_(generator_);
  }
    
  void Random::seed(int seed) {
    if (seed < 0) {
      std::random_device rd;
      generator_.seed(rd());
    }
    else {
      generator_.seed(seed);
    }
  }
} // namespace sct
