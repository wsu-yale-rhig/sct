// sct/utils/random.hh

#ifndef SCT_UTILS_RANDOM_H
#define SCT_UTILS_RANDOM_H


// implements all of the random number distributions we need
// to sample during MC glauber generation, as well as a globally
// unique counter for root object naming.

// to use the generator inside the sct library, call the singleton
// Random::instance().

#include <random>
#include <atomic>

namespace sct {
    
  // singleton Random implementation built
  // on the standard library random
  class Random {
  public:
    
    static Random& instance();
    virtual ~Random();
    
    // samples [0, 1] uniformly
    double uniform();
      
    // samples [0, 2] uniformly
    double uniform2();
    
    // samples [-1, 1] uniformly
    double centeredUniform();
      
    // samples [0, pi] uniformly
    double zeroToPi();
      
    // samples [0, 1] with a linearly increasing probability distribution
    // proportional to x
    double linear();
    
    // returns a globally unique counter, useful for histogram
    // naming where identical names can cause problems... thanks ROOT
    inline unsigned counter() {return counter_++;}
      
    // can reset the seed for the random number engine - does not reset
    // counter
    void seed(int seed = -1);
      
  private:
    
    std::mt19937 generator_;
    
    std::uniform_real_distribution<> unit_uniform_;
    std::uniform_real_distribution<> two_unit_uniform_;
    std::uniform_real_distribution<> two_unit_centered_uniform_;
    std::uniform_real_distribution<> zero_to_pi_;
    
    std::piecewise_linear_distribution<> unit_linear_;
    
    std::atomic<unsigned> counter_;
    
    Random();
    Random(const Random&);
      
  };
    
} // namespace sct

#endif // SCT_UTILS_RANDOM_H
