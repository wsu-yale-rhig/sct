#ifndef SCT_UTILS_NEGATIVE_BINOMIAL_H
#define SCT_UTILS_NEGATIVE_BINOMIAL_H

// implement histogram for NBD distribution

#include "sct/core/base.h"

#include "TH1D.h"

namespace sct {
  class NegativeBinomial {
  public:
  
    // creates a negative binomial distribution: 
    NegativeBinomial(double npp = 2.38, double k = 2.00);
  
    virtual ~NegativeBinomial();
  
    // returns two component multiplicity (1-x)* npart/2 + (x) * ncoll
    double twoComponentMultiplicity(double npart, double ncoll) const;
    // returns convolution with NBD
    double multiplicity(double npart, double ncoll) const;
  
    // randomly sample from the NBD
    unsigned random() const {return static_cast<unsigned>(nbd_->GetRandom());}
    
    // evaluate NBD(npp*m, k*m; n)
    double evaluateNBD(int i, double m = 1.0) const;
  
    void setParameters(double npp, double k);
  
    inline double npp() const {return npp_;}
    inline double k() const {return k_;}
    TH1D* distribution() const {return nbd_.get();}
  
  private:
    
    void initNBD();
    
    double npp_;                  // average multiplicity in pp
    double k_;                    // 1/k deviation from poisson
  
    unique_ptr<TH1D> nbd_;         // negative binomial distribution
  };
} // namespace sct

#endif // SCT_UTILS_NEGATIVE_BINOMIAL_H
