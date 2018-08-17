// sct/systematics/systematic_variable.hh

#ifndef SCT_SYSTEMATICS_SYSTEMATIC_VARIABLE_HH
#define SCT_SYSTEMATICS_SYSTEMATIC_VARIABLE_HH

#include "sct/core/base.hh"
#include "sct/systematics/histogram_collection.hh"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

namespace sct {
  
  class SystematicVariable {
  public:
    
    SystematicVariable() {}
    
    ~SystematicVariable() {}
    
    virtual void Fill1D() {}
    virtual void Fill2D() {}
    virtual void FillProfile() {}
    
  private:
    
    std::string name_; // variable name (impact_parameter, etc)
    std::string title_; // histogram title
    std::string x_title_; // x axis title (if empty, uses name_)
    std::string y_title_; // y axis title
    
    // histogram collections
    HistogramCollection<TH1D> th1_;
    HistogramCollection<TH2D> th2_;
    HistogramCollection<TProfile> tprof_;
    HistogramCollection<TProfile> weights_;
  };
  
}

#endif // SCT_SYSTEMATICS_SYSTEMATIC_VARIABLE_HH
