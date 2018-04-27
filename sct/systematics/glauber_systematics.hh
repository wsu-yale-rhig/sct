// sct/systematics/glauber_systematics.hh

#ifndef SCT_SYSTEMATICS_GLAUBER_SYSTEMATICS_HH
#define SCT_SYSTEMATICS_GLAUBER_SYSTEMATICS_HH

#include "sct/core/base.hh"
#include "sct/systematics/histogram_collection.hh"
#include "sct/centrality/multiplicity_model.hh"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"

namespace sct {
  
  class GlauberSystematics {
  public:
    
    GlauberSystematics();
    GlauberSystematics(const string& nominal, const vector<string>& modifiers);
    
    ~GlauberSystematics();
    
    bool Load(const string& nominal, const vector<string>& modifiers);
    
    void SetMultiplicityModel(double npp, double k, double x, double ppEff,
                              double aaEff, double aaCent, double trigEff = 1.0,
                              bool constEff = false);
    
  private:
    
    std::unique_ptr<TFile> nominal_;
    vector<std::unique_ptr<TFile>> variations_;
    
    std::unique_ptr<MultiplicityModel> mult_model_;
    
    
  };
  
} // namespace sct

#endif // SCT_SYSTEMATICS_GLAUBER_SYSTEMATICS_HH
