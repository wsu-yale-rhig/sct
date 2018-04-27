// sct/systematics/glauber_systematics.cc

#include "sct/systematics/glauber_systematics.hh"

#include "sct/core/logging.hh"

namespace sct {
  
  GlauberSystematics::GlauberSystematics() {
    
  }
  
  GlauberSystematics::GlauberSystematics(const string& nominal, const vector<string>& modifiers) {
    Load(nominal, modifiers);
  }
  
  GlauberSystematics::~GlauberSystematics() {
    
  }
  
  bool GlauberSystematics::Load(const string& nominal, const vector<string>& modifiers) {
    nominal_ = make_unique<TFile>(nominal.c_str(), "READ");
    
    if (nominal_.get() == nullptr || !nominal_->IsOpen()) {
      LOG(ERROR) << "Could not open default glauber file: " << nominal << ": file loading failed";
      return false;
    }
    
    for (auto& file : modifiers) {
      variations_.push_back(make_unique<TFile>(file.c_str(), "READ"));
      if (variations_.back() == nullptr || !variations_.back()->IsOpen()) {
        LOG(ERROR) << "Could not open glauber file: " << file << ": file loading failed";
        variations_.clear();
        return false;
      }
    }
    return true;
  }
  
  void GlauberSystematics::SetMultiplicityModel(double npp, double k, double x, double ppEff,
                                                double aaEff, double aaCent, double trigEff,
                                                bool constEff) {
    mult_model_ = make_unique<MultiplicityModel>(npp, k, x, ppEff, aaEff, aaCent, trigEff, constEff);
  }
  
} // namespace sct
