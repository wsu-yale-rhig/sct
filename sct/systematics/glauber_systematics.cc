// sct/systematics/glauber_systematics.cc

#include "sct/systematics/glauber_systematics.hh"

#include "sct/core/logging.hh"

namespace sct {
  
  GlauberSystematics::GlauberSystematics(const string& out_dir) {
    out_dir_ = out_dir;
  }
  
  GlauberSystematics::~GlauberSystematics() {
    
  }
  
  bool GlauberSystematics::Load(const vector<string>& glauber_variations) {
    
    for (auto& file : glauber_variations) {
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
