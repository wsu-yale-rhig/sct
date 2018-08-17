// sct/systematics/glauber_systematics.hh

#ifndef SCT_SYSTEMATICS_GLAUBER_SYSTEMATICS_HH
#define SCT_SYSTEMATICS_GLAUBER_SYSTEMATICS_HH

#include "sct/core/base.hh"
#include "sct/centrality/multiplicity_model.hh"
#include "sct/systematics/systematic_cumulant.hh"
#include "sct/systematics/systematic_variable.hh"

#include <array>

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFile.h"

namespace sct {
  
  class GlauberSystematics {
  public:
    
    GlauberSystematics(const string& out_dir = "./");
    
    ~GlauberSystematics();

    void SetOutputDirectory(const string& out_dir = "./") {out_dir_ = out_dir;}
    
    // specify all the input files - each file should be a separate
    // variation of the glauber model
    bool Load(const vector<string>& glauber_variations);
    
    // need to set the multiplicity model parameters to the best-fit to data
    void SetMultiplicityModel(double npp, double k, double x, double ppEff,
                              double aaEff, double aaCent, double trigEff = 1.0,
                              bool constEff = false);
    
    // runs systematic analysis over all inputs
    bool Run();
    
    // writes results to disk in directory out_dir_
    bool Write();

  private:
    
    bool Initialize();
    
    string out_dir_;
    
    vector<unique_ptr<TFile>> variations_;
    
    unique_ptr<MultiplicityModel> mult_model_;
    
    // histogram containers
    SystematicVariable impact_parameter_;
    SystematicVariable n_part_;
    SystematicVariable n_coll_;
    SystematicVariable multiplicity_;
    SystematicVariable area_rp_;
    SystematicVariable area_pp_;
    SystematicCumulant ecc_rp_; // npart weight
    SystematicCumulant ecc_rp_mult_; // multiplicity weight
    std::array<SystematicCumulant, 3> ecc_pp_; // npart weight
    std::array<SystematicCumulant, 3> ecc_pp_mult_; // multiplicity weight
    
  };
  
} // namespace sct

#endif // SCT_SYSTEMATICS_GLAUBER_SYSTEMATICS_HH
