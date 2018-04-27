// sct/centrality/multiplicity_model.cc

#include "sct/centrality/multiplicity_model.hh"

#include <limits>

#include "sct/utils/random.hh"
#include "sct/core/logging.hh"

#include "TMath.h"

namespace sct {
  
  MultiplicityModel::MultiplicityModel(double npp, double k,
                                       double x, double ppEff,
                                       double centEff, double centMult,
                                       double triggerBias, bool const_eff)
  : ppEfficiency_(ppEff), centralEfficiency_(centEff), centMult_(centMult),
  triggerBias_(triggerBias), x_(x), constEfficiency_(const_eff) {
  setNBD(npp, k);
  }
  
  MultiplicityModel::MultiplicityModel(const MultiplicityModel& rhs)
  : ppEfficiency_(rhs.ppEfficiency_), centralEfficiency_(rhs.centralEfficiency_),
  centMult_(rhs.centMult_), triggerBias_(rhs.triggerBias_), x_(rhs.x_),
  constEfficiency_(rhs.constEfficiency_) {
    setNBD(rhs.npp(), rhs.k());
  }
  
  MultiplicityModel::~MultiplicityModel() { }
  
  double MultiplicityModel::twoComponentMultiplicity(double npart,
                                                    double ncoll) const {
    return (x_ == 0) ? npart : ((1 - x_) * npart / 2.0) + (x_ * ncoll);
  }
  
  double MultiplicityModel::multiplicity(double npart, double ncoll) const {
    // taking into account trigger & TPC efficiency,
    // get multiplicity from NBD
    
    double nchPP = twoComponentMultiplicity(npart, ncoll);
    double nChSampled = nchPP;
    
    // do the sampling
    unsigned ideal_mult = 0;
    for (int i = 0; i < TMath::Nint(nChSampled); ++i) {
      ideal_mult += random();
    }
    
    // get the efficiency and modify if multiplicity dependent
    double eff = evalEfficiency(ideal_mult);
    
    unsigned mult = 0;
    for (int i = 0; i < ideal_mult; ++i) {
      if (Random::instance().uniform() < eff)
      mult++;
    }
    
    if (triggerBias_ == 1.0)
    return mult;
    
    int count = mult;
    for (int i = 0; i < count; ++i) {
      if (Random::instance().uniform() < triggerBias_)
      mult++;
    }
    return mult;
  }
  
  TH1D* MultiplicityModel::multiplicity(double npart, double ncoll,
                                        double weight) const {
    double nchPP = twoComponentMultiplicity(npart, ncoll);
    
    // include the trigger bias
    int nch = TMath::Nint(nchPP * triggerBias_);
    
    // get the efficiency and modify if multiplicity dependent
    double eff = evalEfficiency(nch);
    int nSampled = TMath::Nint(nch * eff);
    
    unsigned nBins = 1000;
    TH1D* h = new TH1D(MakeString("mult_tmp_", Random::instance().counter()).c_str(),
                                  "", nBins, 0, nBins);
    
    for (int i = 0; i < nBins; ++i) {
      double prob = evaluateNBD(i, nSampled);
      
      if (prob > 0.0 && prob < std::numeric_limits<double>::max()) {
        h->Fill(i+0.5, prob);
      }
    }
    
    return h;
  }
  
  double MultiplicityModel::evalEfficiency(unsigned mult) const {
    if (constEfficiency_)
    return centralEfficiency_;
    
    double d = ppEfficiency_ - centralEfficiency_;
    return ppEfficiency_ * (1.0 - mult * d / centMult_);
  }
  
  void MultiplicityModel::setNBD(double npp, double k) {
    setParameters(npp, k);
  }
  
} // namespace sct
