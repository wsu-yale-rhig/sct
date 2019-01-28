#include "sct/centrality/multiplicity_model.h"

#include "sct/lib/logging.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/random.h"

#include <limits>

#include "TMath.h"

namespace sct {

MultiplicityModel::MultiplicityModel(double npp, double k, double x,
                                     double ppEff, double centEff,
                                     double centMult, double triggerBias,
                                     bool const_eff)
    : pp_efficiency_(ppEff),
      central_efficiency_(centEff),
      cent_mult_(centMult),
      trigger_bias_(triggerBias),
      x_(x),
      const_efficiency_(const_eff) {
  setNBD(npp, k);
}

MultiplicityModel::MultiplicityModel(const MultiplicityModel& rhs)
    : pp_efficiency_(rhs.pp_efficiency_),
      central_efficiency_(rhs.central_efficiency_),
      cent_mult_(rhs.cent_mult_),
      trigger_bias_(rhs.trigger_bias_),
      x_(rhs.x_),
      const_efficiency_(rhs.const_efficiency_) {
  setNBD(rhs.npp(), rhs.k());
}

MultiplicityModel::~MultiplicityModel() {}

double MultiplicityModel::twoComponentMultiplicity(double npart,
                                                   double ncoll) const {
  return (x_ == 0) ? npart : ((1 - x_) * npart / 2.0) + (x_ * ncoll);
}

double MultiplicityModel::multiplicity(double npart, double ncoll) const {
  // taking into account trigger & TPC efficiency,
  // get multiplicity from NBD

  double nch_pp = twoComponentMultiplicity(npart, ncoll);
  double nch_sampled = nch_pp;

  // do the sampling
  unsigned ideal_mult = 0;
  for (int i = 0; i < TMath::Nint(nch_sampled); ++i) {
    ideal_mult += random();
  }

  // get the efficiency and modify if multiplicity dependent
  double eff = evalEfficiency(ideal_mult);

  unsigned mult = 0;
  for (int i = 0; i < ideal_mult; ++i) {
    if (Random::instance().uniform() < eff) mult++;
  }

  if (trigger_bias_ == 1.0) return mult;

  int count = mult;
  for (int i = 0; i < count; ++i) {
    if (Random::instance().uniform() < trigger_bias_) mult++;
  }
  return mult;
}

TH1D* MultiplicityModel::multiplicity(double npart, double ncoll,
                                      double weight) const {
  double nch_pp = twoComponentMultiplicity(npart, ncoll);

  // include the trigger bias
  int nch = TMath::Nint(nch_pp * trigger_bias_);

  // get the efficiency and modify if multiplicity dependent
  double eff = evalEfficiency(nch);
  int n_sampled = TMath::Nint(nch * eff);

  unsigned n_bins = 1000;
  TH1D* h =
      new TH1D(MakeString("mult_tmp_", Random::instance().counter()).c_str(),
               "", n_bins, 0, n_bins);

  for (int i = 0; i < n_bins; ++i) {
    double prob = evaluateNBD(i, n_sampled);

    if (prob > 0.0 && prob < std::numeric_limits<double>::max()) {
      h->Fill(i + 0.5, prob);
    }
  }

  return h;
}

double MultiplicityModel::evalEfficiency(unsigned mult) const {
  if (const_efficiency_) return central_efficiency_;

  double d = pp_efficiency_ - central_efficiency_;
  return pp_efficiency_ * (1.0 - mult * d / cent_mult_);
}

void MultiplicityModel::setNBD(double npp, double k) { setParameters(npp, k); }

}  // namespace sct
