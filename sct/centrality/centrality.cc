#include "sct/centrality/centrality.h"

#include "sct/lib/logging.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/random.h"

#include "TF1.h"

namespace sct {

// helper function to try and avoid rounding errors, when doing
// the fractional calculations in the centrality definitions
template <typename T>
T Round(T t, int digits) {
  if (t == 0.0)  // otherwise it will return 'nan' due to the log10() of zero
    return 0.0;

  double factor = pow(10.0, digits - ceil(log10(fabs(t))));
  return round(t * factor) / factor;
}

Centrality::Centrality(TH1D* data, TH1D* simu) {
  if (data != nullptr) setDataRefmult(data);
  if (simu != nullptr) setSimuRefmult(simu);
}

Centrality::~Centrality() {}

void Centrality::setDataRefmult(TH1D* histogram) {
  data_ = make_unique<TH1D>(*histogram);
  data_->SetName(MakeString("centrality_internal_data_refmult_",
                            Random::instance().counter())
                     .c_str());
  data_->SetDirectory(0);
}

void Centrality::setSimuRefmult(TH1D* histogram) {
  simu_ = make_unique<TH1D>(*histogram);
  simu_->SetName(MakeString("centrality_internal_simu_refmult_",
                            Random::instance().counter())
                     .c_str());
  simu_->SetDirectory(0);
}

std::vector<unsigned> Centrality::centralityBins(XSecMod mod) {
  if (simu_ == nullptr) {
    LOG(ERROR) << "No simulation refmult distribution set";
    LOG(ERROR) << "Centrality calculation failed";
    return std::vector<unsigned>();
  }

  // calculate centrality from MC multiplicity distribution

  // Perform the calculation twice, first from 0 -> 100%, then from 100 -> 0%
  // as a sanity check
  std::vector<unsigned> bounds_forward =
      integrate(simu_.get(), Integral::Forward, mod);
  std::vector<unsigned> bounds_backward =
      integrate(simu_.get(), Integral::Backward, mod);

  if (bounds_forward != bounds_backward) {
    LOG(ERROR)
        << "Centrality bounds not equal when integrating forwards & backwards";
    LOG(ERROR) << "make sure the bin width of the refmult histogram is 1.";
    LOG(ERROR) << "Centrality calculation failed";
    string mod_string;
    switch (mod) {
      case XSecMod::None:
        mod_string = "None";
        break;
      case XSecMod::Plus5:
        mod_string = "+5%";
        break;
      case XSecMod::Minus5:
        mod_string = "-5%";
        break;
    }
    LOG(ERROR) << "xsec mod: " << mod_string;
    LOG(ERROR) << "bounds forward: " << bounds_forward;
    LOG(ERROR) << "bounds backward: " << bounds_backward;
    return std::vector<unsigned>();
  }

  return bounds_forward;
}

std::vector<unsigned> Centrality::integrate(TH1D* h, Integral direction,
                                            XSecMod mod) {
  // returns a std::vector of multiplicity boundaries
  std::vector<unsigned> boundaries;

  // first get xsec modification
  double xsec_percent = 1.0;
  switch (mod) {
    case XSecMod::Plus5:
      xsec_percent = 1.05;
      break;
    case XSecMod::Minus5:
      xsec_percent = 0.95;
      break;
    case XSecMod::None:
      xsec_percent = 1.0;
      break;
  }

  // copy our centrality bounds from sct/lib/enumeration.h
  std::vector<double> cent_min = Centrality16BinLowerBound;
  std::vector<double> cent_max = Centrality16BinUpperBound;

  // get the number of cuts
  unsigned n_cent_bins = cent_min.size();
  boundaries.resize(n_cent_bins);

  // now define our boundary cuts depending on if we are integrating forward or
  // backward, and weighted by the xsec modification
  std::vector<double> cent_cuts;
  switch (direction) {
    case Integral::Forward:
      for (int i = 0; i < n_cent_bins; ++i)
        cent_cuts.push_back(cent_max[i] / 100.0 * xsec_percent);
      break;
    case Integral::Backward:
      for (int i = 0; i < n_cent_bins; ++i)
        cent_cuts.push_back(cent_max[n_cent_bins - 1 - i] / 100.0 *
                            xsec_percent);
      break;
  }
  std::reverse(cent_cuts.begin(), cent_cuts.end());

  // now do the integration - goal is to integrate until we find the bin such
  // that the cumulative integral up to that point, normalized by the total
  // integral, passes one of the centrality bin cuts
  unsigned cent_bin = 0;
  unsigned nbins_hist = h->GetNbinsX();
  double norm = h->Integral();

  for (int i = 1; i <= nbins_hist; ++i) {
    // define the bins to integrate over, depending on if we
    // integrate forward or backwards, and the multiplicity
    unsigned bin_low, bin_high, mult;
    switch (direction) {
      case Integral::Forward:
        bin_low = 1;
        bin_high = i;
        mult = bin_high - 1;
        break;
      case Integral::Backward:
        bin_low = nbins_hist + 1 - i;
        bin_high = nbins_hist;
        mult = bin_low - 1;
        break;
    }

    // get the integral
    double integral = h->Integral(bin_low, bin_high);

    // and the fraction of the total integral
    double fraction = integral / norm;

    // we need to round the fraction to avoid some rounding errors in the
    // floating point arithmetic, which of course involves more errors...
    fraction = Round(fraction, 12);

    // now loop to find if it sits on a bin edge
    switch (direction) {
      case Integral::Forward:
        for (int bin = cent_bin; bin < n_cent_bins; ++bin) {
          if (1.0 - fraction < cent_cuts[bin]) {
            boundaries[cent_bin++] = mult;
          }
        }
        break;
      case Integral::Backward:
        for (int bin = cent_bin; bin < n_cent_bins; ++bin) {
          if (fraction >= cent_cuts[bin]) {
            boundaries[cent_bin++] = mult;
          }
        }
        break;
    }
  }

  // if we are doing backwards integration, we have to sort the entries
  if (direction == Integral::Backward)
    std::sort(boundaries.begin(), boundaries.end());

  return boundaries;
}

std::pair<std::vector<double>, unique_ptr<TH1D>> Centrality::weights(
    unsigned fit_boundary) {
  // make sure both histograms exist
  if (data_.get() == nullptr) {
    LOG(ERROR)
        << "no data refmult distribution loaded, can't calculate weights";
    return {std::vector<double>(), unique_ptr<TH1D>()};
  }

  if (simu_.get() == nullptr) {
    LOG(ERROR)
        << "no simulated refmult distribution loaded, can't calculate weights";
    return {std::vector<double>(), unique_ptr<TH1D>()};
  }

  // now make a copy of the simulated refmult and take the ratio
  unique_ptr<TH1D> ratio = make_unique<TH1D>(*simu_.get());
  ratio->Divide(data_.get());

  // define the fit function
  string fxn_def =
      "[0] + [1]/([2]*x + [3]) + [4]*([2]*x + [3]) + [5]/([2]*x + [3])^2 + "
      "[6]*([2]*x + [3])^2";
  unique_ptr<TF1> ratio_fit =
      make_unique<TF1>("ratio_fit", fxn_def.c_str(), 0, fit_boundary);

  // define some reasonable defaults to help the fit
  ratio_fit->SetParameters(1.5, -20, 1, 6, 1e-3, 600, 5e-6);

  // and fit
  ratio->Fit(ratio_fit.get(), "EMR");
  std::vector<double> weights(7);
  for (int par = 0; par < ratio_fit->GetNpar(); ++par)
    weights[par] = ratio_fit->GetParameter(par);

  return {weights, std::move(ratio)};
}

}  // namespace sct
