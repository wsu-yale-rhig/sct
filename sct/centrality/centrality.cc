// sct/centrality/centrality.cc

#include "sct/centrality/centrality.hh"

#include "sct/utils/random.hh"
#include "sct/core/logging.hh"

#include "TF1.h"

namespace sct {
  
  // helper function to try and avoid rounding errors, when doing
  // the fractional calculations in the centrality definitions
  template<typename T>
  T Round(T t, int digits) {
    if (t == 0.0) // otherwise it will return 'nan' due to the log10() of zero
      return 0.0;
    
    double factor = pow(10.0, digits - ceil(log10(fabs(t))));
    return round(t * factor) / factor;
  }
  
  Centrality::Centrality(TH1D* data, TH1D* simu) {
    if (data != nullptr)
      setDataRefmult(data);
    if (simu != nullptr)
      setSimuRefmult(simu);
  }
  
  Centrality::~Centrality() {
    
  }
  
  void Centrality::setDataRefmult(TH1D* histogram) {
    data_ = make_unique<TH1D>(*histogram);
    data_->SetName(MakeString("centrality_internal_data_refmult_",
                             Random::instance().counter()).c_str());
    data_->SetDirectory(0);
  }
  
  void Centrality::setSimuRefmult(TH1D* histogram) {
    simu_ = make_unique<TH1D>(*histogram);
    simu_->SetName(MakeString("centrality_internal_simu_refmult_",
                             Random::instance().counter()).c_str());
    simu_->SetDirectory(0);
  }
  
  vector<unsigned> Centrality::centralityBins(XSecMod mod) {
    if (simu_ == nullptr) {
      LOG(ERROR) << "No simulation refmult distribution set";
      LOG(ERROR) << "Centrality calculation failed";
      return vector<unsigned>();
    }
    
    // calculate centrality from MC multiplicity distribution
    
    // Perform the calculation twice, first from 0 -> 100%, then from 100 -> 0%
    // as a sanity check
    vector<unsigned> bounds_forward = integrate(simu_.get(), Integral::Forward, mod);
    vector<unsigned> bounds_backward = integrate(simu_.get(), Integral::Backward, mod);
    
    if (bounds_forward != bounds_backward) {
      LOG(ERROR) << "Centrality bounds not equal when integrating forwards & backwards";
      LOG(ERROR) << "make sure the bin width of the refmult histogram is 1.";
      LOG(ERROR) << "Centrality calculation failed";
      string mod_string;
      switch (mod) {
        case XSecMod::None :
          mod_string = "None";
          break;
        case XSecMod::Plus5 :
          mod_string = "+5%";
          break;
        case XSecMod::Minus5 :
          mod_string = "-5%";
          break;
          
      }
      LOG(ERROR) << "xsec mod: " << mod_string;
      LOG(ERROR) << "bounds forward: " << bounds_forward;
      LOG(ERROR) << "bounds backward: " << bounds_backward;
      return vector<unsigned>();
    }
    
    return bounds_forward;
  }
  
  vector<unsigned> Centrality::integrate(TH1D* h, Integral direction, XSecMod mod) {
    // returns a vector of multiplicity boundaries
    vector<unsigned> boundaries;
    
    // first get xsec modification
    double xSecPercent = 1.0;
    switch (mod) {
      case XSecMod::Plus5 :
        xSecPercent = 1.05;
        break;
      case XSecMod::Minus5 :
        xSecPercent = 0.95;
        break;
      case XSecMod::None :
        xSecPercent = 1.0;
        break;
    }
    
    // copy our centrality bounds from sct/core/base.hh
    vector<double> centMin = Centrality16BinLowerBound;
    vector<double> centMax = Centrality16BinUpperBound;
    
    // get the number of cuts
    unsigned nCentBins = centMin.size();
    boundaries.resize(nCentBins);
    
    // now define our boundary cuts depending on if we are integrating forward or backward,
    // and weighted by the xsec modification
    vector<double> centCuts;
    switch (direction) {
      case Integral::Forward :
        for (int i = 0; i < nCentBins; ++i)
          centCuts.push_back(centMax[i] / 100.0 * xSecPercent);
        break;
      case Integral::Backward :
        for (int i = 0; i < nCentBins; ++i)
          centCuts.push_back(centMax[nCentBins-1-i] / 100.0 * xSecPercent);
        break;
    }
    std::reverse(centCuts.begin(), centCuts.end());
    
    // now do the integration - goal is to integrate until we find the bin such that the
    // cumulative integral up to that point, normalized by the total integral, passes one
    // of the centrality bin cuts
    unsigned centBin = 0;
    unsigned nBinsHist = h->GetNbinsX();
    double norm = h->Integral();
  
    for (int i = 1; i <= nBinsHist; ++i) {
      // define the bins to integrate over, depending on if we
      // integrate forward or backwards, and the multiplicity
      unsigned binLow, binHigh, mult;
      switch (direction) {
        case Integral::Forward :
          binLow = 1;
          binHigh = i;
          mult = binHigh - 1;
          break;
        case Integral::Backward :
          binLow = nBinsHist + 1 - i;
          binHigh = nBinsHist;
          mult = binLow - 1;
          break;
      }
      
      // get the integral
      double integral = h->Integral(binLow, binHigh);
      
      // and the fraction of the total integral
      double fraction = integral / norm;
      
      // we need to round the fraction to avoid some rounding errors in the
      // floating point arithmetic, which of course involves more errors...
      fraction = Round(fraction, 12);
      
      // now loop to find if it sits on a bin edge
      switch (direction) {
        case Integral::Forward :
          for (int bin = centBin; bin < nCentBins; ++bin) {
            if (1.0 - fraction < centCuts[bin]) {
              boundaries[centBin++] = mult;
            }
          }
          break;
        case Integral::Backward :
          for (int bin = centBin; bin < nCentBins; ++bin) {
            if (fraction >= centCuts[bin]) {
              boundaries[centBin++] = mult;
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
  
  std::pair<vector<double>, unique_ptr<TH1D>>
  Centrality::weights(unsigned fit_boundary) {
    // make sure both histograms exist
    if (data_.get() == nullptr) {
      LOG(ERROR) << "no data refmult distribution loaded, can't calculate weights";
      return {vector<double>(), unique_ptr<TH1D>()};
    }
    
    if (simu_.get() == nullptr) {
      LOG(ERROR) << "no simulated refmult distribution loaded, can't calculate weights";
      return {vector<double>(), unique_ptr<TH1D>()};
    }
    
    // now make a copy of the simulated refmult and take the ratio
    unique_ptr<TH1D> ratio = make_unique<TH1D>(*simu_.get());
    ratio->Divide(data_.get());
    
    // define the fit function
    string fxn_def = "[0] + [1]/([2]*x + [3]) + [4]*([2]*x + [3]) + [5]/([2]*x + [3])^2 + [6]*([2]*x + [3])^2";
    unique_ptr<TF1> ratio_fit = make_unique<TF1>("ratio_fit", fxn_def.c_str(), 0, fit_boundary);
    
    // define some reasonable defaults to help the fit
    ratio_fit->SetParameters(1.5, -20, 1, 6, 1e-3, 600, 5e-6);
    
    // and fit
    ratio->Fit(ratio_fit.get(), "EMR");
    vector<double> weights(7);
    for (int par = 0; par < ratio_fit->GetNpar(); ++par)
      weights[par] = ratio_fit->GetParameter(par);
    
    return {weights, std::move(ratio)};
  }
  
} // namespace sct
