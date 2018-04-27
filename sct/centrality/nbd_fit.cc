// sct/centrality/nbd_fit.cc

#include "sct/centrality/nbd_fit.hh"

#include <iomanip>

#include "sct/core/logging.hh"
#include "sct/core/enumerations.hh"
#include "sct/utils/negative_binomial.hh"
#include "sct/utils/random.hh"

#include "TFile.h"

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
  
  
  NBDFit::NBDFit(TH1D* data, TH2D* glauber)
  : multiplicityModel_(nullptr), refMultData_(nullptr), nPartnColl_(nullptr),
  refMultSim_(nullptr), minMultFit_(100) {
  
    if (data != nullptr)
      loadData(*data);
    
    if (glauber != nullptr)
      loadGlauber(*glauber);
  }
  
  NBDFit::~NBDFit() {
    
  }
  
  void NBDFit::loadData(const TH1D& data) {
    // first clear the old histogram
    refMultData_.reset();
    
    // copy the data histogram
    refMultData_ = std::make_shared<TH1D>(data);
    refMultData_->SetName(MakeString("nbdfit_internal_data_",
                          Random::instance().counter()).c_str());
    refMultData_->SetDirectory(0);
  }
  
  void NBDFit::loadGlauber(const TH2D& glauber) {
    // first clear the old histogram
    nPartnColl_.reset();
    
    // copy the data histogram
    nPartnColl_ = std::make_shared<TH2D>(glauber);
    nPartnColl_->SetName(MakeString("nbdfit_internal_data_",
                         Random::instance().counter()).c_str());
    nPartnColl_->SetDirectory(0);
  }
  
  // When using Fit(...) must set the NBD parameters beforehand
  void NBDFit::setParameters(double npp, double k, double x, double ppEff, double AAEff,
                             double centMult, double triggerBias, bool constEfficiency) {
    multiplicityModel_ = make_unique<MultiplicityModel>(npp, k, x, ppEff, AAEff, centMult,
                                                        triggerBias, constEfficiency);
  }
  
  std::shared_ptr<FitResult> NBDFit::fit(unsigned nevents, string name) {
    // fit real data w/ simulated multiplicity distribution
    
    // create output FitResults
    std::shared_ptr<FitResult> result = std::make_shared<FitResult>();
    
    // first make sure refmult & npartncoll have been loaded
    if (refMultData_ == nullptr) {
      LOG(ERROR) << "no data refmult distribution has been loaded: Fit failure";
      return result;
    }
    
    if (nPartnColl_ == nullptr) {
      LOG(ERROR) << "no glauber nPartnColl distribution has been loaded: Fit failure";
      return result;
    }
    
    // make the simulated refmult histogram with the same bin edges as our
    // data refmult distribution
    refMultSim_.reset();
    string histName;
    if (name == "" || name == "refmultsim")
      histName = MakeString("refmultsim", Random::instance().counter());
    else
      histName = name;
    refMultSim_ = std::make_shared<TH1D>(name.c_str(), "",
                                         refMultData_->GetXaxis()->GetNbins(),
                                         refMultData_->GetXaxis()->GetXmin(),
                                         refMultData_->GetXaxis()->GetXmax());
    refMultSim_->SetDirectory(0);
    refMultSim_->Sumw2();
    
    // now fill the simulated refmult distribution from the negative binomial,
    // with the MC glauber npart x ncoll distribution, nevents times
    for (int event = 0; event < nevents; ++event) {
      // first sample from the npart ncoll distribution
      double npart, ncoll;
      nPartnColl_->GetRandom2(npart, ncoll);
      
      // check if any collisions took place
      if (npart < 2 || ncoll < 1) continue;
      
      unsigned mult = multiplicityModel_->multiplicity(static_cast<int>(npart), static_cast<int>(ncoll));
      refMultSim_->Fill(mult);
    }
    
    // normalize
    double norm_ = norm(refMultData_.get(), refMultSim_.get());
    refMultSim_->Scale(norm_);
    
    // get chi2
    std::pair<double, int> chi2_res = chi2(refMultData_.get(), refMultSim_.get());
    
    // fill in the fit results
    result->chi2 = chi2_res.first;
    result->ndf  = chi2_res.second;
    result->data = refMultData_;
    result->simu = refMultSim_;
    
    return result;
  }
  
  unordered_map<string, shared_ptr<FitResult>>
  NBDFit::scan(unsigned nevents, unsigned npp_bins, double npp_min,
               double npp_max, unsigned k_bins, double k_min,
               double k_max, unsigned x_bins, double x_min, double x_max,
               double ppEff, double AAEff, double centMult, double triggerBias,
               bool constEfficiency, bool saveAllHist) {
    // Perform the fit routine over a grid of NBD values (Npp, K, X)
    // and return a dictionary of results
    unordered_map<string, shared_ptr<FitResult>> result_map;
    
    // total number of bins
    unsigned nBins = npp_bins * k_bins * x_bins;
    
    // define the step widths in all three parameters
    double dNpp = (npp_max - npp_min) / static_cast<double>(npp_bins);
    double dK   = (k_max - k_min) / static_cast<double>(k_bins);
    double dX   = (x_max - x_min) / static_cast<double>(x_bins);
    
    // for book-keeping
    double best_chi2 = 0.0;
    string best_chi2_key;
    
    for (int bin_npp = 0; bin_npp < npp_bins; ++bin_npp) {
      for (int bin_k = 0; bin_k < k_bins; ++bin_k) {
        for (int bin_x = 0; bin_x < x_bins; ++bin_x) {
          
          // get the current values
          double npp = npp_min + dNpp * bin_npp;
          double k   = k_min + dK * bin_k;
          double x   = x_min + dX * bin_x;
          
          setParameters(npp, k, x, ppEff, AAEff, centMult, triggerBias, constEfficiency);
          
          string key = MakeString("npp_", npp, "_k_", k, "_x_", x);
          
          shared_ptr<FitResult> result = fit(nevents, key);
          
          // if we only save the best fit, we will check if this fit is better
          // than the current result and update
          if (saveAllHist == false) {
            if (best_chi2_key.empty() || best_chi2 <= 0.0) {
              best_chi2 = result->chi2 / result->ndf;
              best_chi2_key = key;
            }
            else {
              double current_chi2 = result->chi2 / result->ndf;
              if (current_chi2 < best_chi2) {
                result_map[best_chi2_key]->simu.reset();
                best_chi2 = current_chi2;
                best_chi2_key = key;
              }
            }
          }
          
          // calculate current entry
          unsigned current_bin = bin_x + bin_k * x_bins + bin_npp * x_bins * k_bins;
          if (current_bin % 10 == 0) {
            LOG(INFO) << "Scan " << std::setprecision(2) << std::fixed
                      << (double)current_bin/nBins * 100.0 <<"% complete: current entry: ";
            LOG(INFO) << "Npp: " << npp << " k: " << k << " x: " << x;
            LOG(INFO) << "chi2/ndf: " << result->chi2 / result->ndf;
          }
          
          // add the result from this (npp, k, x) set to the dictionary
          result_map[key] = result;
        }
      }
    }
    return result_map;
  }
  
  double NBDFit::norm(TH1D* h1, TH1D* h2) {
    // get the normalization between two histograms in the region from
    // minMultFit_ to h1->GetXaxis()->GetXmax()
    double min = minMultFit_;
    int minBin = h1->GetXaxis()->FindBin(min);
    double max = h1->GetXaxis()->GetXmax();
    int maxBin = h1->GetXaxis()->FindBin(max);
    
    // calculate numerator and denominator
    double numerator = 0;
    double denominator = 0;
    for (int i = minBin; i <= maxBin; ++i) {
      double n1 = h1->GetBinContent(i);
      double n1Error = h1->GetBinContent(i);
      double n2 = h2->GetBinContent(i);
      
      if (n1 == 0.0 || n1Error == 0.0) continue;
      
      numerator += n1 * n2 / pow(n1Error, 2.0);
      denominator += n2 * n2 / pow(n1Error, 2.0);
    }
    return (denominator == 0.0 ? 1.0 : numerator / denominator);
  }
  
  std::pair<double, int> NBDFit::chi2(TH1D* h1, TH1D* h2) {
    
    //calculate the chi2 from data and simulation - return a pair of (chi2, ndf)
    // note: only calculates from minMultFit_ to max bin, same range used for the
    // normalization. We'll use ROOT's chi2 test by using set range
    std::pair<double, int> result;
    
    double min = minMultFit_;
    int minBin = h1->GetXaxis()->FindBin(min+0.001);
    int maxBin = h2->GetXaxis()->GetNbins();
    
    // set bin range
    h1->GetXaxis()->SetRange(minBin, maxBin);
    h2->GetXaxis()->SetRange(minBin, maxBin);
    
    double chi2;
    int ndf, good;
    string opts = "UU NORM";
    
    h1->Chi2TestX(h2, chi2, ndf, good, opts.c_str());
    
    // reset bin ranges
    h1->GetXaxis()->SetRange();
    h2->GetXaxis()->SetRange();
    
    result.first = chi2;
    result.second = ndf;
    return result;
  }
  
} // namespace sct
