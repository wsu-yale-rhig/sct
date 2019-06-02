#ifndef SCT_CENTRALITY_NBD_FIT_H
#define SCT_CENTRALITY_NBD_FIT_H

/* Given a data reference multiplicity, a sct NPart vs NColl
 * distribution, and a set of NBD parameters, will generate a simulated
 * refmult distribution, and fit the simulated distribution to the data:
 * NBDFit fitter(data_file_name, glauber_file_name, data_hist_name);
 * fitter.setParameters(...);
 * fitter.fit(n_events_to_simulate);
 *
 * Can also scan a grid of Npp, K, X values to find the minimum chi2/NBD of
 * the set.
 * NBDFit fitter(data_file_name, glauber_file_name);
 * fitter.scan(...);
 */

#include "sct/centrality/multiplicity_model.h"
#include "sct/lib/enumerations.h"
#include "sct/lib/map.h"
#include "sct/lib/memory.h"

#include "TH1.h"
#include "TH2.h"

namespace sct {
class NegativeBinomial;

struct FitResult {
  double chi2;
  int ndf;
  TH1D* data;
  unique_ptr<TH1D> simu;

  FitResult() : chi2(0.0), ndf(0), data(nullptr), simu(nullptr){};
};

class NBDFit {
 public:
  // If no data file or glauber file is specified during construction, must
  // use loadData() and loadGlauber() before calling fit() or scan().
  NBDFit(TH1D* data = nullptr, TH2D* glauber = nullptr);

  virtual ~NBDFit();

  // To load in individually the data or glauber histograms, clears any
  // older histogram that was already loaded.
  void loadData(const TH1D& data);
  void loadGlauber(const TH2D& glauber);

  // Can perform centrality definition calculation
  void makeCentDefs(bool flag = true);

  // When using fit(...) must set the NBD parameters beforehand
  void setParameters(double npp, double k, double x, double pp_eff, double aa_eff,
                     double cent_mult, double trigger_bias, bool const_efficiency);

  // From the given NPart x NColl distribution, samples nevents times,
  // and generates a refmult distribution (with name name) from a negative
  // binomial. Then normalizes the simulated distribution to the data, and
  // reports the chi2 & NDF for that fit.
  unique_ptr<FitResult> fit(unsigned nevents = 10000,
                            string name = "refmultsim");

  // When using scan(...), not necessary to call setParameters(...).
  // scan() will perform the simulation and fit for a 3D grid of Npp,
  // K, X values passed by the user and will return all results
  // If saveAllHist == true, then every simulated refmult distribution
  // is saved - otherwise, only the best fit is saved
  sct_map<string, unique_ptr<FitResult>> scan(
      unsigned nevents, unsigned npp_bins, double npp_min, double npp_max,
      unsigned k_bins, double k_min, double k_max, unsigned x_bins,
      double x_min, double x_max, double pp_eff, double aa_eff, double cent_mult,
      double trigger_bias, bool const_efficiency, bool save_all_hist = false);

  // to restrict the fits to multiplicity > minmult_fit_, which will have an
  // effect on the chi2, since the low multiplicity regime is where the data
  // will deviate from the glauber simulation.
  void minimumMultiplicityCut(unsigned min) { minmult_fit_ = min; }
  inline unsigned minimumMultiplicityCut() const { return minmult_fit_; }

  // get normalization between two histograms in range
  // (minMultFit < x < h1->GetXaxis()->GetXmax());
  double norm(TH1D* h1, TH1D* h2);

  // get chi2 difference between h1 & h2
  std::pair<double, int> chi2(TH1* h1, TH1* h2);

  // use homebrewed chi2 copied from StGlauber library, instead of ROOTs
  // NOTE: StGlauber Chi2 is not symmetric and only uses errors from h1
  // in the denominator
  void useStGlauberChi2(bool flag = true) { use_stglauber_chi2_ = flag; }
  inline bool usingStGlauberChi2() { return use_stglauber_chi2_; }

 private:
  std::pair<double, int> chi2_root(TH1* h1, TH1* h2);
  std::pair<double, int> chi2_stglauber(TH1* h1, TH1* h2);

  // multiplicity model
  unique_ptr<MultiplicityModel> multiplicity_model_;

  // we will make local copies, outside of the ROOT files, so
  // we have to explicity manage the memory...
  unique_ptr<TH1D> refmult_data_;
  unique_ptr<TH2D> npart_ncoll_;

  // minimum multiplicity for fitting range
  double minmult_fit_;

  // flag for using StGlauber chi2 or ROOT's own chi-square minimization
  bool use_stglauber_chi2_;
};
}  // namespace sct

#endif  // SCT_CENTRALITY_NBD_FIT_H
