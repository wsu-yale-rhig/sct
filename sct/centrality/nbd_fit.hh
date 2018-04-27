// sct/centrality/nbd_fit.hh

#ifndef SCT_CENTRALITY_NBD_FIT_HH
#define SCT_CENTRALITY_NBD_FIT_HH

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
 
#include "sct/core/base.hh"
#include "sct/core/enumerations.hh"
#include "sct/centrality/multiplicity_model.hh"

#include "TH1.h"
#include "TH2.h"

namespace sct {
  class NegativeBinomial;
  
  struct FitResult {
    double chi2;
    int ndf;
    std::shared_ptr<TH1D> data;
    std::shared_ptr<TH1D> simu;
    
    FitResult() : chi2(0.0), ndf(0),
      data(nullptr), simu(nullptr) {};
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
    void setParameters(double npp, double k, double x, double ppEff, double AAEff,
                       double centMult, double triggerBias, bool constEfficiency);
    
    // From the given NPart x NColl distribution, samples nevents times,
    // and generates a refmult distribution (with name name) from a negative binomial.
    // Then normalizes the simulated distribution to the data, and
    // reports the chi2 & NDF for that fit.
    std::shared_ptr<FitResult> fit(unsigned nevents = 1000, string name = "refmultsim");
    
    // When using scan(...), not necessary to call setParameters(...).
    // scan() will perform the simulation and fit for a 3D grid of Npp,
    // K, X values passed by the user and will return all results
    // If saveAllHist == true, then every simulated refmult distribution
    // is saved - otherwise, only the best fit is saved
    unordered_map<string, shared_ptr<FitResult>>
    scan(unsigned nevents, unsigned npp_bins, double npp_min, double npp_max,
         unsigned k_bins, double k_min, double k_max, unsigned x_bins, double x_min,
         double x_max, double ppEff, double AAEff, double centMult, double triggerBias,
         bool constEfficiency, bool saveAllHist = false);
    
    // to restrict the fits to multiplicity > minMultFit_, which will have an effect
    // on the chi2, since the low multiplicity regime is where the data will deviate
    // from the glauber simulation.
    void minimumMultiplicityCut(unsigned min) {minMultFit_ = min;}
    inline unsigned minimumMultiplicityCut() const {return minMultFit_;}

    // get normalization between two histograms in range
    // (minMultFit < x < h1->GetXaxis()->GetXmax());
    double norm(TH1D* h1, TH1D* h2);
    
    // get chi2 difference between h1 & h2
    std::pair<double, int> chi2(TH1D* h1, TH1D* h2);
    
  private:
    
    // multiplicity model
    unique_ptr<MultiplicityModel> multiplicityModel_;
    
    // we will make local copies, outside of the ROOT files, so
    // we have to explicity manage the memory...
    shared_ptr<TH1D> refMultData_;
    shared_ptr<TH2D> nPartnColl_;
    
    // simulated refmult using NBD
    shared_ptr<TH1D> refMultSim_;
    
    // minimum multiplicity for fitting range
    double minMultFit_;
  };
} // namespace sct
  
#endif // SCT_CENTRALITY_NBD_FIT_HH
