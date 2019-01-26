#ifndef SCT_CENTRALITY_CENTRALITY_H
#define SCT_CENTRALITY_CENTRALITY_H

/* Given a refmult distribution, finds the
 * 5% bin boundaries defining the normal centrality
 * boundaries (0-5%, 5-10%, etc).
 */

#include "sct/core/base.h"
#include "sct/core/enumerations.h"

#include "TH1D.h"

namespace sct {
  class Centrality {
  public:
    
    // simulation and data must have some normalization between
    // them, or weights will not have any meaning 
    Centrality(TH1D* data = nullptr, TH1D* simu = nullptr);
    ~Centrality();
    
    void setDataRefmult(TH1D* histogram);
    void setSimuRefmult(TH1D* histogram);
    inline TH1D* dataRefmult() const {return data_.get();}
    inline TH1D* simuRefmult() const {return simu_.get();}
    
    // get the centrality bin edges, calculated from the simulated
    // refmult distribution by integrating bins of 5%
    vector<unsigned> centralityBins(XSecMod mod = XSecMod::None);
    
    // get the relative weight between the simulation and data refmult distribution
    // defined using the functional form:
    // [0] + [1]/([2]*x + [3]) + [4]*([2]*x + [3]) + [5]/([2]*x + [3])^2 + [6]*([2]*x + [3])^2"
    // returns the parameter values for the best fit for parameters 0-6
    std::pair<vector<double>, unique_ptr<TH1D>> weights(unsigned fit_cutoff = 400);
    
  private:
    
    vector<unsigned> integrate(TH1D* h, Integral direction,
                               XSecMod mod = XSecMod::None);
    
    unique_ptr<TH1D> data_;
    unique_ptr<TH1D> simu_;
    
  };
} // namespace sct

#endif // SCT_CENTRALITY_CENTRALITY_H
