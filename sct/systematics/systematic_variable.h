#ifndef SCT_SYSTEMATICS_SYSTEMATIC_VARIABLE_H
#define SCT_SYSTEMATICS_SYSTEMATIC_VARIABLE_H

// a histogram container for GlauberSystematics,
// that allows similar sets of histograms to be made
// for each systematic variation, and provides a 
// single fill interface for each event

#include "sct/lib/enumerations.h"
#include "sct/systematics/histogram_collection.h"
#include "sct/lib/memory.h"

#include <set>
#include <vector>

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"

namespace sct {

  template<class T>
  using EventDict = sct_map<GlauberObservable, T, EnumClassHash>;
  
  class SystematicVariable {
  public:

    SystematicVariable();
    
    // calls Init(obs) on initialization
    SystematicVariable(GlauberObservable obs);
    
    ~SystematicVariable() {};

    void init(GlauberObservable obs);
    
    // builds extra histograms to calculate the even cumulants (2, 4, 6)
    // should be called before Add()
    // CAUTION: this has been copied from the old STAR maker without validation of the
    // correctness of the math or the algorithm. It has only been confirmed to produce 
    // the same results as the old maker.
    void calculateCumulants(bool flag = true) {}

    // initializes histograms for a specific glauber variation
    void add(GlauberMod mod);

    // delete all histograms
    void clear();
    
    bool fillEvent(GlauberMod mod, EventDict<double>& dict, double weight);
    
    // write histograms to specified TFile (will not change current TDirectory to file)
    void write(TFile* file);

  private:

    // returns a name for a histogram given the specifics
    string getHistogramName(GlauberMod mod, GlauberObservable x_axis, string tag, unsigned cumulant_order = 0);

    // used to build the properly weighted histograms before writing 
    void reweight(TProfile* source, TProfile* weight, TH1D* target);

    // calculates the nth order cumulant from profile values
    // these functions have essentially been copied from the old maker for completeness -
    // the math has not been verified. 
    double nthOrderCumulant(std::vector<double> moments, unsigned order);
    double nthOrderCumulantError(std::vector<double> moments, std::vector<double> errors, unsigned order);
    
    GlauberObservable y_; // variable flag
    string name_; // variable name (ncoll, npart, etc)
    string label_; // axis label for plotting (N_{part}, N_{coll}, etc)
    unsigned bins_;
    double low_;
    double high_;

    bool cumulant_flag_;
    std::set<unsigned> cumulant_order_ = {2, 4, 6}; // order of cumulants to calculate - odd cumulants will be zero

    // observables to calculate the SystematicVariable with respect to - e.g. for x_ in x_axis_variables_ : plot(x_, y_)
    std::set<GlauberObservable> x_axis_variables_ = {GlauberObservable::B, GlauberObservable::Ncoll, GlauberObservable::Npart, GlauberObservable::Multiplicity}; 

    std::set<GlauberMod> modifications_;
    
    // histogram collections
    HistogramCollection<TH1D> th1_;
    HistogramCollection<TH2D> th2_;
    HistogramCollection<TProfile> tprof_;
    HistogramCollection<TProfile> weights_;

    // cumulant histogram collections
    HistogramCollection<TH1D> moment_th1_ ;
    HistogramCollection<TProfile> moment_tprof_;
    HistogramCollection<TH1D> cumulant_th1_;
  };
  
} // namespace sct

#endif // SCT_SYSTEMATICS_SYSTEMATIC_VARIABLE_H
