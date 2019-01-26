#ifndef SCT_UTILS_HISTOGRAM_INFO_HH
#define SCT_UTILS_HISTOGRAM_INFO_HH

#include "sct/core/base.hh"
#include "sct/core/enumerations.hh"

// defines constants used in the creation of histograms - 
// number of bins, lower and upper bounds for common observables
// such as multiplicity, impact parameter, etc

namespace sct {

  class HistogramInfo {
  public:
    
    static HistogramInfo& instance();
    virtual ~HistogramInfo();

    unsigned bins(GlauberObservable obs) {return bins_[obs];}
    double   lowEdge(GlauberObservable obs) {return low_edge_[obs];}
    double   highEdge(GlauberObservable obs) {return high_edge_[obs];}
    string   name(GlauberObservable obs) {return name_[obs];}
    string   label(GlauberObservable obs) {return label_[obs];}
    
  private:
    
    unordered_map<GlauberObservable, unsigned, EnumClassHash> bins_;
    unordered_map<GlauberObservable, double, EnumClassHash> low_edge_;
    unordered_map<GlauberObservable, double, EnumClassHash> high_edge_;
    unordered_map<GlauberObservable, string, EnumClassHash> name_;
    unordered_map<GlauberObservable, string, EnumClassHash> label_;

    vector<double> CentralityLowerBound = {0, 5, 10, 15, 20, 25, 30, 
                       35, 40, 45, 50, 55, 60, 65, 70, 75, 80};
    vector<double> CentralityUpperBound = {5, 10, 15, 20, 25, 30, 35, 
                       40, 45, 50, 55, 60, 65, 70, 75, 80, 100};
    
    HistogramInfo();
    HistogramInfo(const HistogramInfo&);
    
  };

} // namespace sct


#endif // SCT_UTILS_HISTOGRAM_INFO_HH