// sct/systematics/histogram_collection.cc

#include "sct/systematics/histogram_collection.hh"

namespace sct {
  
  // instantiate for common histograms
  template class HistogramCollection<string, TH1D>;
  template class HistogramCollection<string, TH2D>;
  template class HistogramCollection<string, TH3D>;
  template class HistogramCollection<string, TProfile>;
  template class HistogramCollection<string, TProfile2D>;
  
} // namespace sct
