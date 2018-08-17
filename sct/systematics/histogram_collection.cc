// sct/systematics/histogram_collection.cc

#include "sct/systematics/histogram_collection.hh"

namespace sct {
  
  // instantiate for common histograms
  template class HistogramCollection<TH1D>;
  template class HistogramCollection<TH2D>;
  template class HistogramCollection<TH3D>;
  template class HistogramCollection<TProfile>;
  template class HistogramCollection<TProfile2D>;
  
} // namespace sct
