#include "sct/utils/histogram_info.h"

namespace sct {

HistogramInfo& HistogramInfo::instance() {
  static HistogramInfo instance_;
  return instance_;
}

HistogramInfo::HistogramInfo() {
  // initialize each variable

  // impact parameter
  bins_.insert({GlauberObservable::B, 400});
  low_edge_.insert({GlauberObservable::B, 0.0});
  high_edge_.insert({GlauberObservable::B, 20.0});
  name_.insert({GlauberObservable::B, "b"});
  label_.insert({GlauberObservable::B, "b [fm]"});

  // number of collisions
  bins_.insert({GlauberObservable::Ncoll, 1800});
  low_edge_.insert({GlauberObservable::Ncoll, 0.0});
  high_edge_.insert({GlauberObservable::Ncoll, 1800.0});
  name_.insert({GlauberObservable::Ncoll, "ncoll"});
  label_.insert({GlauberObservable::Ncoll, "N_{coll}"});

  // number of participants
  bins_.insert({GlauberObservable::Npart, 500});
  low_edge_.insert({GlauberObservable::Npart, 0.0});
  high_edge_.insert({GlauberObservable::Npart, 500.0});
  name_.insert({GlauberObservable::Npart, "npart"});
  label_.insert({GlauberObservable::Npart, "N_{part}"});

  // centrality
  CentralityLowerBound = {0,  5,  10, 15, 20, 25, 30, 35, 40,
                          45, 50, 55, 60, 65, 70, 75, 80};
  CentralityUpperBound = {5,  10, 15, 20, 25, 30, 35, 40, 45,
                          50, 55, 60, 65, 70, 75, 80, 100};
  bins_.insert({GlauberObservable::Centrality, CentralityLowerBound.size()});
  low_edge_.insert({GlauberObservable::Centrality, 0});
  high_edge_.insert(
      {GlauberObservable::Centrality, CentralityLowerBound.size()});
  name_.insert({GlauberObservable::Centrality, "cent"});
  label_.insert({GlauberObservable::Centrality, "centrality"});

  // number of participants
  bins_.insert({GlauberObservable::Multiplicity, 2000});
  low_edge_.insert({GlauberObservable::Multiplicity, 0.0});
  high_edge_.insert({GlauberObservable::Multiplicity, 2000.0});
  name_.insert({GlauberObservable::Multiplicity, "mult"});
  label_.insert({GlauberObservable::Multiplicity, "multiplicity"});

  // area reaction plane
  bins_.insert({GlauberObservable::RPArea, 100});
  low_edge_.insert({GlauberObservable::RPArea, 0.0});
  high_edge_.insert({GlauberObservable::RPArea, 50.0});
  name_.insert({GlauberObservable::RPArea, "rparea"});
  label_.insert({GlauberObservable::RPArea, "area_{RP} [fm^2]"});

  // area participant plane
  bins_.insert({GlauberObservable::PPArea, 100});
  low_edge_.insert({GlauberObservable::PPArea, 0.0});
  high_edge_.insert({GlauberObservable::PPArea, 50.0});
  name_.insert({GlauberObservable::PPArea, "pparea"});
  label_.insert({GlauberObservable::PPArea, "area_{PP} [fm^2]"});

  // eccentricity reaction plane
  bins_.insert({GlauberObservable::RP2Ecc, 100});
  low_edge_.insert({GlauberObservable::RP2Ecc, -1.0});
  high_edge_.insert({GlauberObservable::RP2Ecc, 1.0});
  name_.insert({GlauberObservable::RP2Ecc, "rpecc"});
  label_.insert({GlauberObservable::RP2Ecc, "#epsilon_{RP}"});

  // eccentricity participant plane 2nd order
  bins_.insert({GlauberObservable::PP2Ecc, 100});
  low_edge_.insert({GlauberObservable::PP2Ecc, -1.0});
  high_edge_.insert({GlauberObservable::PP2Ecc, 1.0});
  name_.insert({GlauberObservable::PP2Ecc, "pp2ecc"});
  label_.insert({GlauberObservable::PP2Ecc, "#epsilon_{2,PP}"});

  // eccentricity participant plane 3rd order
  bins_.insert({GlauberObservable::PP3Ecc, 100});
  low_edge_.insert({GlauberObservable::PP3Ecc, -1.0});
  high_edge_.insert({GlauberObservable::PP3Ecc, 1.0});
  name_.insert({GlauberObservable::PP3Ecc, "pp3ecc"});
  label_.insert({GlauberObservable::PP3Ecc, "#epsilon_{3,PP}"});

  // eccentricity participant plane 4th order
  bins_.insert({GlauberObservable::PP4Ecc, 100});
  low_edge_.insert({GlauberObservable::PP4Ecc, -1.0});
  high_edge_.insert({GlauberObservable::PP4Ecc, 1.0});
  name_.insert({GlauberObservable::PP4Ecc, "pp4ecc"});
  label_.insert({GlauberObservable::PP4Ecc, "#epsilon_{4,PP}"});
}

HistogramInfo::~HistogramInfo() {}

}  // namespace sct
