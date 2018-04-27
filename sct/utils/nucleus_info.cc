// sct/utils/nucleus_info.cc

#include "sct/utils/nucleus_info.hh"

namespace sct {
  
  NucleusInfo& NucleusInfo::instance() {
    static NucleusInfo instance_;
    return instance_;
  }
  
  NucleusInfo::NucleusInfo() {
    
    // initialize Au197
    massNumber_.insert({GlauberSpecies::Au197, 197});
    chargeNumber_.insert({GlauberSpecies::Au197, 79});
    radius_.insert({GlauberSpecies::Au197, 6.38});
    skinDepth_.insert({GlauberSpecies::Au197, 0.535});
    radiusError_.insert({GlauberSpecies::Au197, 0.12});
    skinDepthError_.insert({GlauberSpecies::Au197, 0.054});
    beta2_.insert({GlauberSpecies::Au197, -0.131});
    beta4_.insert({GlauberSpecies::Au197, -0.031});
    name_.insert({GlauberSpecies::Au197, "Au197"});
    
    // initialize Sm154
    // no reported error on Sm154 radius, assign +- 0.005
    massNumber_.insert({GlauberSpecies::Sm154, 154});
    chargeNumber_.insert({GlauberSpecies::Sm154, 62});
    radius_.insert({GlauberSpecies::Sm154, 5.9387});
    skinDepth_.insert({GlauberSpecies::Sm154, 0.522});
    radiusError_.insert({GlauberSpecies::Sm154, 0.01});
    skinDepthError_.insert({GlauberSpecies::Sm154, 0.030});
    beta2_.insert({GlauberSpecies::Sm154, 0.270});
    beta4_.insert({GlauberSpecies::Sm154, 0.113});
    name_.insert({GlauberSpecies::Sm154, "Sm154"});
    
    // initialize U238
    // no error on radius is reported for uranium, instead,
    // use the difference of the two measurements in the referenced
    // paper as our 1 sigma error, and double it for 2 sigma
    // R = 6.8054, d = 0.605 +/- 0.016
    // R = 6.874,  d = 0.556
    // dR = 0.0686 da = 0.049
    massNumber_.insert({GlauberSpecies::U238, 238});
    chargeNumber_.insert({GlauberSpecies::U238, 92});
    radius_.insert({GlauberSpecies::U238, 6.8054});
    skinDepth_.insert({GlauberSpecies::U238, 0.605});
    radiusError_.insert({GlauberSpecies::U238, 0.137});
    skinDepthError_.insert({GlauberSpecies::U238, 0.098});
    beta2_.insert({GlauberSpecies::U238, 0.215});
    beta4_.insert({GlauberSpecies::U238, 0.093});
    name_.insert({GlauberSpecies::U238, "U238"});
    
    // initialize Pb208
    // here we use the Pb207 radius & skin depth...
    massNumber_.insert({GlauberSpecies::Pb208, 208});
    chargeNumber_.insert({GlauberSpecies::Pb208, 82});
    radius_.insert({GlauberSpecies::Pb208, 6.62});
    skinDepth_.insert({GlauberSpecies::Pb208, 0.546});
    radiusError_.insert({GlauberSpecies::Pb208, 0.120});
    skinDepthError_.insert({GlauberSpecies::Pb208, 0.020});
    beta2_.insert({GlauberSpecies::Pb208, 0.000});
    beta4_.insert({GlauberSpecies::Pb208, 0.000});
    name_.insert({GlauberSpecies::Pb208, "Pb208"});
    
    // initialize Cu63
    massNumber_.insert({GlauberSpecies::Cu63, 63});
    chargeNumber_.insert({GlauberSpecies::Cu63, 29});
    radius_.insert({GlauberSpecies::Cu63, 4.218});
    skinDepth_.insert({GlauberSpecies::Cu63, 0.596});
    radiusError_.insert({GlauberSpecies::Cu63, 0.028});
    skinDepthError_.insert({GlauberSpecies::Cu63, 0.010});
    beta2_.insert({GlauberSpecies::Cu63, 0.162});
    beta4_.insert({GlauberSpecies::Cu63, -0.006});
    name_.insert({GlauberSpecies::Cu63, "Cu63"});
    
  }
  
  NucleusInfo::~NucleusInfo() {}
  
} // namespace sct
