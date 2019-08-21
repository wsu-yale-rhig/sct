#include "sct/utils/nucleus_info.h"

namespace sct {

NucleusInfo &NucleusInfo::instance() {
  static NucleusInfo instance_;
  return instance_;
}

NucleusInfo::NucleusInfo() {
  // initialize Au197
  mass_number_.insert({GlauberSpecies::Au197, 197});
  charge_number_.insert({GlauberSpecies::Au197, 79});
  radius_.insert({GlauberSpecies::Au197, 6.38});
  skin_depth_.insert({GlauberSpecies::Au197, 0.535});
  radius_error_.insert({GlauberSpecies::Au197, 0.12});
  skin_depth_error_.insert({GlauberSpecies::Au197, 0.054});
  beta2_.insert({GlauberSpecies::Au197, -0.131});
  beta4_.insert({GlauberSpecies::Au197, -0.031});
  hulthen_a_.insert({GlauberSpecies::Au197, 0.0});
  hulthen_b_.insert({GlauberSpecies::Au197, 0.0});
  name_.insert({GlauberSpecies::Au197, "Au197"});

  // initialize Sm154
  // no reported error on Sm154 radius, assign +- 0.005
  mass_number_.insert({GlauberSpecies::Sm154, 154});
  charge_number_.insert({GlauberSpecies::Sm154, 62});
  radius_.insert({GlauberSpecies::Sm154, 5.9387});
  skin_depth_.insert({GlauberSpecies::Sm154, 0.522});
  radius_error_.insert({GlauberSpecies::Sm154, 0.01});
  skin_depth_error_.insert({GlauberSpecies::Sm154, 0.030});
  beta2_.insert({GlauberSpecies::Sm154, 0.270});
  beta4_.insert({GlauberSpecies::Sm154, 0.113});
  hulthen_a_.insert({GlauberSpecies::Sm154, 0.0});
  hulthen_b_.insert({GlauberSpecies::Sm154, 0.0});
  name_.insert({GlauberSpecies::Sm154, "Sm154"});

  // initialize U238
  // no error on radius is reported for uranium, instead,
  // use the difference of the two measurements in the referenced
  // paper as our 1 sigma error, and double it for 2 sigma
  // R = 6.8054, d = 0.605 +/- 0.016
  // R = 6.874,  d = 0.556
  // dR = 0.0686 da = 0.049
  mass_number_.insert({GlauberSpecies::U238, 238});
  charge_number_.insert({GlauberSpecies::U238, 92});
  radius_.insert({GlauberSpecies::U238, 6.8054});
  skin_depth_.insert({GlauberSpecies::U238, 0.605});
  radius_error_.insert({GlauberSpecies::U238, 0.137});
  skin_depth_error_.insert({GlauberSpecies::U238, 0.098});
  beta2_.insert({GlauberSpecies::U238, 0.215});
  beta4_.insert({GlauberSpecies::U238, 0.093});
  hulthen_a_.insert({GlauberSpecies::U238, 0.0});
  hulthen_b_.insert({GlauberSpecies::U238, 0.0});
  name_.insert({GlauberSpecies::U238, "U238"});

  // initialize Pb208
  // here we use the Pb207 radius & skin depth...
  mass_number_.insert({GlauberSpecies::Pb208, 208});
  charge_number_.insert({GlauberSpecies::Pb208, 82});
  radius_.insert({GlauberSpecies::Pb208, 6.62});
  skin_depth_.insert({GlauberSpecies::Pb208, 0.546});
  radius_error_.insert({GlauberSpecies::Pb208, 0.120});
  skin_depth_error_.insert({GlauberSpecies::Pb208, 0.020});
  beta2_.insert({GlauberSpecies::Pb208, 0.000});
  beta4_.insert({GlauberSpecies::Pb208, 0.000});
  hulthen_a_.insert({GlauberSpecies::Pb208, 0.0});
  hulthen_b_.insert({GlauberSpecies::Pb208, 0.0});
  name_.insert({GlauberSpecies::Pb208, "Pb208"});

  // initialize Cu63
  mass_number_.insert({GlauberSpecies::Cu63, 63});
  charge_number_.insert({GlauberSpecies::Cu63, 29});
  radius_.insert({GlauberSpecies::Cu63, 4.218});
  skin_depth_.insert({GlauberSpecies::Cu63, 0.596});
  radius_error_.insert({GlauberSpecies::Cu63, 0.028});
  skin_depth_error_.insert({GlauberSpecies::Cu63, 0.010});
  beta2_.insert({GlauberSpecies::Cu63, 0.162});
  beta4_.insert({GlauberSpecies::Cu63, -0.006});
  hulthen_a_.insert({GlauberSpecies::Cu63, 0.0});
  hulthen_b_.insert({GlauberSpecies::Cu63, 0.0});
  name_.insert({GlauberSpecies::Cu63, "Cu63"});

  // initialize proton
  // radius is set to a small non-zero number to simulate a delta function
  // distribution
  mass_number_.insert({GlauberSpecies::p1, 1});
  charge_number_.insert({GlauberSpecies::p1, 1});
  radius_.insert({GlauberSpecies::p1, 1e-5});
  skin_depth_.insert({GlauberSpecies::p1, 0.0});
  radius_error_.insert({GlauberSpecies::p1, 0.0});
  skin_depth_error_.insert({GlauberSpecies::p1, 0.0});
  beta2_.insert({GlauberSpecies::p1, 0.0});
  beta4_.insert({GlauberSpecies::p1, 0.0});
  hulthen_a_.insert({GlauberSpecies::p1, 0.0});
  hulthen_b_.insert({GlauberSpecies::p1, 0.0});
  name_.insert({GlauberSpecies::p1, "p1"});

  // initialize proton
  // deuteron does not use the normal parameters - instead, it uses the Hulthen
  // form, which can be found in sct/utils/functions.h
  mass_number_.insert({GlauberSpecies::d2, 2});
  charge_number_.insert({GlauberSpecies::d2, 1});
  radius_.insert({GlauberSpecies::d2, 0.0});
  skin_depth_.insert({GlauberSpecies::d2, 0.0});
  radius_error_.insert({GlauberSpecies::d2, 0.0});
  skin_depth_error_.insert({GlauberSpecies::d2, 0.0});
  beta2_.insert({GlauberSpecies::d2, 0.0});
  beta4_.insert({GlauberSpecies::d2, 0.0});
  hulthen_a_.insert({GlauberSpecies::d2, 0.228});
  hulthen_b_.insert({GlauberSpecies::d2, 1.177});
  name_.insert({GlauberSpecies::d2, "d2"});
}

NucleusInfo::~NucleusInfo() {}

} // namespace sct
