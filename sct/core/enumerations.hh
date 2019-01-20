// sct/core/enumerations.hh

#ifndef SCT_CORE_ENUMERATIONS_HH
#define SCT_CORE_ENUMERATIONS_HH

#include "sct/core/base.hh"

#include <set>

namespace  sct {
  
  // collision systems that we have centrality tables & NBD parameters for
  enum class CollisionSystem {
    AuAu200,                   // Au+Au at 200 GeV
    AuAu62,                    // Au+Au at 62 GeV
    AuAu39,                    // Au+Au at 39 GeV
    AuAu27,                    // Au+Au at 27 GeV
    AuAu19,                    // Au+Au at 19 GeV
    AuAu11,                    // Au+Au at 11 GeV
    dAu200                     // d+Au  at 200 GeV
  };
  
  // nuclei for which we have default parameters for radius, skin depth
  // and eccentricity
  enum class GlauberSpecies {
    Au197,       // gold 197
    Sm154,       // samarium 154
    U238,        // uranium 238
    Pb208,       // lead 208
    Cu63         // copper 63
  };
  
  static std::unordered_map<GlauberSpecies, string, EnumClassHash> speciesString {
    {GlauberSpecies::Au197, "au197"}, {GlauberSpecies::Sm154, "sm154"},
    {GlauberSpecies::U238, "u238"}, {GlauberSpecies::Pb208, "pb208"},
    {GlauberSpecies::Cu63, "cu63"} };
  
  // collision energies for which we have hard coded nucleon-nucleon inelastic xsections
  enum class CollisionEnergy {
    E2760 = 2760,   // 2.76 TeV
    E200 = 200,     // 200 GeV
    E62 = 62,       // 62.4 GeV
    E39 = 39,       // 39 GeV
    E27 = 27,       // 27 GeV
    E19 = 19,       // 19.6 GeV
    E14 = 14,       // 14.5 GeV
    E11 = 11,       // 11.5 GeV
    E7 = 7          // 7.7 GeV
  };
  
  // collision profile can either be a hard core, or gaussian
  enum class CollisionProfile {
    HardCore,
    Gaussian
  };
  
  // positional smearing of nucleons during glauber modeling
  // default is off (pure woods-saxon), can either smear in a hard sphere with a
  // constant probability, or with a 3D gaussian distribution
  enum class NucleonSmearing {
    None,
    HardCore,
    Gaussian
  };
  
  // For systematics: vary the settings of the glauber model
  enum class GlauberMod {
    Nominal,      // nominal settings
    Large,        // Large radius(+2%), small skin depth(-10%)
    Small,        // Small radius(-2%), large skin depth(+10%)
    LargeXSec,    // Large inelastic NN cross section (+1mb)
    SmallXSec,    // Small inelastic NN cross section (-1mb)
    Gauss         // Gaussian collision profile
  };

  // set of all glauber modifications as strings
  static std::set<string> glauberModString{"nominal", "large", "small", "largexsec",
    "smallxsec", "gauss"};
  
  // dictionaries for glaubermod <==> string
  static std::unordered_map<sct::GlauberMod, std::string, sct::EnumClassHash> modToString {
    {GlauberMod::Nominal, "nominal"}, {GlauberMod::Large, "large"}, {GlauberMod::Small, "small"},
    {GlauberMod::LargeXSec, "largexsec"}, {GlauberMod::SmallXSec, "smallxsec"},
    {GlauberMod::Gauss, "gauss"} };
  static std::unordered_map<std::string, sct::GlauberMod> stringToMod {
    {"nominal", GlauberMod::Nominal}, {"large", GlauberMod::Large}, {"small", GlauberMod::Small},
    {"largexsec", GlauberMod::LargeXSec}, {"smallxsec", GlauberMod::SmallXSec},
    {"gauss", GlauberMod::Gauss} };
  
  // For systematics: vary the NBD Npp by 5%
  // also has options for 39 GeV AuAu to use reweighting (see centrality_table.hh)
  enum class NppMod {
    Nominal,                   // The nominal centrality
    Low,                       // low npp, high x ( -5% npp )
    High,                      // high npp, low x ( +5% npp )
    LowReweighted,             // Nominal - reweighted (used for 39 GeV AuAu only)
    HighReweighted,            // Nominal + reweighted (used for 39 GeV AuAu only)
  };
  
  // For systematics: vary the NN cross section by 5%
  enum class XSecMod {
    None,                      // default
    Minus5,                    // -5% total cross section
    Plus5                      // +5% total cross section
  };

  // glauber model observables - used in systematic analysis
  enum class GlauberObservable {
    npart,
    ncoll,
    multiplicity,
    b,
    centrality,
    areaRP,
    areaPP,
    eccRP,
    eccPP2,
    eccPP3,
    eccPP4
  };
  
  // For sanity checks on code & histograms: when performing a centrality calculation,
  // it can be performed forwards & backwards to compare results
  enum class Integral {
    Forward,
    Backward
  };
  
  
  static const unsigned nGlauberWeights = 3; // 0: nPart, 1: ncoll, 3: spectators
  enum class GlauberWeights {
    NPart = 0,        // weight by number of participants
    NColl = 1,        // weight by number of collisions
    Spectators = 2    // weight by number of spectators
  };
  
} // namespace sct

#endif // SCT_CORE_ENUMERATIONS_HH
