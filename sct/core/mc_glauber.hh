// sct/core/mc_glauber.hh

#ifndef SCT_CORE_MC_GLAUBER_HH
#define SCT_CORE_MC_GLAUBER_HH

/* Monte Carlo Glauber generator
 * For an overview of Glauber modeling, see
 * https://arxiv.org/abs/nucl-ex/0701025
 *
 * Impact parameter is drawn assuming the differential cross section
 * dsigma / dB = const * B
 * impact parameter range can be set using setMinB() & setMaxB()
 *
 */

#include "sct/core/base.hh"
#include "sct/core/enumerations.hh"
#include "sct/core/glauber_tree.hh"
#include "sct/core/nucleon.hh"
#include "sct/core/nucleus.hh"
#include "sct/core/collision.hh"

#include "TH1D.h"

namespace sct {
  class MCGlauber {
  public:
    
    // Constructor that provides access to default settings
    // for specific nuclei (Au, Pb, Sm, Cu, U)
    // default construction is for AuAu at 200 GeV, spherical
    // with no modification to the glauber parameters
    MCGlauber(GlauberSpecies speciesA = GlauberSpecies::Au197,
              GlauberSpecies speciesB = GlauberSpecies::Au197,
              CollisionEnergy energy = CollisionEnergy::E200,
              GlauberMod mod = GlauberMod::Nominal,
              bool deformationA = false, bool deformationB = false);
    
    // generic symmetric collision constructor
    MCGlauber(unsigned massNumber,    // mass number of nucleus (number of nucleons)
              double radius,          // radius of nucleus
              double skinDepth,       // skin depth of nucleus
              double beta2,           // second order deformation parameter
              double beta4,           // 4th order deformation parameter
              double inelasticXSec,   // inelastic NN cross section
              double energy);         // energy (sqrt(sNN))
    
    // generic asymmetric collision cross section
    MCGlauber(unsigned massNumberA,   // mass number nucleus A
              double radiusA,         // radius nucleus A
              double skinDepthA,      // skin depth nucleus A
              double beta2A,          // second order deformation parameter nucleus A
              double beta4A,          // fourth order deformation parameter nucleus A
              unsigned massNumberB,   // mass number nucleus B
              double radiusB,         // radius nucleus B
              double skinDepthB,      // skin depth nucleus B
              double beta2B,          // second order deformation parameter nucleus B
              double beta4B,          // second order deformation parameter nucleus B
              double inelasticXSec,   // inelastic NN cross section
              double energy);         // energy (sqrt(sNN))
    
    virtual ~MCGlauber();
    
    // set minimum & maximum impact parameter
    void setImpactParameterRange(double b_min, double b_max);
    inline double minB() const {return b_min_;}
    inline double maxB() const {return b_max_;}
    
    // add repulsion distance to nucleons (default is 0 fm)
    // forces generated nucleons to be at least repulsionDistance() away from
    // each other inside a generated nucleus
    void setRepulsionDistance(double fm);
    inline double repulsionDistance() const {return nucleusA_->repulsionDistance();}
    
    // nucleon smearing (default is NucleonSmearing::None)
    // smears the positions of the nuclei from the initial Woods-Saxon distribution
    // can smear position in a uniform sphere, or with a Gaussian distribution
    void setSmearing(NucleonSmearing smear);
    NucleonSmearing smearing() const {return nucleusA_->nucleonSmearing();}
    
    // collision profile (default is CollisionProfile::HardCore)
    // using a HardCore profile, if the nucleons overlap (radius defined with the
    // nucleon-nucleon cross section), then a collision is counted. Using a
    // Gaussian profile, there is a gaussian probability for collision as a
    // function of the radial distance between nucleons (sigma defined w/ the
    // nucleon-nucleon cross section
    void setCollisionProfile(CollisionProfile profile);
    CollisionProfile collisionProfile() const {return collision_.collisionProfile();}
    
    // if one already has parameters for the two-part multiplicity model, either
    // from fits or some other source, they can be used to generate multiplicity
    // estimates in the glauber trees. This also allows observables like eccentricity
    // to be estimated with a multiplicity weighting
    void setMultiplicityModel(double npp, double k, double x, double ppEff,
                              double aaEff, double aaCent, double trigEff = 1.0,
                              bool constEff = false);
    
    // run the glauber MC for N events
    void run(unsigned N = 1000);
    
    // lookup table for known NN cross sections (as a function of energy)
    double lookupXSec(CollisionEnergy energy);
    
    // generates a single event: returns true if nColl > 0
    bool generate();
    
    // write the summary into header tree (called automatically by run())
    void writeHeader();
    
    // gives access to tree
    GlauberTree* results() const {return tree_.get();}
    
    // gives access to generated parameters
    TH2D* woodsSaxonA() {return nucleusA_->generatedRCosTheta();}
    TH2D* woodsSaxonB() {return nucleusB_->generatedRCosTheta();}
    TH1D* generatedImpactParameter() {return generated_ip_.get();}
    TH1D* acceptedImpactParameter()  {return accepted_ip_.get();}
    
    // clears current tree entry & nuclei
    void clear();
    
  private:
    
    void init(unsigned massNumberA, double radiusA, double skinDepthA,
              double beta2A, double beta4A, unsigned massNumberB,
              double radiusB, double skinDepthB, double beta2B, double beta4B);
    
    void initOutput();
    void initQA();
    
    // output tree
    unique_ptr<GlauberTree> tree_;
    
    // containers for nucleus A & B
    unique_ptr<Nucleus> nucleusA_;
    unique_ptr<Nucleus> nucleusB_;
    
    Collision collision_;
    
    // counters
    unsigned eventsGenerated_;
    unsigned eventsAccepted_;
    
    // impact parameter
    double b_min_;
    double b_max_;
    
    // collision parameters
    double  energy_;            // center of mass energy
    GlauberMod modification_;   // Used for systematics: define a chance in the glauber parameters
                                // default is GlauberMod::Nominal

    // QA histograms for impact parameter
    unique_ptr<TH1D> generated_ip_;
    unique_ptr<TH1D> accepted_ip_;
    
  };
} // namespace sct

#endif // SCT_CORE_MC_GLAUBER_HH
