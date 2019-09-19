#ifndef SCT_GLAUBER_MC_GLAUBER_H
#define SCT_GLAUBER_MC_GLAUBER_H

/* Monte Carlo Glauber generator
 * For an overview of Glauber modeling, see
 * https://arxiv.org/abs/nucl-ex/0701025
 *
 * Impact parameter is drawn assuming the differential cross section
 * dsigma / dB = const * B
 * impact parameter range can be set using setMinB() & setMaxB()
 *
 */

#include "sct/glauber/collision.h"
#include "sct/glauber/glauber_tree.h"
#include "sct/glauber/nucleon.h"
#include "sct/glauber/nucleus.h"
#include "sct/lib/enumerations.h"

#include "TH1D.h"

namespace sct {

// the return status of a single generated event - Error for a failure to
// generate an event, Miss for an event with zero binary collisions, Hit
// otherwise
enum class EventStatus { Error, Miss, Hit };

class MCGlauber {
public:
  // Constructor that provides access to default settings
  // for specific nuclei (Au, Pb, Sm, Cu, U, p, d)
  // default construction is for AuAu at 200 GeV, spherical
  // with no modification to the glauber parameters
  MCGlauber(GlauberSpecies species_A = GlauberSpecies::Au197,
            GlauberSpecies species_B = GlauberSpecies::Au197,
            CollisionEnergy energy = CollisionEnergy::E200,
            GlauberMod mod = GlauberMod::Nominal, bool deformation_A = false,
            bool deformation_B = false);

  // generic symmetric collision constructor
  MCGlauber(
      unsigned mass_number,  // mass number of nucleus (number of nucleons)
      NucleonPDF::PDF pdf,   // the PDF to use for the nucleon distribution.
                             // Available options are defined in nucleon_pdf.h
      parameter_list params, // relevant parameters for the
      double inelastic_xsec, // inelastic NN cross section
      double energy);        // energy (sqrt(sNN))

  // generic asymmetric collision cross section
  MCGlauber(unsigned mass_number_A,  // mass number nucleus A
            NucleonPDF::PDF pdf_A,   // PDF for nucleus A
            parameter_list params_A, // relevant parameters for the PDF for A
            unsigned mass_number_B,  // mass number nucleus B
            NucleonPDF::PDF pdf_B,   // the PDF for nucleus B
            parameter_list params_B, // relevant parameters for the PDF for B
            double inelastic_xsec,   // inelastic NN cross section
            double energy);          // energy (sqrt(sNN))

  virtual ~MCGlauber();

  // set minimum & maximum impact parameter
  void setImpactParameterRange(double b_min, double b_max);
  inline double minB() const { return b_min_; }
  inline double maxB() const { return b_max_; }

  // add repulsion distance to nucleons (default is 0 fm)
  // forces generated nucleons to be at least repulsionDistance() away from
  // each other inside a generated nucleus
  void setRepulsionDistance(double fm);
  inline double repulsionDistance() const {
    return nucleusA_->repulsionDistance();
  }

  // nucleon smearing (default is NucleonSmearing::None)
  // smears the positions of the nuclei from the initial Woods-Saxon
  // distribution can smear position in a uniform sphere, or with a Gaussian
  // distribution
  void setSmearing(NucleonSmearing smear);
  NucleonSmearing smearing() const { return nucleusA_->nucleonSmearing(); }

  // collision profile (default is CollisionProfile::HardCore)
  // using a HardCore profile, if the nucleons overlap (radius defined with the
  // nucleon-nucleon cross section), then a collision is counted. Using a
  // Gaussian profile, there is a gaussian probability for collision as a
  // function of the radial distance between nucleons (sigma defined w/ the
  // nucleon-nucleon cross section
  void setCollisionProfile(CollisionProfile profile);
  CollisionProfile collisionProfile() const {
    return collision_.collisionProfile();
  }

  // if one already has parameters for the two-part multiplicity model, either
  // from fits or some other source, they can be used to generate multiplicity
  // estimates in the glauber trees. This also allows observables like
  // eccentricity to be estimated with a multiplicity weighting
  void setMultiplicityModel(double npp, double k, double x, double pp_eff,
                            double aa_eff, double aa_cent,
                            double trig_eff = 1.0, bool const_eff = false);

  // run the glauber MC for N events
  void run(unsigned N = 1000);

  // lookup table for known NN cross sections (as a function of energy)
  double lookupXSec(CollisionEnergy energy);

  // generates a single event: returns Error if nuclei could not be generated,
  // Miss if there is no binary collision, and Hit if at least one binary
  // collision occured
  EventStatus generate();

  // write the summary into header tree (called automatically by run())
  void writeHeader();

  // gives access to tree and individual nuclei
  GlauberTree *results() const { return tree_.get(); }
  Nucleus *nucleusA() const { return nucleusA.get(); }
  Nucleus *nucleusB() const { return nucleusB.get(); }

  // gives access to generated parameters
  TH2D *nuclearPDFA() { return nucleusA_->generatedRCosTheta(); }
  TH2D *nuclearPDFB() { return nucleusB_->generatedRCosTheta(); }
  TH1D *generatedImpactParameter() { return generated_ip_.get(); }
  TH1D *acceptedImpactParameter() { return accepted_ip_.get(); }

  // clears current tree entry & nuclei
  void clear();

private:
  void init(GlauberSpecies species_A, GlauberSpecies species_B, GlauberMod mod,
            bool deformed_a, bool deformed_b);
  void init(unsigned mass_number_A, NucleonPDF::PDF pdf_A,
            parameter_list params_A, unsigned mass_number_B,
            NucleonPDF::PDF pdf_B, parameter_list params_B);

  void initOutput();
  void initQA();

  // output tree
  unique_ptr<GlauberTree> tree_;

  // containers for nucleus A & B
  unique_ptr<Nucleus> nucleusA_;
  unique_ptr<Nucleus> nucleusB_;

  Collision collision_;

  // counters
  unsigned events_generated_;
  unsigned events_accepted_;

  // impact parameter
  double b_min_;
  double b_max_;

  // collision parameters
  double energy_;           // center of mass energy
  GlauberMod modification_; // Used for systematics: define a chance in the
                            // glauber parameters default is GlauberMod::Nominal

  // QA histograms for impact parameter
  unique_ptr<TH1D> generated_ip_;
  unique_ptr<TH1D> accepted_ip_;
};
} // namespace sct

#endif // SCT_GLAUBER_MC_GLAUBER_H
