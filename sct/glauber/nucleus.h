#ifndef SCT_GLAUBER_NUCLEUS_H
#define SCT_GLAUBER_NUCLEUS_H

/*  Nucleus: used internally by MCGlauber to store nucleus parameters,
 *  and to generate & store nucleons for each generated event.
 *
 *  Default constructor gives an invalid state to avoid ambiguity - must
 *  set the parameters by hand.
 *
 *  To generate a new nucleus from the Woods-Saxon density profile, call
 *  Generate(b), where b is the impact parameter (constant offset along the
 *  x axis for all generated nucleons)
 */
#include "sct/glauber/nucleon.h"
#include "sct/lib/enumerations.h"
#include "sct/lib/memory.h"

#include <string>
#include <vector>

#include "TF1.h"
#include "TF3.h"
#include "TH2.h"
#include "TH3.h"
#include "TVector3.h"

namespace sct {
class Nucleus {
public:
  // default constructor initializes to an invalid state,
  // to avoid ambiguity
  Nucleus();

  // generic constructor
  Nucleus(unsigned mass_number, // number of nucleons
          double radius,        // radius of the nucleus (in fm)
          double skin_depth,    // skin depth of the nucleus (in fm)
          double beta2,         // 2nd order deformation parameter
          double beta4);        // 4th order deformation parameter

  Nucleus(const Nucleus &rhs);

  virtual ~Nucleus();

  // clears the current set of generated nucleons
  void clear();

  // allows the user to (re)set the parameters after construction, if resetting,
  // will delete any currently generated event
  void setParameters(unsigned mass_number, double radius, double skin_depth,
                     double beta2, double beta4);

  // creates a new set of nucleons according to the nucleus parameters,
  // where b is the impact parameter of the generated event. If b = 0,
  // the center of the nucleus sits at (0, 0, 0), otherwise, the nucleus
  // center will sit at (b, 0, 0)
  bool generate(double b = 0.0);

  // turn on nucleon smearing away from the Woods-Saxon profile, using either
  // a hard core (flat probability with radius sigmaNN)  or gaussian profile
  // (with sigma of 0.79/sqrt(3), and max of sigmaNN * 5)
  void setNucleonSmearing(NucleonSmearing smear, double sigmaNN);

  // sets a minimum distance between generated nucleons
  void setRepulsionDistance(double fm);

  // set string identifier for nucleus
  void setName(string name) { name_ = name; }
  string name() const { return name_; }

  // if set to true, the nuclei will be rotated randomly along polar & azimuthal
  // angles. Default is set to true, but only used if the nucleus is not
  // spherical (ie, beta2 or beta4 != 0.0). The option to turn this off is
  // really only included for testing, since it needs to be on during glauber
  // simulations.
  void setRandomOrientation(bool flag) { random_orientation_ = flag; }
  inline bool randomOrientation() const { return random_orientation_; }

  // access to nucleons
  const Nucleon &operator[](unsigned idx) const;
  Nucleon &operator[](unsigned idx);

  // access to parameters
  inline unsigned massNumber() const { return mass_number_; }
  inline unsigned size() const { return nucleons_.size(); }
  inline double radius() const { return radius_; }
  inline double skinDepth() const { return skin_depth_; }
  inline double beta2() const { return beta2_; }
  inline double beta4() const { return beta4_; }
  inline NucleonSmearing nucleonSmearing() const { return smear_; }
  inline double nnCrossSection() const { return sigmaNN_; }
  inline double repulsionDistance() const { return repulsion_distance_; }
  inline double nucleusTheta() const { return nucleus_theta_; }
  inline double nucleusPhi() const { return nucleus_phi_; }
  inline double impactParameter() const { return b_; }

  inline TH2D *generatedRCosTheta() const {
    return (TH2D *)generated_rcos_theta_.get();
  }
  inline TH3D *generatedPosition() const {
    return (TH3D *)generated_position_.get();
  }
  inline TH3D *generatedSmear() const { return (TH3D *)generated_smear_.get(); }
  inline TH3D *smearedPosition() const {
    return (TH3D *)smeared_position_.get();
  }

  inline TF1 *woodsSaxon() const { return woods_saxon_.get(); }

private:
  // generates one random nucleon and rotates it wrt the nucleus orientation, if
  // applicable
  void addNucleon(double b);

  // (re)creates the Woods-Saxon density profile if parameters are changed
  void initWoodsSaxon();

  // (re)creates QA histograms
  void initHistograms();

  // functions to generate a random nucleon position, and to generate a smearing
  // factor
  TVector3 generateNucleonPosition();
  TVector3 smear();

  std::vector<Nucleon> nucleons_; // container for nucleons

  string name_; // string identifier

  unsigned mass_number_; // nuclear mass number
  double radius_;        // radius of nucleus
  double skin_depth_;    // skin depth of nucleus
  double beta2_;         // 2nd order deformation parameter
  double beta4_;         // 4th order deformation parameter

  NucleonSmearing smear_;     // flag for nucleon position smearing
  double sigmaNN_;            // nucleon-nucleon inellastic cross section
  double repulsion_distance_; // force nucleons to be minimum
                              // repulsionDistance_ away from each other

  bool random_orientation_; // if set to true, the nucleus will be oriented in
                            // a random direction (only useful if beta2 or
                            // beta4 are non-zero)
  double nucleus_theta_; // for deformed nuclei, specifies the polar & azimuthal
  double nucleus_phi_;   // angles in the collision frame for a specific event
  double b_;             // impact parameter for a specific event

  unique_ptr<TF1> woods_saxon_; // density profile for nucleons
  unique_ptr<TF3>
      smearing_profile_; // Used to smear nucleon position if requested

  unique_ptr<TH2> generated_rcos_theta_; // histogram recording sampled r/theta
                                         // from woods-saxon
  unique_ptr<TH3> generated_position_;   // histogram recording the generated
                                         // nucleon position
  unique_ptr<TH3>
      generated_smear_; // histogram recording the generated smearing
  unique_ptr<TH3>
      smeared_position_; // histogram recording the final smeared position
};
} // namespace sct

#endif // SCT_GLAUBER_NUCLEUS_H