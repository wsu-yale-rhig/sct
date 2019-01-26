#ifndef SCT_GLAUBER_COLLISION_H
#define SCT_GLAUBER_COLLISION_H

/* Calculates number of binary collisions, number
 * of participants, average event quantities, and
 * eccentricities, given two nuclei.
 */

#include <array>

#include "sct/lib/enumerations.h"
#include "sct/lib/memory.h"
#include "sct/glauber/nucleus.h"
#include "sct/glauber/nucleon.h"
#include "sct/centrality/multiplicity_model.h"

namespace sct {
  class Collision {
  public:
    Collision();
    ~Collision();
    
    // defines the cross section for a collision between two nucleons
    inline void setNNCrossSection(double NNxSec) {inelasticXSec_ = NNxSec;}
    inline double NNCrossSection() const {return inelasticXSec_;}
    
    // defines the method for calculating nucleon-nucleon collisions:
    // 1) CollisionProfile::HardCore (default)
    //       if the two nucleons are less than a distance R = sqrt(NNCrossSection/pi)
    //       a collision occurs
    // 2) CollisionProfile::Gaussian
    //       gaussian probability distribution for collision: exp(-x^2/(2*sig^2))
    //       x = distance between nucleons in XY plane
    //       sig = sqrt(NNCrossSection/pi)
    inline void setCollisionProfile(CollisionProfile profile) {profile_ = profile;}
    inline CollisionProfile collisionProfile() const {return profile_;}
    
    // estimates a multiplicity for each binary collision, allowing
    // observables to be calculated weighted by the estimated multiplicity
    // generally, these parameters come from a fit of the glauber model to data,
    // so it may be required to run the model twice. Once to generate
    // Npart Ncoll pairs for fitting, and again once the multiplicity model
    // has been optimized
    void setMultiplicityModel(double npp, double k, double x, double ppEff,
                              double aaEff, double aaCent, double trigEff = 1.0,
                              bool constEff = false);

    // collide two lists of nucleons - return true if at least one
    // nucleon-nucleon collision occurred. This is instantiated for both
    // Nucleus & std::vector<Nucleon> in the library. For other containers,
    // you'll need to instantiate and recompile :) 
    template<typename Container>
    bool collide(Container& nucleusA, Container& nucleusB);
    
    // checks if two nucleons collide using the specified collision profile
    bool nucleonCollision(Nucleon nucleonA, Nucleon nucleonB);
    
    // zeros collision statistics
    void clear();
    
    // returns number of binary collisions
    inline unsigned nColl() const {return nColl_;}
    inline unsigned nPart() const {return count_[0];}
    inline unsigned spectators() const {return count_[2];}
    
    // access to arrays
    std::array<unsigned, nGlauberWeights>& countArray() {return count_;}
    std::array<double, nGlauberWeights>& averageX() {return avgX_;}
    std::array<double, nGlauberWeights>& averageY() {return avgY_;}
    std::array<double, nGlauberWeights>& averageX2() {return avgX2_;}
    std::array<double, nGlauberWeights>& averageY2() {return avgY2_;}
    std::array<double, nGlauberWeights>& averageXY() {return avgXY_;}
    std::array<double, nGlauberWeights>& reactionPlane2Ecc() {return eccRP2_;}
    std::array<double, nGlauberWeights>& partPlane2Ecc() {return eccPP2_;}
    std::array<double, nGlauberWeights>& partPlane3Ecc() {return eccPP3_;}
    std::array<double, nGlauberWeights>& partPlane4Ecc() {return eccPP4_;}
    std::array<double, nGlauberWeights>& partPlane2() {return pp2_;}
    std::array<double, nGlauberWeights>& partPlane3() {return pp3_;}
    std::array<double, nGlauberWeights>& partPlane4() {return pp4_;}
    
  private:
    
    // used by collide() to calculate participant plane eccentricity
    template<typename Container>
    std::pair<double, double> participantPlaneEcc(Container& nucleusA,
                                                  Container& nucleusB,
                                                  int order, int weightIdx);
    
    double inelasticXSec_;        // nucleon-nucleon cross section
    
    CollisionProfile profile_;    // either normal hard-core collision profile
                                  // or gaussian profile
    
    // multiplicity model
    unique_ptr<MultiplicityModel> mult_model_;

    unsigned nColl_;
    std::array<unsigned, nGlauberWeights> count_;
    std::array<double, nGlauberWeights> avgX_;
    std::array<double, nGlauberWeights> avgY_;
    std::array<double, nGlauberWeights> avgX2_;
    std::array<double, nGlauberWeights> avgY2_;
    std::array<double, nGlauberWeights> avgXY_;
    std::array<double, nGlauberWeights> eccRP2_;
    std::array<double, nGlauberWeights> eccPP2_;
    std::array<double, nGlauberWeights> eccPP3_;
    std::array<double, nGlauberWeights> eccPP4_;
    std::array<double, nGlauberWeights> pp2_;
    std::array<double, nGlauberWeights> pp3_;
    std::array<double, nGlauberWeights> pp4_;
  };
} // namespace sct

#endif // SCT_GLAUBER_COLLISION_H
