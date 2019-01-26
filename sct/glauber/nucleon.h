#ifndef SCT_GLAUBER_NUCLEON_H
#define SCT_GLAUBER_NUCLEON_H

#include "TVector3.h"

namespace sct {
  
  class Nucleon {
  public:
    
    Nucleon();
    Nucleon(const Nucleon& rhs);
    virtual ~Nucleon();
    
    // resets to default
    void clear();
    
    // sets r, theta, phi, and rotates wrt the orientation
    // of the nucleus (given by nucleusTheta, nucleusPhi),
    // and modifies position along X given the impact parameter
    // rotation is done along theta axis first if thetaFirst is set
    // to true
    void set(TVector3 position, double b = 0, double nucleusTheta = 0.0,
             double nucleusPhi = 0.0, bool thetaFirst = true);
    
    inline TVector3 position() const  {return position_;}
    inline double phi() const         {return position_.Phi();}
    inline double theta() const       {return position_.Theta();}
    inline double r() const           {return position_.Mag();}
    
    inline double x() const {return position_.X();}
    inline double y() const {return position_.Y();}
    inline double z() const {return position_.Z();}
    inline double x2() const {return position_.X() * position_.X();}
    inline double y2() const {return position_.Y() * position_.Y();}
    inline double z2() const {return position_.Z() * position_.Z();}
    inline double xy() const {return position_.X() * position_.Y();}
    inline double yz() const {return position_.Y() * position_.Z();}
    inline double zx() const {return position_.Z() * position_.X();}
    inline double x3() const {return position_.X() * position_.X() * position_.X();}
    inline double y3() const {return position_.Y() * position_.Y() * position_.Y();}
    inline double z3() const {return position_.Z() * position_.Z() * position_.Z();}
    inline double x2y() const {return position_.X() * position_.X() * position_.Y();}
    inline double xy2() const {return position_.X() * position_.Y() * position_.Y();}
    inline double y2z() const {return position_.Y() * position_.Y() * position_.Z();}
    inline double yz2() const {return position_.Y() * position_.Z() * position_.Z();}
    inline double z2x() const {return position_.Z() * position_.Z() * position_.X();}
    inline double zx2() const {return position_.Z() * position_.X() * position_.X();}
    
    // three dimensional delta R
    double deltaR(const Nucleon& rhs) const;
    
    // two dimensional delta R (in XY)
    double deltaXY(const Nucleon& rhs) const;
    
    inline bool     participant() const  {return n_coll_ ? true : false;}
    inline unsigned nColl() const        {return n_coll_;}
    inline double   multiplicity() const {return multiplicity_;}
    
    inline void incrementNColl() {++n_coll_;}
    inline void setMultiplicity(double mult) {multiplicity_ = mult;}
    
  private:
    TVector3 position_;
    unsigned n_coll_;
    unsigned multiplicity_;
  };
  
} // namespace sct

#endif // SCT_GLAUBER_NUCLEON_H
