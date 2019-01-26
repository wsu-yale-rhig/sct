#include "sct/core/nucleon.h"

#include "sct/core/base.h"
#include "sct/core/logging.h"

namespace sct {
  Nucleon::Nucleon()
  : position_(0,0,0),
    n_coll_(0),
    multiplicity_(0.0) {}
  
  Nucleon::Nucleon(const Nucleon& rhs)
  : position_(rhs.position()),
    n_coll_(rhs.nColl()),
    multiplicity_(rhs.multiplicity()) { }

  Nucleon::~Nucleon() {}

  void Nucleon::clear() {
    position_.SetXYZ(0, 0, 0);
    n_coll_ = 0;
    multiplicity_ = 0;
  }

  void Nucleon::set(TVector3 position, double b, double nucleusTheta,
                    double nucleusPhi, bool thetaFirst) {
    position_ = position;
  
    // rotate the nucleus 
    if (thetaFirst) {
      position_.RotateY(nucleusTheta);
      position_.RotateZ(nucleusPhi);
    }
    else {
      position_.RotateZ(nucleusPhi);
      position_.RotateY(nucleusTheta);
    }
  
    // offset by impact parameter
    position_.SetX(position_.X() + b);
    n_coll_ = 0;
    multiplicity_ = 0.0;
  }
  
  double Nucleon::deltaR(const Nucleon& rhs) const {
    double dX = x() - rhs.x();
    double dY = y() - rhs.y();
    double dZ = z() - rhs.z();
    return sqrt(pow(dX, 2) + pow(dY, 2) + pow(dZ, 2));
  }
  
  double Nucleon::deltaXY(const Nucleon& rhs) const {
    double dX = x() - rhs.x();
    double dY = y() - rhs.y();
    return sqrt(pow(dX, 2) + pow(dY, 2));
  }
} // namespace sct
