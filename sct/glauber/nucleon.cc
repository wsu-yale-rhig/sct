#include "sct/glauber/nucleon.h"

#include "sct/lib/logging.h"

namespace sct {

Nucleon::Nucleon() : position_(0, 0, 0), n_coll_(0), multiplicity_(0.0) {}

Nucleon::Nucleon(const Nucleon& rhs)
    : position_(rhs.position()),
      n_coll_(rhs.nColl()),
      multiplicity_(rhs.multiplicity()) {}

Nucleon::~Nucleon() {}

void Nucleon::clear() {
  position_.SetXYZ(0, 0, 0);
  n_coll_ = 0;
  multiplicity_ = 0;
}

void Nucleon::set(TVector3 position, double b, double nucleus_theta,
                  double nucleus_phi, bool theta_first) {
  position_ = position;

  // rotate the nucleus
  if (theta_first) {
    position_.RotateY(nucleus_theta);
    position_.RotateZ(nucleus_phi);
  } else {
    position_.RotateZ(nucleus_phi);
    position_.RotateY(nucleus_theta);
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

}  // namespace sct
