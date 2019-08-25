#include "sct/glauber/nucleon.h"
#include "sct/lib/enumerations.h"
#include "sct/lib/logging.h"
#include "sct/lib/math.h"
#include "sct/utils/random.h"

#include <random>

#include "gtest/gtest.h"

// first, a few specific examples of increasing complexity
// to make sure I understand how its working
TEST(nucleon, specificExample1) {
  sct::Nucleon nucleon;

  // define a position
  double x = 2.0;
  double y = 0.0;
  double z = 0.0;

  // translate to r, theta phi
  double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  double theta = acos(z / r);
  double phi = atan2(y, x);

  TVector3 pos;
  pos.SetMagThetaPhi(r, theta, phi);

  nucleon.set(pos);

  // check r, theta, phi, x, y, z
  EXPECT_NEAR(2.0, nucleon.x(), 1e-8);
  EXPECT_NEAR(0, nucleon.y(), 1e-8);
  EXPECT_NEAR(0, nucleon.z(), 1e-8);
  EXPECT_NEAR(2.0, nucleon.r(), 1e-8);
  EXPECT_NEAR(sct::pi / 2.0, nucleon.theta(), 1e-8);
  EXPECT_NEAR(0.0, nucleon.phi(), 1e-8);
}

// Now we'll do random cases as a simple fuzzing

TEST(nucleon, position) {
  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_real_distribution<> dis(-2.0, 2.0);

  sct::Nucleon nucleon;
  for (int i = 0; i < 1e7; ++i) {
    // define a position
    double x = dis(generator);
    double y = dis(generator);
    double z = dis(generator);

    // translate to r, theta phi
    double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    double theta = acos(z / r);
    double phi = atan2(y, x);

    TVector3 pos;
    pos.SetMagThetaPhi(r, theta, phi);

    nucleon.set(pos);

    // check r, theta, phi, x, y, z
    EXPECT_NEAR(x, nucleon.x(), 1e-8);
    EXPECT_NEAR(y, nucleon.y(), 1e-8);
    EXPECT_NEAR(z, nucleon.z(), 1e-8);
    EXPECT_NEAR(r, nucleon.r(), 1e-8);
    EXPECT_NEAR(theta, nucleon.theta(), 1e-8);
    EXPECT_NEAR(phi, nucleon.phi(), 1e-8);
  }
}

// test position wrt the angles
TEST(nucleon, angles) {
  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_real_distribution<> dis(-2.0, 2.0);

  sct::Nucleon nucleon;
  for (int i = 0; i < 1e7; ++i) {
    // define a position
    double x = dis(generator);
    double y = dis(generator);
    double z = dis(generator);

    // translate to r, theta phi
    double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    double theta = acos(z / r);
    double phi = atan2(y, x);

    TVector3 pos;
    pos.SetMagThetaPhi(r, theta, phi);

    nucleon.set(pos);

    // check r, theta, phi, x, y, z
    EXPECT_NEAR(x, nucleon.x(),
                1e-8);  // impact parameter is always a shift along X
    EXPECT_NEAR(y, nucleon.y(), 1e-8);
    EXPECT_NEAR(z, nucleon.z(), 1e-8);

    // now we have to recalculate r, theta, phi
    double rprime = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    double thetaprime = acos(z / rprime);
    double phiprime = atan2(y, x);

    EXPECT_NEAR(rprime, nucleon.r(), 1e-8);
    EXPECT_NEAR(thetaprime, nucleon.theta(), 1e-8);
    EXPECT_NEAR(phiprime, nucleon.phi(), 1e-8);
  }
}

TEST(nucleon, deltaR) {
  sct::Nucleon nucleonA;
  sct::Nucleon nucleonB;

  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_real_distribution<> dis(-2.0, 2.0);

  for (int i = 0; i < 1e7; ++i) {
    // define a position
    double x_a = dis(generator);
    double y_a = dis(generator);
    double z_a = dis(generator);
    double x_b = dis(generator);
    double y_b = dis(generator);
    double z_b = dis(generator);

    // translate to r, theta phi
    double r_a = sqrt(pow(x_a, 2) + pow(y_a, 2) + pow(z_a, 2));
    double theta_a = acos(z_a / r_a);
    double phi_a = atan2(y_a, x_a);
    double r_b = sqrt(pow(x_b, 2) + pow(y_b, 2) + pow(z_b, 2));
    double theta_b = acos(z_b / r_b);
    double phi_b = atan2(y_b, x_b);

    TVector3 pos_a;
    pos_a.SetMagThetaPhi(r_a, theta_a, phi_a);

    TVector3 pos_b;
    pos_b.SetMagThetaPhi(r_b, theta_b, phi_b);

    nucleonA.set(pos_a);
    nucleonB.set(pos_b);

    double delta_R =
        sqrt(pow(x_a - x_b, 2) + pow(y_a - y_b, 2) + pow(z_a - z_b, 2));
    EXPECT_NEAR(nucleonA.deltaR(nucleonB), delta_R, 1e-8);
  }
}

TEST(nucleon, specificDeltaR) {
  sct::Nucleon nucleonA;
  sct::Nucleon nucleonB;

  double x_a = 2.0;
  double y_a = 0.0;
  double z_a = 0.0;
  double x_b = 1.0;
  double y_b = 0.0;
  double z_b = 0.0;

  // translate to r, theta phi
  double r_a = sqrt(pow(x_a, 2) + pow(y_a, 2) + pow(z_a, 2));
  double theta_a = acos(z_a / r_a);
  double phi_a = atan2(y_a, x_a);
  double r_b = sqrt(pow(x_b, 2) + pow(y_b, 2) + pow(z_b, 2));
  double theta_b = acos(z_b / r_b);
  double phi_b = atan2(y_b, x_b);

  TVector3 pos_a;
  pos_a.SetMagThetaPhi(r_a, theta_a, phi_a);

  TVector3 pos_b;
  pos_b.SetMagThetaPhi(r_b, theta_b, phi_b);

  nucleonA.set(pos_a);
  nucleonB.set(pos_b);

  EXPECT_NEAR(nucleonA.deltaR(nucleonB), 1.0, 1e-8);
}
