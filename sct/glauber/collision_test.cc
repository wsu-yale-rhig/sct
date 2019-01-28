#include "sct/glauber/collision.h"
#include "sct/glauber/nucleon.h"
#include "sct/glauber/nucleus.h"
#include "sct/lib/enumerations.h"
#include "sct/lib/logging.h"
#include "sct/lib/math.h"

#include <vector>

#include "gtest/gtest.h"

// testing that parameters are properly set
TEST(Collision, initialization) {
  sct::Collision collision;
  collision.setCollisionProfile(sct::CollisionProfile::HardCore);
  collision.setNNCrossSection(1.0);

  EXPECT_EQ(collision.collisionProfile(), sct::CollisionProfile::HardCore);
  EXPECT_EQ(collision.NNCrossSection(), 1.0);

  collision.setCollisionProfile(sct::CollisionProfile::Gaussian);
  collision.setNNCrossSection(2.0);

  EXPECT_EQ(collision.collisionProfile(), sct::CollisionProfile::Gaussian);
  EXPECT_EQ(collision.NNCrossSection(), 2.0);
}

// test that a single proton-proton collision occurs as expected
TEST(Collision, singleProtonCollision) {
  sct::Collision collision;
  collision.setCollisionProfile(sct::CollisionProfile::HardCore);
  collision.setNNCrossSection(1.0 * sct::pi);

  sct::Nucleon nucleonA;
  nucleonA.set(TVector3(0, 0, 0));
  sct::Nucleon nucleonB;
  nucleonB.set(TVector3(0, 0, 0));

  std::vector<sct::Nucleon> nucleusA;
  nucleusA.push_back(nucleonA);
  std::vector<sct::Nucleon> nucleusB;
  nucleusB.push_back(nucleonB);

  collision.collide(nucleusA, nucleusB);

  EXPECT_EQ(collision.nColl(), 1);
  EXPECT_EQ(nucleusA[0].nColl(), 1);
  EXPECT_EQ(nucleusB[0].nColl(), 1);
}

// test that multiple proton-proton collisions occur as expected
TEST(Collision, multipleProtonCollision) {
  sct::Collision collision;
  collision.setCollisionProfile(sct::CollisionProfile::HardCore);
  collision.setNNCrossSection(1.0 * sct::pi);

  sct::Nucleon nucleonA;
  nucleonA.set(TVector3(0, 0, 0));
  sct::Nucleon nucleonB;
  nucleonB.set(TVector3(0, 0, 0));

  std::vector<sct::Nucleon> nucleusA;
  nucleusA.push_back(nucleonA);
  nucleusA.push_back(nucleonB);
  std::vector<sct::Nucleon> nucleusB;
  nucleusB.push_back(nucleonB);
  nucleusB.push_back(nucleonB);

  collision.collide(nucleusA, nucleusB);

  EXPECT_EQ(collision.nColl(), 4);
  EXPECT_EQ(nucleusA[0].nColl(), 2);
  EXPECT_EQ(nucleusB[0].nColl(), 2);
  EXPECT_EQ(nucleusA[1].nColl(), 2);
  EXPECT_EQ(nucleusB[1].nColl(), 2);
}

// test that a single proton-proton miss occurs as expected
TEST(Collision, singleProtonMiss) {
  sct::Collision collision;
  collision.setCollisionProfile(sct::CollisionProfile::HardCore);
  collision.setNNCrossSection(1.0 * sct::pi);

  sct::Nucleon nucleonA;
  nucleonA.set(TVector3(0, 0, 0));
  sct::Nucleon nucleonB;
  nucleonB.set(TVector3(1.01, 0, 0));

  std::vector<sct::Nucleon> nucleusA;
  nucleusA.push_back(nucleonA);
  std::vector<sct::Nucleon> nucleusB;
  nucleusB.push_back(nucleonB);

  collision.collide(nucleusA, nucleusB);

  EXPECT_EQ(collision.nColl(), 0);
  EXPECT_EQ(nucleusA[0].nColl(), 0);
  EXPECT_EQ(nucleusB[0].nColl(), 0);
}

// test that a single proton-proton produces correct event averages
TEST(Collision, singleProtonCollisionAverages) {
  sct::Collision collision;
  collision.setCollisionProfile(sct::CollisionProfile::HardCore);
  collision.setNNCrossSection(1.0 * sct::pi);

  sct::Nucleon nucleonA;
  nucleonA.set(TVector3(0.5, 0, 0));
  sct::Nucleon nucleonB;
  nucleonB.set(TVector3(1, 0, 0));

  std::vector<sct::Nucleon> nucleusA;
  nucleusA.push_back(nucleonA);
  std::vector<sct::Nucleon> nucleusB;
  nucleusB.push_back(nucleonB);

  collision.collide(nucleusA, nucleusB);

  EXPECT_NEAR(collision.averageX()[0], 0.75, 1e-8);
  EXPECT_NEAR(collision.averageX()[1], 0.75, 1e-8);
  EXPECT_NEAR(collision.averageX()[2], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageY()[0], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageY()[1], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageY()[2], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageX2()[0], 0.625, 1e-8);
  EXPECT_NEAR(collision.averageX2()[1], 0.625, 1e-8);
  EXPECT_NEAR(collision.averageX2()[2], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageY2()[0], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageY2()[1], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageY2()[2], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageXY()[0], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageXY()[1], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageXY()[2], 0.0, 1e-8);
}

// test that a single proton-proton produces correct event averages,
// rotated to a new direction
TEST(Collision, singleProtonCollisionAveragesRotated) {
  sct::Collision collision;
  collision.setCollisionProfile(sct::CollisionProfile::HardCore);
  collision.setNNCrossSection(1.0 * sct::pi);

  sct::Nucleon nucleonA;
  nucleonA.set(TVector3(0.5, 0.5, 0));
  sct::Nucleon nucleonB;
  nucleonB.set(TVector3(1.0, 1.0, 0));

  std::vector<sct::Nucleon> nucleusA;
  nucleusA.push_back(nucleonA);
  std::vector<sct::Nucleon> nucleusB;
  nucleusB.push_back(nucleonB);

  collision.collide(nucleusA, nucleusB);

  EXPECT_NEAR(collision.averageX()[0], 0.75, 1e-8);
  EXPECT_NEAR(collision.averageX()[1], 0.75, 1e-8);
  EXPECT_NEAR(collision.averageX()[2], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageY()[0], 0.75, 1e-8);
  EXPECT_NEAR(collision.averageY()[1], 0.75, 1e-8);
  EXPECT_NEAR(collision.averageY()[2], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageX2()[0], 0.625, 1e-8);
  EXPECT_NEAR(collision.averageX2()[1], 0.625, 1e-8);
  EXPECT_NEAR(collision.averageX2()[2], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageX2()[0], 0.625, 1e-8);
  EXPECT_NEAR(collision.averageX2()[1], 0.625, 1e-8);
  EXPECT_NEAR(collision.averageY2()[2], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageXY()[0], 0.625, 1e-8);
  EXPECT_NEAR(collision.averageXY()[1], 0.625, 1e-8);
  EXPECT_NEAR(collision.averageXY()[2], 0.0, 1e-8);
}

// test that a single proton-proton produces correct event averages,
// rotated to a new direction
TEST(Collision, singleProtonCollisionAveragesRotated2) {
  sct::Collision collision;
  collision.setCollisionProfile(sct::CollisionProfile::HardCore);
  collision.setNNCrossSection(1.0 * sct::pi);

  sct::Nucleon nucleonA;
  nucleonA.set(TVector3(0.0, 1.0, 0.0));
  sct::Nucleon nucleonB;
  nucleonB.set(TVector3(0.0, 0.5, 0.0));

  std::vector<sct::Nucleon> nucleusA;
  nucleusA.push_back(nucleonA);
  std::vector<sct::Nucleon> nucleusB;
  nucleusB.push_back(nucleonB);

  collision.collide(nucleusA, nucleusB);

  EXPECT_NEAR(collision.averageX()[0], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageX()[1], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageX()[2], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageY()[0], 0.75, 1e-8);
  EXPECT_NEAR(collision.averageY()[1], 0.75, 1e-8);
  EXPECT_NEAR(collision.averageY()[2], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageX2()[0], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageX2()[1], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageX2()[2], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageY2()[0], 0.625, 1e-8);
  EXPECT_NEAR(collision.averageY2()[1], 0.625, 1e-8);
  EXPECT_NEAR(collision.averageY2()[2], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageXY()[0], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageXY()[1], 0.0, 1e-8);
  EXPECT_NEAR(collision.averageXY()[2], 0.0, 1e-8);
}

// test that the participant plane is calculated correctly
TEST(Collision, participantPlane2) {
  sct::Collision collision;
  collision.setCollisionProfile(sct::CollisionProfile::HardCore);
  collision.setNNCrossSection(2.0 * sct::pi);

  double x = 0.2;
  double y = 0.2;
  double h = sqrt(x * x + y * y);

  sct::Nucleon nucleonA;
  nucleonA.set(TVector3(x, -y, 0.0));
  sct::Nucleon nucleonB;
  nucleonB.set(TVector3(0.0, 0.0, 0.0));
  sct::Nucleon nucleonC;
  nucleonC.set(TVector3(-x, y, 0.0));

  std::vector<sct::Nucleon> nucleusA;
  nucleusA.push_back(nucleonA);
  std::vector<sct::Nucleon> nucleusB;
  nucleusB.push_back(nucleonB);

  collision.collide(nucleusA, nucleusB);

  EXPECT_NEAR(-999, collision.partPlane2()[0], 1e-8);

  nucleusA.push_back(nucleonC);

  collision.collide(nucleusA, nucleusB);

  double phi1 = 3.0 / 4.0 * sct::pi;
  double phi2 = phi1 + sct::pi;
  double weight = h * h;

  double qx = weight * (cos(2.0 * phi1) + cos(2.0 * phi2));
  double qy = weight * (sin(2.0 * phi1) + sin(2.0 * phi2));

  EXPECT_NEAR(atan2(qy, -qx), collision.partPlane2()[0], 1e-8);
}

// test that the participant plane is calculated correctly, rotated to
// a new frame
TEST(Collision, participantPlane2Rotated) {
  sct::Collision collision;
  collision.setCollisionProfile(sct::CollisionProfile::HardCore);
  collision.setNNCrossSection(2.0 * sct::pi);

  double x = 0.1;
  double y = 0.15;
  double h = sqrt(x * x + y * y);

  sct::Nucleon nucleonA;
  nucleonA.set(TVector3(-x - 0.06, 0, 0.0));
  sct::Nucleon nucleonB;
  nucleonB.set(TVector3(0.0, y, 0.0));
  sct::Nucleon nucleonC;
  nucleonC.set(TVector3(x, 0, 0.0));

  std::vector<sct::Nucleon> nucleusA;
  nucleusA.push_back(nucleonA);
  std::vector<sct::Nucleon> nucleusB;
  nucleusB.push_back(nucleonC);

  collision.collide(nucleusA, nucleusB);

  EXPECT_NEAR(-999, collision.partPlane2()[0], 1e-8);

  nucleusA.push_back(nucleonB);

  collision.collide(nucleusA, nucleusB);

  double avg_y = y / 3.0;
  double avg_x = (0.1 - 0.1 - 0.06 + 0) / 3.0;

  TVector3 v1(-x - 0.06 - avg_x, -avg_y, 0);
  TVector3 v2(0.0 - avg_x, y - avg_y, 0);
  TVector3 v3(x - avg_x, -avg_y, 0);

  double qx = v1.Perp() * v1.Perp() * cos(3.0 * v1.Phi()) +
              v2.Perp() * v2.Perp() * cos(3.0 * v2.Phi()) +
              v3.Perp() * v3.Perp() * cos(3.0 * v3.Phi());
  double qy = v1.Perp() * v1.Perp() * sin(3.0 * v1.Phi()) +
              v2.Perp() * v2.Perp() * sin(3.0 * v2.Phi()) +
              v3.Perp() * v3.Perp() * sin(3.0 * v3.Phi());

  EXPECT_NEAR(atan2(qy / 3.0, -qx / 3.0), collision.partPlane3()[0], 1e-8);
}

// test that the participant plane eccentricity is calculated correctly
// for a circular distribution
TEST(Collision, participantPlaneEccentricity2) {
  sct::Collision collision;
  collision.setCollisionProfile(sct::CollisionProfile::HardCore);
  collision.setNNCrossSection(2.0 * sct::pi);

  double x = 0.1;
  double y = 0.1;

  sct::Nucleon nucleonA;
  nucleonA.set(TVector3(-x, 0, 0.0));
  sct::Nucleon nucleonB;
  nucleonB.set(TVector3(x, 0.0, 0.0));
  sct::Nucleon nucleonC;
  nucleonC.set(TVector3(0.0, y, 0.0));
  sct::Nucleon nucleonD;
  nucleonD.set(TVector3(0.0, -y, 0.0));

  std::vector<sct::Nucleon> nucleusA;
  nucleusA.push_back(nucleonA);
  nucleusA.push_back(nucleonB);
  std::vector<sct::Nucleon> nucleusB;
  nucleusB.push_back(nucleonC);
  nucleusB.push_back(nucleonD);

  collision.collide(nucleusA, nucleusB);

  EXPECT_NEAR(0.0, collision.partPlane2Ecc()[0], 1e-8);
}

// test that the participant plane eccentricity is calculated correctly
// for an elongated distribution
TEST(Collision, reactionPlaneEccentricity2) {
  sct::Collision collision;
  collision.setCollisionProfile(sct::CollisionProfile::HardCore);
  collision.setNNCrossSection(2.0 * sct::pi);

  double x = 0.1;
  double y = 0.1;

  sct::Nucleon nucleonA;
  nucleonA.set(TVector3(-x, 0, 0.0));
  sct::Nucleon nucleonB;
  nucleonB.set(TVector3(x, 0.0, 0.0));
  sct::Nucleon nucleonC;
  nucleonC.set(TVector3(0.0, y, 0.0));
  sct::Nucleon nucleonD;
  nucleonD.set(TVector3(0.0, -y, 0.0));

  std::vector<sct::Nucleon> nucleusA;
  nucleusA.push_back(nucleonA);
  nucleusA.push_back(nucleonB);
  std::vector<sct::Nucleon> nucleusB;
  nucleusB.push_back(nucleonC);
  nucleusB.push_back(nucleonD);

  collision.collide(nucleusA, nucleusB);

  EXPECT_NEAR(0.0, collision.reactionPlane2Ecc()[0], 1e-8);

  y *= 2;

  nucleonC.set(TVector3(0.0, y, 0.0));
  nucleonD.set(TVector3(0.0, y, 0.0));

  nucleusB[0] = nucleonC;
  nucleusB[1] = nucleonD;

  collision.collide(nucleusA, nucleusB);

  // calculate averages
  double avgX = (nucleonA.x() + nucleonB.x()) / 2.0;
  double avgY = (nucleonC.y() + nucleonD.y()) / 2.0;
  double avgX2 = (pow(nucleonA.x(), 2.0) + pow(nucleonB.x(), 2.0));
  double avgY2 = (pow(nucleonC.y(), 2.0) + pow(nucleonD.y(), 2.0));
  double sigmaX2 = avgX2 - pow(avgX, 2.0);
  double sigmaY2 = avgY2 - pow(avgY, 2.0);

  // calculate eccentricity
  double e = (sigmaY2 - sigmaX2) / (sigmaX2 + sigmaY2);

  EXPECT_NEAR(e, collision.reactionPlane2Ecc()[0], 1e-8);
}

// test that a more complex system gives proper nColl, nPart, and spectators
TEST(Collision, multipleProtonCollisionSpectators) {
  sct::Collision collision;
  collision.setCollisionProfile(sct::CollisionProfile::HardCore);
  collision.setNNCrossSection(1.0 * sct::pi);

  sct::Nucleon nucleon1;
  sct::Nucleon nucleon2;
  sct::Nucleon nucleon3;
  sct::Nucleon nucleon4;
  sct::Nucleon nucleon5;
  sct::Nucleon nucleon6;

  nucleon1.set(TVector3(0.0, 0.0, 0.0));
  nucleon2.set(TVector3(0.0, 1.5, 0.0));
  nucleon3.set(TVector3(1.5, 0.0, 0.0));
  nucleon4.set(TVector3(0.0, 1.0, 0.0));
  nucleon5.set(TVector3(0.0, -1.0, 0.0));
  nucleon6.set(TVector3(0.0, 8.0, 0.0));

  std::vector<sct::Nucleon> nucleusA = {nucleon1, nucleon2, nucleon3};
  std::vector<sct::Nucleon> nucleusB = {nucleon4, nucleon5, nucleon6};

  collision.collide(nucleusA, nucleusB);

  EXPECT_EQ(collision.nColl(), 3);
  EXPECT_EQ(collision.countArray()[2], 2);
  EXPECT_EQ(collision.countArray()[1] / 2, 3);
  EXPECT_EQ(collision.countArray()[0], 4);
}
