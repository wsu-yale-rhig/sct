#include "sct/utils/nucleus_info.h"
#include "sct/lib/enumerations.h"

#include "gtest/gtest.h"

TEST(NucleusInfo, Au197) {
  sct::NucleusInfo& info = sct::NucleusInfo::instance();

  EXPECT_EQ(info.massNumber(sct::GlauberSpecies::Au197), 197);
  EXPECT_EQ(info.chargeNumber(sct::GlauberSpecies::Au197), 79);
  EXPECT_NEAR(info.radius(sct::GlauberSpecies::Au197), 6.38, 1e-5);
  EXPECT_NEAR(info.skinDepth(sct::GlauberSpecies::Au197), 0.535, 1e-5);
  EXPECT_NEAR(info.radiusError(sct::GlauberSpecies::Au197), 0.12, 1e-5);
  EXPECT_NEAR(info.skinDepthError(sct::GlauberSpecies::Au197), 0.054, 1e-5);
  EXPECT_NEAR(info.beta2(sct::GlauberSpecies::Au197), -0.131, 1e-5);
  EXPECT_NEAR(info.beta4(sct::GlauberSpecies::Au197), -0.031, 1e-5);
  EXPECT_NEAR(info.hulthenA(sct::GlauberSpecies::Au197), 0.0, 1e-5);
  EXPECT_NEAR(info.hulthenB(sct::GlauberSpecies::Au197), 0.0, 1e-5);
  EXPECT_EQ(info.name(sct::GlauberSpecies::Au197), "Au197");
}

TEST(NucleusInfo, Sm154) {
  sct::NucleusInfo& info = sct::NucleusInfo::instance();

  EXPECT_EQ(info.massNumber(sct::GlauberSpecies::Sm154), 154);
  EXPECT_EQ(info.chargeNumber(sct::GlauberSpecies::Sm154), 62);
  EXPECT_NEAR(info.radius(sct::GlauberSpecies::Sm154), 5.9387, 1e-5);
  EXPECT_NEAR(info.skinDepth(sct::GlauberSpecies::Sm154), 0.522, 1e-5);
  EXPECT_NEAR(info.radiusError(sct::GlauberSpecies::Sm154), 0.010, 1e-5);
  EXPECT_NEAR(info.skinDepthError(sct::GlauberSpecies::Sm154), 0.030, 1e-5);
  EXPECT_NEAR(info.beta2(sct::GlauberSpecies::Sm154), 0.270, 1e-5);
  EXPECT_NEAR(info.beta4(sct::GlauberSpecies::Sm154), 0.113, 1e-5);
  EXPECT_NEAR(info.hulthenA(sct::GlauberSpecies::Sm154), 0.0, 1e-5);
  EXPECT_NEAR(info.hulthenB(sct::GlauberSpecies::Sm154), 0.0, 1e-5);
  EXPECT_EQ(info.name(sct::GlauberSpecies::Sm154), "Sm154");
}

TEST(NucleusInfo, U238) {
  sct::NucleusInfo& info = sct::NucleusInfo::instance();

  EXPECT_EQ(info.massNumber(sct::GlauberSpecies::U238), 238);
  EXPECT_EQ(info.chargeNumber(sct::GlauberSpecies::U238), 92);
  EXPECT_NEAR(info.radius(sct::GlauberSpecies::U238), 6.8054, 1e-5);
  EXPECT_NEAR(info.skinDepth(sct::GlauberSpecies::U238), 0.605, 1e-5);
  EXPECT_NEAR(info.radiusError(sct::GlauberSpecies::U238), 0.137, 1e-5);
  EXPECT_NEAR(info.skinDepthError(sct::GlauberSpecies::U238), 0.098, 1e-5);
  EXPECT_NEAR(info.beta2(sct::GlauberSpecies::U238), 0.215, 1e-5);
  EXPECT_NEAR(info.beta4(sct::GlauberSpecies::U238), 0.093, 1e-5);
  EXPECT_NEAR(info.hulthenA(sct::GlauberSpecies::U238), 0.0, 1e-5);
  EXPECT_NEAR(info.hulthenB(sct::GlauberSpecies::U238), 0.0, 1e-5);
  EXPECT_EQ(info.name(sct::GlauberSpecies::U238), "U238");
}

TEST(NucleusInfo, Pb208) {
  sct::NucleusInfo& info = sct::NucleusInfo::instance();

  EXPECT_EQ(info.massNumber(sct::GlauberSpecies::Pb208), 208);
  EXPECT_EQ(info.chargeNumber(sct::GlauberSpecies::Pb208), 82);
  EXPECT_NEAR(info.radius(sct::GlauberSpecies::Pb208), 6.62, 1e-5);
  EXPECT_NEAR(info.skinDepth(sct::GlauberSpecies::Pb208), 0.546, 1e-5);
  EXPECT_NEAR(info.radiusError(sct::GlauberSpecies::Pb208), 0.12, 1e-5);
  EXPECT_NEAR(info.skinDepthError(sct::GlauberSpecies::Pb208), 0.020, 1e-5);
  EXPECT_NEAR(info.beta2(sct::GlauberSpecies::Pb208), 0.00, 1e-5);
  EXPECT_NEAR(info.beta4(sct::GlauberSpecies::Pb208), 0.00, 1e-5);
  EXPECT_NEAR(info.hulthenA(sct::GlauberSpecies::Pb208), 0.0, 1e-5);
  EXPECT_NEAR(info.hulthenB(sct::GlauberSpecies::Pb208), 0.0, 1e-5);
  EXPECT_EQ(info.name(sct::GlauberSpecies::Pb208), "Pb208");
}

TEST(NucleusInfo, Cu63) {
  sct::NucleusInfo& info = sct::NucleusInfo::instance();

  EXPECT_EQ(info.massNumber(sct::GlauberSpecies::Cu63), 63);
  EXPECT_EQ(info.chargeNumber(sct::GlauberSpecies::Cu63), 29);
  EXPECT_NEAR(info.radius(sct::GlauberSpecies::Cu63), 4.218, 1e-5);
  EXPECT_NEAR(info.skinDepth(sct::GlauberSpecies::Cu63), 0.596, 1e-5);
  EXPECT_NEAR(info.radiusError(sct::GlauberSpecies::Cu63), 0.028, 1e-5);
  EXPECT_NEAR(info.skinDepthError(sct::GlauberSpecies::Cu63), 0.01, 1e-5);
  EXPECT_NEAR(info.beta2(sct::GlauberSpecies::Cu63), 0.162, 1e-5);
  EXPECT_NEAR(info.beta4(sct::GlauberSpecies::Cu63), -0.006, 1e-5);
  EXPECT_NEAR(info.hulthenA(sct::GlauberSpecies::Cu63), 0.0, 1e-5);
  EXPECT_NEAR(info.hulthenB(sct::GlauberSpecies::Cu63), 0.0, 1e-5);
  EXPECT_EQ(info.name(sct::GlauberSpecies::Cu63), "Cu63");
}

TEST(NucleusInfo, p1) {
  sct::NucleusInfo& info = sct::NucleusInfo::instance();

  EXPECT_EQ(info.massNumber(sct::GlauberSpecies::p1), 1);
  EXPECT_EQ(info.chargeNumber(sct::GlauberSpecies::p1), 1);
  EXPECT_NEAR(info.radius(sct::GlauberSpecies::p1), 1e-5, 1e-5);
  EXPECT_NEAR(info.skinDepth(sct::GlauberSpecies::p1), 0, 1e-5);
  EXPECT_NEAR(info.radiusError(sct::GlauberSpecies::p1), 0, 1e-5);
  EXPECT_NEAR(info.skinDepthError(sct::GlauberSpecies::p1), 0, 1e-5);
  EXPECT_NEAR(info.beta2(sct::GlauberSpecies::p1), 0, 1e-5);
  EXPECT_NEAR(info.beta4(sct::GlauberSpecies::p1), 0, 1e-5);
  EXPECT_NEAR(info.hulthenA(sct::GlauberSpecies::p1), 0.0, 1e-5);
  EXPECT_NEAR(info.hulthenB(sct::GlauberSpecies::p1), 0.0, 1e-5);
  EXPECT_EQ(info.name(sct::GlauberSpecies::p1), "p1");
}

TEST(NucleusInfo, d2) {
  sct::NucleusInfo& info = sct::NucleusInfo::instance();

  EXPECT_EQ(info.massNumber(sct::GlauberSpecies::d2), 2);
  EXPECT_EQ(info.chargeNumber(sct::GlauberSpecies::d2), 1);
  EXPECT_NEAR(info.radius(sct::GlauberSpecies::d2), 0, 1e-5);
  EXPECT_NEAR(info.skinDepth(sct::GlauberSpecies::d2), 0, 1e-5);
  EXPECT_NEAR(info.radiusError(sct::GlauberSpecies::d2), 0, 1e-5);
  EXPECT_NEAR(info.skinDepthError(sct::GlauberSpecies::d2), 0, 1e-5);
  EXPECT_NEAR(info.beta2(sct::GlauberSpecies::d2), 0, 1e-5);
  EXPECT_NEAR(info.beta4(sct::GlauberSpecies::d2), 0, 1e-5);
  EXPECT_NEAR(info.hulthenA(sct::GlauberSpecies::d2), 0.228, 1e-5);
  EXPECT_NEAR(info.hulthenB(sct::GlauberSpecies::d2), 1.177, 1e-5);
  EXPECT_EQ(info.name(sct::GlauberSpecies::d2), "d2");
}