#include "sct/glauber/mc_glauber.h"
#include "sct/glauber/glauber_tree.h"
#include "sct/lib/logging.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/nucleus_info.h"

#include "gtest/gtest.h"

#include "TH2.h"
#include "TProfile.h"

TEST(MCGlauber, eventCounter) {
  int nEvents = 5;

  sct::MCGlauber generator;
  generator.run(nEvents);
  sct::GlauberTree* tree = generator.results();

  EXPECT_EQ(tree->getEntries(), nEvents);
}

TEST(MCGlauber, defaultHeader) {
  int nEvents = 1;

  sct::MCGlauber generator;
  generator.run(nEvents);
  sct::GlauberTree* tree = generator.results();

  EXPECT_EQ(tree->getEntries(), nEvents);
  tree->getHeaderEntry(0);
  EXPECT_EQ(tree->nameNucleusA(), "Au197");
  EXPECT_EQ(tree->nameNucleusB(), "Au197");
  EXPECT_EQ(tree->massNumberA(), 197);
  EXPECT_EQ(tree->massNumberB(), 197);
  EXPECT_NEAR(tree->radiusA(),
              sct::NucleusInfo::instance().radius(sct::GlauberSpecies::Au197),
              1e-8);
  EXPECT_NEAR(tree->radiusB(),
              sct::NucleusInfo::instance().radius(sct::GlauberSpecies::Au197),
              1e-8);
  EXPECT_NEAR(
      tree->skinDepthA(),
      sct::NucleusInfo::instance().skinDepth(sct::GlauberSpecies::Au197), 1e-8);
  EXPECT_NEAR(
      tree->skinDepthB(),
      sct::NucleusInfo::instance().skinDepth(sct::GlauberSpecies::Au197), 1e-8);
  EXPECT_NEAR(tree->beta2A(), 0.0, 1e-8);
  EXPECT_NEAR(tree->beta2B(), 0.0, 1e-8);
  EXPECT_NEAR(tree->beta4A(), 0.0, 1e-8);
  EXPECT_NEAR(tree->beta4B(), 0.0, 1e-8);
  EXPECT_NEAR(tree->sigmaNN(), 4.2, 1e-8);
  EXPECT_NEAR(tree->sqrtSNN(), 200.0, 1e-8);
  EXPECT_EQ(tree->smearHardCore(), false);
  EXPECT_EQ(tree->smearGaussian(), false);
  EXPECT_EQ(tree->collisionHardCore(), true);
  EXPECT_EQ(tree->collisionGaussian(), false);
  EXPECT_NEAR(tree->BMax(), 20, 1e-8);
  EXPECT_NEAR(tree->BMin(), 0.0, 1e-8);
}

TEST(MCGlauber, nonDefaultHeader) {
  int nEvents = 1;

  unsigned massNumberA = 195;
  unsigned massNumberB = 200;
  double radiusA = 6.0;
  double radiusB = 6.0;
  double skinDepthA = 0.4;
  double skinDepthB = 1.0;
  double beta2A = 0.01;
  double beta2B = 0.02;
  double beta4A = 0.03;
  double beta4B = 0.04;
  double inelasticXsec = 4.2;
  double energy = 100;
  double bMin = 0.5;
  double bMax = 10.0;

  sct::parameter_list paramsA;
  paramsA["radius"] = radiusA;
  paramsA["skin_depth"] = skinDepthA;
  paramsA["beta2"] = beta2A;
  paramsA["beta4"] = beta4A;
  sct::parameter_list paramsB;
  paramsB["radius"] = radiusB;
  paramsB["skin_depth"] = skinDepthB;
  paramsB["beta2"] = beta2B;
  paramsB["beta4"] = beta4B;

  sct::MCGlauber generator(massNumberA, sct::NucleonPDF::PDF::WoodsSaxon2D, paramsA,
                           massNumberB, sct::NucleonPDF::PDF::WoodsSaxon2D, paramsB,
                           inelasticXsec, energy);
  generator.setSmearing(sct::NucleonSmearing::HardCore);
  generator.setCollisionProfile(sct::CollisionProfile::Gaussian);
  generator.setImpactParameterRange(bMin, bMax);
  generator.run(nEvents);
  sct::GlauberTree* tree = generator.results();

  EXPECT_EQ(tree->getEntries(), nEvents);
  tree->getHeaderEntry(0);
  EXPECT_EQ(tree->nameNucleusA(), sct::MakeString(massNumberA));
  EXPECT_EQ(tree->nameNucleusB(), sct::MakeString(massNumberB));
  EXPECT_EQ(tree->massNumberA(), massNumberA);
  EXPECT_EQ(tree->massNumberB(), massNumberB);
  EXPECT_NEAR(tree->radiusA(), radiusA, 1e-8);
  EXPECT_NEAR(tree->radiusB(), radiusB, 1e-8);
  EXPECT_NEAR(tree->skinDepthA(), skinDepthA, 1e-8);
  EXPECT_NEAR(tree->skinDepthB(), skinDepthB, 1e-8);
  EXPECT_NEAR(tree->beta2A(), beta2A, 1e-8);
  EXPECT_NEAR(tree->beta2B(), beta2B, 1e-8);
  EXPECT_NEAR(tree->beta4A(), beta4A, 1e-8);
  EXPECT_NEAR(tree->beta4B(), beta4B, 1e-8);
  EXPECT_NEAR(tree->sigmaNN(), inelasticXsec, 1e-8);
  EXPECT_NEAR(tree->sqrtSNN(), energy, 1e-8);
  EXPECT_EQ(tree->smearHardCore(), true);
  EXPECT_EQ(tree->smearGaussian(), false);
  EXPECT_EQ(tree->collisionHardCore(), false);
  EXPECT_EQ(tree->collisionGaussian(), true);
  EXPECT_NEAR(tree->BMax(), bMax, 1e-8);
  EXPECT_NEAR(tree->BMin(), bMin, 1e-8);
}

TEST(MCGlauber, nonDefaultSymmetricHeader) {
  int nEvents = 1;

  unsigned massNumber = 195;
  double radius = 6.0;
  double skinDepth = 0.4;
  double beta2 = 0.01;
  double beta4 = 0.03;
  double inelasticXsec = 4.2;
  double energy = 100;
  double bMin = 0.5;
  double bMax = 10.0;

  sct::parameter_list params;
  params["radius"] = radius;
  params["skin_depth"] = skinDepth;
  params["beta2"] = beta2;
  params["beta4"] = beta4;

  sct::MCGlauber generator(massNumber, sct::NucleonPDF::PDF::WoodsSaxon2D, params,
                           inelasticXsec, energy);
  generator.setSmearing(sct::NucleonSmearing::Gaussian);
  generator.setCollisionProfile(sct::CollisionProfile::HardCore);
  generator.setImpactParameterRange(bMin, bMax);
  generator.run(nEvents);
  sct::GlauberTree* tree = generator.results();

  EXPECT_EQ(tree->getEntries(), nEvents);
  tree->getHeaderEntry(0);
  EXPECT_EQ(tree->nameNucleusA(), sct::MakeString(massNumber));
  EXPECT_EQ(tree->nameNucleusB(), sct::MakeString(massNumber));
  EXPECT_EQ(tree->massNumberA(), massNumber);
  EXPECT_EQ(tree->massNumberB(), massNumber);
  EXPECT_NEAR(tree->radiusA(), radius, 1e-8);
  EXPECT_NEAR(tree->radiusB(), radius, 1e-8);
  EXPECT_NEAR(tree->skinDepthA(), skinDepth, 1e-8);
  EXPECT_NEAR(tree->skinDepthB(), skinDepth, 1e-8);
  EXPECT_NEAR(tree->beta2A(), beta2, 1e-8);
  EXPECT_NEAR(tree->beta2B(), beta2, 1e-8);
  EXPECT_NEAR(tree->beta4A(), beta4, 1e-8);
  EXPECT_NEAR(tree->beta4B(), beta4, 1e-8);
  EXPECT_NEAR(tree->sigmaNN(), inelasticXsec, 1e-8);
  EXPECT_NEAR(tree->sqrtSNN(), energy, 1e-8);
  EXPECT_EQ(tree->smearHardCore(), false);
  EXPECT_EQ(tree->smearGaussian(), true);
  EXPECT_EQ(tree->collisionHardCore(), true);
  EXPECT_EQ(tree->collisionGaussian(), false);
  EXPECT_NEAR(tree->BMax(), bMax, 1e-8);
  EXPECT_NEAR(tree->BMin(), bMin, 1e-8);
}
