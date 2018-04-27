// sct/core/nucleus_test.cc

// *notice* some of these tests are statistical in nature,
// so failure does not necessarily imply a bug

#include "gtest/gtest.h"
#include "sct/core/nucleus.hh"

#include "sct/core/base.hh"
#include "sct/core/logging.hh"
#include "sct/core/enumerations.hh"
#include "sct/core/nucleon.hh"

#include "TH3.h"
#include "TH2.h"
#include "TH1.h"


TEST(nucleus, defaultGold) {

  sct::Nucleus nucleus;
  nucleus.setParameters(197, 6.38, 0.535, 0.0, 0.0);
  nucleus.setRepulsionDistance(0.01);

  nucleus.generate();

  EXPECT_EQ(nucleus.size(), 197);
  EXPECT_EQ(nucleus.massNumber(), 197);
  EXPECT_EQ(nucleus.beta2(), 0.0);
  EXPECT_EQ(nucleus.beta4(), 0.0);
  EXPECT_EQ(nucleus.radius(), 6.38);
  EXPECT_EQ(nucleus.skinDepth(), 0.535);
  EXPECT_EQ(nucleus.nnCrossSection(), 0.0);
  EXPECT_EQ(nucleus.repulsionDistance(), 0.01);
  EXPECT_LE(abs(nucleus.nucleusTheta()), sct::pi); // [0, pi]
  EXPECT_LE(abs(nucleus.nucleusPhi()), sct::pi); // [-pi, pi]

}

TEST(nucleus, repulsionDistance) {

  sct::Nucleus nucleus;
  nucleus.setParameters(197, 6.38, 0.535, 0.0, 0.0);
  nucleus.setNucleonSmearing(sct::NucleonSmearing::None, 0.0);
  nucleus.setRepulsionDistance(0.05);

  for (int event = 0; event < 1e4; ++event) {
    nucleus.generate();
    for (int i = 0; i < nucleus.size(); ++i) {
      for (int j = i + 1; j < nucleus.size(); ++j) {
        EXPECT_GE(nucleus[i].deltaR(nucleus[j]), 0.05);
      }
    }
  }
}

TEST(nucleus, smearingHardCore) {
  // this will be a statistical test of the width of the generated nucleon
  // distribution in x, y, z
  sct::Nucleus nucleus;
  nucleus.setParameters(197, 6.38, 0.535, 0.0, 0.0);
  nucleus.setRepulsionDistance(0.00);

  sct::Nucleus smearedNucleus;
  smearedNucleus.setParameters(197, 6.38, 0.535, 0.0, 0.0);
  smearedNucleus.setNucleonSmearing(sct::NucleonSmearing::HardCore, 30.0);
  smearedNucleus.setRepulsionDistance(0.00);

  for (int event = 0; event < 1e4; ++event) {
    nucleus.generate();
    smearedNucleus.generate();
  }

  TH3D* position = nucleus.smearedPosition();
  TH3D* smearedPosition = smearedNucleus.smearedPosition();

  TH1D* ProjectionX_1 = position->ProjectionX();
  TH1D* ProjectionX_2 = smearedPosition->ProjectionX();

  TH1D* ProjectionY_1 = position->ProjectionY();
  TH1D* ProjectionY_2 = smearedPosition->ProjectionY();

  TH1D* ProjectionZ_1 = position->ProjectionZ();
  TH1D* ProjectionZ_2 = smearedPosition->ProjectionZ();

  EXPECT_LE(ProjectionX_1->KolmogorovTest(ProjectionX_2), 0.01);
  EXPECT_LE(ProjectionY_1->KolmogorovTest(ProjectionY_2), 0.01);
  EXPECT_LE(ProjectionZ_1->KolmogorovTest(ProjectionZ_2), 0.01);
}

TEST(nucleus, smearingGaussian) {
  // this will be a statistical test of the width of the generated nucleon
  // distribution in x, y, z
  sct::Nucleus nucleus;
  nucleus.setParameters(197, 6.38, 0.535, 0.0, 0.0);
  nucleus.setRepulsionDistance(0.00);

  sct::Nucleus smearedNucleus;
  smearedNucleus.setParameters(197, 6.38, 0.535, 0.0, 0.0);
  smearedNucleus.setNucleonSmearing(sct::NucleonSmearing::Gaussian, 5.0);
  smearedNucleus.setRepulsionDistance(0.00);

  for (int event = 0; event < 1e4; ++event) {
    nucleus.generate();
    smearedNucleus.generate();
  }

  TH3D* position = nucleus.smearedPosition();
  TH3D* smearedPosition = smearedNucleus.smearedPosition();

  TH1D* ProjectionX_1 = position->ProjectionX();
  TH1D* ProjectionX_2 = smearedPosition->ProjectionX();

  TH1D* ProjectionY_1 = position->ProjectionY();
  TH1D* ProjectionY_2 = smearedPosition->ProjectionY();

  TH1D* ProjectionZ_1 = position->ProjectionZ();
  TH1D* ProjectionZ_2 = smearedPosition->ProjectionZ();

  EXPECT_LE(ProjectionX_1->KolmogorovTest(ProjectionX_2), 0.01);
  EXPECT_LE(ProjectionY_1->KolmogorovTest(ProjectionY_2), 0.01);
  EXPECT_LE(ProjectionZ_1->KolmogorovTest(ProjectionZ_2), 0.01);
}

TEST(nucleus, sphericalDistribution) {
  // this will be a statistical test of the width of the generated nucleon
  // distribution in x, y, z
  sct::Nucleus nucleus;
  nucleus.setParameters(197, 6.38, 0.535, 0.0, 0.0);
  nucleus.setRepulsionDistance(0.00);

  for (int event = 0; event < 1e4; ++event) {
    nucleus.generate();
  }

  TH3D* position = nucleus.smearedPosition();

  TH1D* ProjectionX = position->ProjectionX();
  TH1D* ProjectionY = position->ProjectionY();
  TH1D* ProjectionZ = position->ProjectionZ();

  EXPECT_GE(ProjectionX->KolmogorovTest(ProjectionY), 0.01);
  EXPECT_GE(ProjectionY->KolmogorovTest(ProjectionZ), 0.01);
  EXPECT_GE(ProjectionZ->KolmogorovTest(ProjectionX), 0.01);
}

TEST(nucleus, woodsSaxon1D) {

  sct::Nucleus nucleus;
  double radius = 6.38;
  double skinDepth = 0.535;
  nucleus.setParameters(197, radius, skinDepth, 0.0, 0.0);
  nucleus.setRepulsionDistance(0.00);

  // set the fitting function for the woods-saxon
  TF1* WS = new TF1("WS", "[2]*x*x/(1.0 + TMath::Exp((x - [0]) / [1]))", 0, 20);
  WS->SetParameter(0, 6.38);
  WS->SetParameter(1, 0.535);
  WS->SetParameter(2, 1);

  TH1D* hWS = new TH1D("h", "h", 100, 0, 20);

  for (int i = 0; i < 5 * 1e4; ++i) {
    nucleus.generate();
    for (int j = 0; j < nucleus.size(); ++j) {
      hWS->Fill(nucleus[j].r());
    }
  }
  hWS->Fit(WS, "Q");

  EXPECT_NEAR(radius, WS->GetParameter(0), 5e-3);
  EXPECT_NEAR(skinDepth, WS->GetParameter(1), 5e-3);

}

double WoodsSaxonDeformed(double *x,
                          double *par) {
  double r = x[0];
  double cosTheta = x[1];
  double R0 = par[0];
  double a = par[1];
  double beta2 = par[2];
  double beta4 = par[3];
  double norm = par[4];

  // powers of cosTheta
  double cosTheta2 = std::pow(cosTheta, 2.0);
  double cosTheta4 = std::pow(cosTheta, 4.0);

  // spherical harmonics Y(l=2, m=0), Y(l=4, m=0)
  double Y20 = std::sqrt(5.0 / sct::pi) / 4.0 * (3.0 * cosTheta2 - 1.0);
  double Y40 = std::sqrt(1.0 / sct::pi) * 3.0 / 16.0 * (35.0 * cosTheta4 - 30.0 * cosTheta2 + 3.0);

  double R = R0 * (1.0 + beta2 * Y20 + beta4 * Y40);

  return norm * std::pow(r, 2) / (1.0 + std::exp((r - R) / a));
}

TEST(nucleus, woodsSaxon2D) {
  
  sct::Nucleus nucleus;
  double radius = 6.38;
  double skinDepth = 0.535;
  double beta2 = 0.2;
  double beta4 = 0;
  nucleus.setParameters(197, radius, skinDepth, beta2, beta4);
  nucleus.setRepulsionDistance(0.00);
  nucleus.setRandomOrientation(false);
  
  TF2* fWSD = new TF2("WSD", WoodsSaxonDeformed, 0, 20, -1, 1, 5);
  fWSD->SetNpx(200);
  fWSD->SetNpy(200);
  fWSD->SetParameter(0, radius + 0.01);
  fWSD->SetParameter(1, skinDepth + 0.01);
  fWSD->SetParameter(2, beta2 + 0.01);
  fWSD->SetParameter(3, beta4 + 0.01);
  fWSD->SetParameter(4, 1);
  
  TH2D* rcostheta = new TH2D("rcostheta", ";R;cos(#theta)", 100, 0, 20, 100, -1, 1);
  
  for (int i = 0; i < 5 * 1e4; ++i) {
    nucleus.generate();
    for (int j = 0; j < nucleus.size(); ++j) {
      rcostheta->Fill(nucleus[j].r(), cos(nucleus[j].theta()));
    }
  }
  
  rcostheta->Fit(fWSD, "Q");
  
  EXPECT_NEAR(fWSD->GetParameter(0), radius, 1e-2);
  EXPECT_NEAR(fWSD->GetParameter(1), skinDepth, 1e-2);
  EXPECT_NEAR(fWSD->GetParameter(2), beta2, 1e-2);
  EXPECT_NEAR(fWSD->GetParameter(3), beta4, 1e-2);
}
