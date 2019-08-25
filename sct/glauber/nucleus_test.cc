#include "sct/glauber/nucleon.h"
#include "sct/glauber/nucleus.h"
#include "sct/lib/enumerations.h"
#include "sct/lib/logging.h"
#include "sct/lib/math.h"
#include "sct/utils/nucleus_info.h"

#include <random>

#include "gtest/gtest.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

// *notice* some of these tests are statistical in nature,
// so failure does not necessarily imply a bug

// first test the nucleon rotation methods
// next, position when we also have a nuclear rotation
TEST(nucleon, nucleusRotation) {
  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_real_distribution<> dis(-2.0, 2.0);
  sct::Nucleus nucleus;

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

    // and zero impact parameter
    double b = 0.0;

    nucleon.set(pos);

    nucleus.rotateAndOffset(nucleon);

    double nuclear_theta = nucleus.nucleusTheta();
    double nuclear_phi = nucleus.nucleusPhi();

    // now create a vector to replicate the rotation
    TVector3 vec(x, y, z);
    vec.RotateY(nuclear_theta);
    vec.RotateZ(nuclear_phi);

    // check r, theta, phi, x, y, z
    EXPECT_NEAR(vec.X(), nucleon.x(), 1e-8);
    EXPECT_NEAR(vec.Y(), nucleon.y(), 1e-8);
    EXPECT_NEAR(vec.Z(), nucleon.z(), 1e-8);
    EXPECT_NEAR(vec.Mag(), nucleon.r(), 1e-8);
    EXPECT_NEAR(vec.Theta(), nucleon.theta(), 1e-8);
    EXPECT_NEAR(vec.Phi(), nucleon.phi(), 1e-8);
  }
}

TEST(nucleus, defaultGold) {
  sct::Nucleus nucleus;
  sct::GlauberSpecies species = sct::GlauberSpecies::Au197;
  nucleus.setParameters(species);
  nucleus.setRepulsionDistance(0.01);

  nucleus.generate();

  sct::parameter_list params = nucleus.nuclearPDF().parameters();

  EXPECT_EQ(nucleus.size(), sct::NucleusInfo::instance().massNumber(species));
  EXPECT_EQ(nucleus.massNumber(),
            sct::NucleusInfo::instance().massNumber(species));
  EXPECT_EQ(params["beta2"], 0.0);
  EXPECT_EQ(params["beta4"], 0.0);
  EXPECT_EQ(params["radius"], sct::NucleusInfo::instance().radius(species));
  EXPECT_EQ(params["skin_depth"],
            sct::NucleusInfo::instance().skinDepth(species));
  EXPECT_EQ(nucleus.repulsionDistance(), 0.01);
  EXPECT_LE(abs(nucleus.nucleusTheta()), sct::pi); // [0, pi]
  EXPECT_LE(abs(nucleus.nucleusPhi()), sct::pi);   // [-pi, pi]
}

TEST(nucleus, repulsionDistance) {
  sct::Nucleus nucleus;
  nucleus.setParameters(sct::GlauberSpecies::Au197);
  nucleus.setNucleonSmearing(sct::NucleonSmearing::None, 0.0);
  double repulsion = 0.05;
  nucleus.setRepulsionDistance(repulsion);

  for (int event = 0; event < 1e4; ++event) {
    nucleus.generate();
    for (int i = 0; i < nucleus.size(); ++i) {
      for (int j = i + 1; j < nucleus.size(); ++j) {
        EXPECT_GE(nucleus[i].deltaR(nucleus[j]), repulsion);
      }
    }
  }
}

TEST(nucleus, sphericalDistribution) {
  // this will be a statistical test of the width of the generated nucleon
  // distribution in x, y, z
  sct::Nucleus nucleus;
  nucleus.setParameters(sct::GlauberSpecies::Au197);
  nucleus.setRepulsionDistance(0.00);

  for (int event = 0; event < 1e4; ++event) {
    nucleus.generate();
  }

  TH3D *position = nucleus.smearedPosition();

  TH1D *projection_x = position->ProjectionX();
  TH1D *projection_y = position->ProjectionY();
  TH1D *projection_z = position->ProjectionZ();

  EXPECT_GE(projection_x->KolmogorovTest(projection_y), 0.01);
  EXPECT_GE(projection_y->KolmogorovTest(projection_z), 0.01);
  EXPECT_GE(projection_z->KolmogorovTest(projection_x), 0.01);
}

TEST(nucleus, woodsSaxon1D) {
  sct::Nucleus nucleus;
  sct::GlauberSpecies species = sct::GlauberSpecies::Au197;
  nucleus.setParameters(species);
  nucleus.setRepulsionDistance(0.00);

  // set the fitting function for the woods-saxon
  TF1 *WS = new TF1("WS", "[2]*x*x/(1.0 + TMath::Exp((x - [0]) / [1]))", 0, 20);
  WS->SetParameter(0, sct::NucleusInfo::instance().radius(species));
  WS->SetParameter(1, sct::NucleusInfo::instance().skinDepth(species));
  WS->SetParameter(2, 1);

  TH1D *hWS = new TH1D("h", "h", 100, 0, 20);

  for (int i = 0; i < 5 * 1e4; ++i) {
    nucleus.generate();
    for (int j = 0; j < nucleus.size(); ++j) {
      hWS->Fill(nucleus[j].r());
    }
  }
  hWS->Fit(WS, "Q");

  EXPECT_NEAR(sct::NucleusInfo::instance().radius(species), WS->GetParameter(0),
              5e-3);
  EXPECT_NEAR(sct::NucleusInfo::instance().skinDepth(species),
              WS->GetParameter(1), 5e-3);
}

double WoodsSaxonDeformed(double *x, double *par) {
  double r = x[0];
  double costheta = x[1];
  double R0 = par[0];
  double a = par[1];
  double beta2 = par[2];
  double beta4 = par[3];
  double norm = par[4];

  // powers of costheta
  double costheta2 = std::pow(costheta, 2.0);
  double costheta4 = std::pow(costheta, 4.0);

  // spherical harmonics Y(l=2, m=0), Y(l=4, m=0)
  double Y20 = std::sqrt(5.0 / sct::pi) / 4.0 * (3.0 * costheta2 - 1.0);
  double Y40 = std::sqrt(1.0 / sct::pi) * 3.0 / 16.0 *
               (35.0 * costheta4 - 30.0 * costheta2 + 3.0);

  double R = R0 * (1.0 + beta2 * Y20 + beta4 * Y40);

  return norm * std::pow(r, 2) / (1.0 + std::exp((r - R) / a));
}

TEST(nucleus, woodsSaxon2D) {
  sct::Nucleus nucleus;
  sct::GlauberSpecies species = sct::GlauberSpecies::Au197;
  nucleus.setParameters(species, sct::GlauberMod::Nominal, true);
  nucleus.setRepulsionDistance(0.00);
  nucleus.setRandomOrientation(false);

  TF2 *fWSD = new TF2("WSD", WoodsSaxonDeformed, 0, 20, -1, 1, 5);
  fWSD->SetNpx(400);
  fWSD->SetNpy(400);
  fWSD->SetParameter(0, sct::NucleusInfo::instance().radius(species));
  fWSD->SetParameter(1, sct::NucleusInfo::instance().skinDepth(species));
  fWSD->SetParameter(2, sct::NucleusInfo::instance().beta2(species));
  fWSD->SetParameter(3, sct::NucleusInfo::instance().beta4(species));
  fWSD->SetParameter(4, 1.0);

  TH2D *rcostheta =
      new TH2D("rcostheta", ";R;cos(#theta)", 100, 0, 20, 100, -1, 1);

  for (int i = 0; i < 5 * 1e4; ++i) {
    nucleus.generate();
    for (int j = 0; j < nucleus.size(); ++j) {
      rcostheta->Fill(nucleus[j].r(), cos(nucleus[j].theta()));
    }
  }

  rcostheta->Fit(fWSD, "QEM");

  EXPECT_NEAR(fWSD->GetParameter(0),
              sct::NucleusInfo::instance().radius(species), 1e-2);
  EXPECT_NEAR(fWSD->GetParameter(1),
              sct::NucleusInfo::instance().skinDepth(species), 1e-2);
  EXPECT_NEAR(fWSD->GetParameter(2),
              sct::NucleusInfo::instance().beta2(species), 1e-2);
  EXPECT_NEAR(fWSD->GetParameter(3),
              sct::NucleusInfo::instance().beta4(species), 1e-2);
}
