#include "sct/glauber/nucleon_pdf.h"
#include "sct/utils/nucleus_info.h"
#include "sct/lib/math.h"

#include "gtest/gtest.h"

TEST(nucleonpdf, default_auau) {
  sct::NucleonPDF pdf;

  pdf.init(sct::GlauberSpecies::Au197, sct::GlauberMod::Nominal, false);
  auto params = pdf.parameters();

  sct::parameter_list expected{
      {"radius",
       sct::NucleusInfo::instance().radius(sct::GlauberSpecies::Au197)},
      {"skin_depth",
       sct::NucleusInfo::instance().skinDepth(sct::GlauberSpecies::Au197)}};

  EXPECT_FALSE(pdf.deformed());
  for (auto &par : params) {
    EXPECT_NEAR(par.second, expected[par.first], 1e-5);
  }

  for (int i = 0; i < 1e5; ++i) {
    double r = 0.0;
    double theta = 0.0;
    double phi = 0.0;
    pdf.sample(r, theta, phi);
    EXPECT_GE(r, 0.0);
    EXPECT_LE(r, 20.0);
    EXPECT_GE(theta, 0);
    EXPECT_LE(theta, sct::pi);
    EXPECT_GE(phi, -sct::pi);
    EXPECT_LE(phi, sct::pi);
  }
}

TEST(nucleonpdf, default_auau_deformed) {
  sct::NucleonPDF pdf;

  pdf.init(sct::GlauberSpecies::Au197, sct::GlauberMod::Nominal, true);
  auto params = pdf.parameters();

  sct::parameter_list expected{
      {"radius",
       sct::NucleusInfo::instance().radius(sct::GlauberSpecies::Au197)},
      {"skin_depth",
       sct::NucleusInfo::instance().skinDepth(sct::GlauberSpecies::Au197)},
      {"beta2", sct::NucleusInfo::instance().beta2(sct::GlauberSpecies::Au197)},
      {"beta4",
       sct::NucleusInfo::instance().beta4(sct::GlauberSpecies::Au197)}};

  EXPECT_TRUE(pdf.deformed());
  for (auto &par : params) {
    EXPECT_NEAR(par.second, expected[par.first], 1e-5);
  }
}

TEST(nucleonpdf, default_auau_clear_deformed) {
  sct::NucleonPDF pdf;

  pdf.init(sct::GlauberSpecies::Au197, sct::GlauberMod::Nominal, false);
  pdf.init(sct::GlauberSpecies::Au197, sct::GlauberMod::Nominal, true);
  auto params = pdf.parameters();

  sct::parameter_list expected{
      {"radius",
       sct::NucleusInfo::instance().radius(sct::GlauberSpecies::Au197)},
      {"skin_depth",
       sct::NucleusInfo::instance().skinDepth(sct::GlauberSpecies::Au197)},
      {"beta2", sct::NucleusInfo::instance().beta2(sct::GlauberSpecies::Au197)},
      {"beta4",
       sct::NucleusInfo::instance().beta4(sct::GlauberSpecies::Au197)}};

  EXPECT_TRUE(pdf.deformed());
  for (auto &par : params) {
    EXPECT_NEAR(par.second, expected[par.first], 1e-5);
  }
}

TEST(nucleonpdf, default_auau_clear_spherical) {
  sct::NucleonPDF pdf;

  pdf.init(sct::GlauberSpecies::Au197, sct::GlauberMod::Nominal, true);
  pdf.init(sct::GlauberSpecies::Au197, sct::GlauberMod::Nominal, false);
  auto params = pdf.parameters();

  sct::parameter_list expected{
      {"radius",
       sct::NucleusInfo::instance().radius(sct::GlauberSpecies::Au197)},
      {"skin_depth",
       sct::NucleusInfo::instance().skinDepth(sct::GlauberSpecies::Au197)},
      {"beta2", 0.0},
      {"beta4", 0.0}};

  EXPECT_FALSE(pdf.deformed());
  EXPECT_EQ(params.size(), 2);
  for (auto &par : expected) {
    EXPECT_NEAR(par.second, params[par.first], 1e-5);
  }
}

TEST(nucleonpdf, proton_default) {
  sct::NucleonPDF pdf;

  pdf.init(sct::GlauberSpecies::p1);
  auto params = pdf.parameters();

  sct::parameter_list expected{
      {"d",
       sct::NucleusInfo::instance().radius(sct::GlauberSpecies::p1)}};

  EXPECT_FALSE(pdf.deformed());
  for (auto &par : params) {
    EXPECT_NEAR(par.second, expected[par.first], 1e-5);
  }

  for (int i = 0; i < 1e5; ++i) {
    double r = 0.0;
    double theta = 0.0;
    double phi = 0.0;
    pdf.sample(r, theta, phi);
    EXPECT_GE(r, 0.0);
    EXPECT_LE(r, 1e-3);
    EXPECT_GE(theta, 0);
    EXPECT_LE(theta, sct::pi);
    EXPECT_GE(phi, -sct::pi);
    EXPECT_LE(phi, sct::pi);
  }
}

TEST(nucleonpdf, deuteron_default) {
  sct::NucleonPDF pdf;

  pdf.init(sct::GlauberSpecies::d2);
  auto params = pdf.parameters();

  sct::parameter_list expected{
      {"a",
       sct::NucleusInfo::instance().hulthenA(sct::GlauberSpecies::d2)},
       {"b",
       sct::NucleusInfo::instance().hulthenB(sct::GlauberSpecies::d2)}};

  EXPECT_FALSE(pdf.deformed());
  for (auto &par : params) {
    EXPECT_NEAR(par.second, expected[par.first], 1e-5);
  }
}