#include "gtest/gtest.h"
#include "sct/centrality/nbd_fit.h"

#include <random>

#include "sct/core/base.h"
#include "sct/core/logging.h"
#include "sct/core/enumerations.h"

#include "TH1D.h"

TEST(NBDFit, norm) {
  
  const double norm_factor = 2.0;
  
  std::unique_ptr<TH1D> h1 = sct::make_unique<TH1D>("h1", "", 150, 0, 150);
  std::unique_ptr<TH1D> h2 = sct::make_unique<TH1D>("h2", "", 150, 0, 150);
  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_int_distribution<> dis(0, 150);
  for (int i = 0; i <= 1e5; ++i) {
    h1->Fill(dis(generator));
    h2->Fill(dis(generator), norm_factor);
  }
  
  sct::NBDFit fitter;
  double calculated_norm = fitter.norm(h1.get(), h2.get());
  
  EXPECT_NEAR(calculated_norm, 1.0/norm_factor, 1e-2);
}

TEST(NBDFit, chi2) {
  
  std::unique_ptr<TH1D> h1 = sct::make_unique<TH1D>("h1", "", 150, 0, 150);
  h1->Sumw2();
  std::unique_ptr<TH1D> h2 = sct::make_unique<TH1D>("h2", "", 150, 0, 150);
  h2->Sumw2();
  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_int_distribution<> dis(0, 150);
  for (int i = 0; i <= 1e6; ++i) {
    h1->Fill(dis(generator));
    h2->Fill(dis(generator));
  }
  
  double chi2;
  int ndf, good;
  sct::string opts = "UU NORM";
  h1->Chi2TestX(h2.get(), chi2, ndf, good, opts.c_str());
  double root_result = chi2 / ndf;
  
  sct::NBDFit fitter;
  fitter.minimumMultiplicityCut(0);
  auto fitter_result = fitter.chi2(h1.get(), h2.get());
  double fitter_chi2 = fitter_result.first / fitter_result.second;
  
  EXPECT_NEAR(root_result, fitter_chi2, 1e-2);
}

TEST(NBDFit, chi2_fail) {
  
  std::unique_ptr<TH1D> h1 = sct::make_unique<TH1D>("h1", "", 150, 0, 150);
  h1->Sumw2();
  std::unique_ptr<TH1D> h2 = sct::make_unique<TH1D>("h2", "", 150, 0, 150);
  h2->Sumw2();
  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_int_distribution<> dis(0, 150);
  std::piecewise_linear_distribution<> unit_linear(1, 0.0, 1.0, [](double x) {return x;});
  for (int i = 0; i <= 1e6; ++i) {
    h1->Fill(dis(generator));
    h2->Fill(unit_linear(generator));
  }
  
  double chi2;
  int ndf, good;
  sct::string opts = "UU NORM";
  h1->Chi2TestX(h2.get(), chi2, ndf, good, opts.c_str());
  double root_result = chi2 / ndf;
  
  sct::NBDFit fitter;
  fitter.minimumMultiplicityCut(0);
  auto fitter_result = fitter.chi2(h1.get(), h2.get());
  double fitter_chi2 = fitter_result.first / fitter_result.second;
  
  EXPECT_NEAR(root_result, fitter_chi2, 1e-2);
  EXPECT_GE(fitter_chi2, 100);
}

TEST(NBDFit, chi2Weighted) {
  
  std::unique_ptr<TH1D> h1 = sct::make_unique<TH1D>("h1", "", 150, 0, 150);
  std::unique_ptr<TH1D> h2 = sct::make_unique<TH1D>("h2", "", 150, 0, 150);
  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_int_distribution<> dis(0, 150);
  for (int i = 0; i <= 1e6; ++i) {
    h1->Fill(dis(generator));
    h2->Fill(dis(generator));
    h2->Fill(dis(generator));
  }
  
  h2->Scale(1.0/2.0);
  
  double chi2;
  int ndf, good;
  sct::string opts = "UU NORM";
  h1->Chi2TestX(h2.get(), chi2, ndf, good, opts.c_str());
  double root_result = chi2 / ndf;
  
  sct::NBDFit fitter;
  fitter.minimumMultiplicityCut(0);
  auto fitter_result = fitter.chi2(h1.get(), h2.get());
  double fitter_chi2 = fitter_result.first / fitter_result.second;
  
  EXPECT_NEAR(root_result, fitter_chi2, 1e-2);
}

TEST(NBDFit, chi2RestrictedRange) {
  
  std::unique_ptr<TH1D> h1 = sct::make_unique<TH1D>("h1", "", 150, 0, 150);
  h1->Sumw2();
  std::unique_ptr<TH1D> h2 = sct::make_unique<TH1D>("h2", "", 150, 0, 150);
  h2->Sumw2();
  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_int_distribution<> dis(0, 150);
  for (int i = 0; i <= 1e6; ++i) {
    h1->Fill(dis(generator));
    h2->Fill(dis(generator));
  }
  
  double chi2;
  int ndf, good;
  sct::string opts = "UU NORM";
  h1->Chi2TestX(h2.get(), chi2, ndf, good, opts.c_str());
  double root_result = chi2 / ndf;
  
  sct::NBDFit fitter;
  fitter.minimumMultiplicityCut(50);
  auto fitter_result = fitter.chi2(h1.get(), h2.get());
  double fitter_chi2 = fitter_result.first / fitter_result.second;
  
  EXPECT_GE(chi2, fitter_result.first);
  EXPECT_GE(ndf, fitter_result.second);
  EXPECT_EQ(fitter_result.second, 99);
}
