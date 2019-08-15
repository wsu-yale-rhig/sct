#include "gtest/gtest.h"

#include "sct/lib/logging.h"
#include "sct/utils/negative_binomial.h"

#include <memory>
#include <random>
#include <vector>

#include "stglauber/StCentralityMaker/StNegativeBinomial.h"

#include "TH1D.h"

TEST(nbd, test_star_distribution_comp_default) {
  sct::NegativeBinomial nbd;
  StNegativeBinomial star_nbd;

  TH1D *my_nbd_hist = nbd.distribution();
  TH1D *star_nbd_hist = star_nbd.distribution();

  EXPECT_EQ(my_nbd_hist->GetNbinsX(), star_nbd_hist->GetNbinsX());
  EXPECT_NEAR(my_nbd_hist->GetMean(), star_nbd_hist->GetMean(), 1e-5);
  EXPECT_NEAR(my_nbd_hist->GetStdDev(), star_nbd_hist->GetStdDev(), 1e-5);
}

TEST(nbd, test_star_distribution_comp_non_default) {
  sct::NegativeBinomial nbd(2.0, 3.0);
  StNegativeBinomial star_nbd(2.0, 3.0);

  TH1D *my_nbd_hist = nbd.distribution();
  TH1D *star_nbd_hist = star_nbd.distribution();

  EXPECT_EQ(my_nbd_hist->GetNbinsX(), star_nbd_hist->GetNbinsX());
  EXPECT_NEAR(my_nbd_hist->GetMean(), star_nbd_hist->GetMean(), 1e-5);
  EXPECT_NEAR(my_nbd_hist->GetStdDev(), star_nbd_hist->GetStdDev(), 1e-5);
  for (int i = 1; i < my_nbd_hist->GetNbinsX(); ++i) {
    EXPECT_NEAR(my_nbd_hist->GetBinContent(i), star_nbd_hist->GetBinContent(i),
                1e-5);
  }
}

TEST(nbd, test_star_distribution_comp_sample) {
  sct::NegativeBinomial nbd(2.0, 3.0);
  StNegativeBinomial star_nbd(2.0, 3.0);

  TH1D *sample_my = new TH1D("my_sample", "", 100, 0, 100);
  sample_my->Sumw2();
  TH1D *sample_star = new TH1D("star_sample", "", 100, 0, 100);
  sample_star->Sumw2();
  for (int i = 0; i < 1e8; ++i) {
    sample_my->Fill(nbd.distribution()->GetRandom());
    sample_star->Fill(star_nbd.distribution()->GetRandom());
  }

  EXPECT_GE(sample_my->Chi2Test(sample_star), 0.2);
  EXPECT_NEAR(sample_my->GetMean(), sample_star->GetMean(), 1e-3);
  EXPECT_NEAR(sample_my->GetStdDev(), sample_star->GetStdDev(), 1e-3);
}
