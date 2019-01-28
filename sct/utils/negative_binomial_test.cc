#include "sct/utils/negative_binomial.h"
#include "gtest/gtest.h"
#include "sct/lib/logging.h"

#include <memory>
#include <random>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"

TEST(nbd, test_nbd_dist) {
  TF1* root_nbd = new TF1(
      "fnbd", "[0]*ROOT::Math::negative_binomial_pdf(x,[1],[2])", 0, 50);
  sct::NegativeBinomial nbd;
  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_real_distribution<> dis(0.0, 5.0);

  for (int i = 0; i < 100; ++i) {
    double npp = dis(generator);
    double k = dis(generator);
    nbd.setParameters(npp, k);
    TH1D* tmp = nbd.distribution();
    root_nbd->SetParameters(1.0, k / (npp + k), k);
    root_nbd->FixParameter(0, 1.0);
    root_nbd->FixParameter(2, k);

    for (int j = 0; j < 50; ++j) {
      EXPECT_NEAR(nbd.evaluateNBD(j), root_nbd->Eval(j), 1e-3);
    }
    EXPECT_LE(tmp->Chisquare(root_nbd), 1e-3);
  }
}

TEST(nbd, test_mean) {
  sct::NegativeBinomial nbd;
  TH1D* tmp = new TH1D("h", "h", 100, 0, 100);
  for (int i = 0; i < 1e6; ++i) {
    tmp->Fill((int)nbd.random());
  }
  EXPECT_NEAR(nbd.npp(), tmp->GetMean(), 1e-3);
}
