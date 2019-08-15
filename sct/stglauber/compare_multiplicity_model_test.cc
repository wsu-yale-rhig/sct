#include "gtest/gtest.h"

#include "sct/centrality/multiplicity_model.h"
#include "sct/lib/logging.h"
#include "sct/utils/negative_binomial.h"
#include "sct/utils/random.h"

#include <memory>
#include <random>
#include <vector>

#include "stglauber/StCentralityMaker/StNegativeBinomial.h"

#include "TH1D.h"
#include "TMath.h"
#include "TRandom.h"

TEST(mult_model, test_star_multiplicity_comp_tcm) {
  sct::MultiplicityModel model(2.38, 2.00, 0.13, 0.98, 0.84, 540, 1.0, false);
  StNegativeBinomial star_nbd(2.38, 2.0, 0.13, 0.14, 1.0, false);

  std::vector<int> npart;
  std::vector<int> ncoll;

  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_int_distribution<> dis(0, 200);

  for (int i = 0; i < 100; ++i) {
    npart.push_back(dis(generator));
    ncoll.push_back(dis(generator));
  }

  for (auto &np : npart) {
    for (auto &nc : ncoll) {
      EXPECT_EQ(model.twoComponentMultiplicity(np, nc),
                star_nbd.GetTwoComponentMultiplicity(np, nc));
    }
  }
}

TEST(mult_model, test_star_multiplicity_comp_mult) {
  sct::MultiplicityModel model(2.38, 2.00, 0.13, 0.98, 0.84, 540, 1.0, false);
  StNegativeBinomial star_nbd(2.38, 2.0, 0.13, 0.14, 1.0, false);

  int npart = 100;
  int ncoll = 300;
  int nevents = 1e6;

  TH1D *h1 = new TH1D("h1", "", 300, 50, 350);
  TH1D *h2 = new TH1D("h2", "", 300, 50, 350);

  for (int i = 0; i < nevents; ++i) {
    double model_mult = model.multiplicity(npart, ncoll);
    double star_mult = star_nbd.GetMultiplicity(npart, ncoll);

    h1->Fill(model_mult);
    h2->Fill(star_mult);
  }

  double test_stat = h1->AndersonDarlingTest(h2);
  EXPECT_GE(test_stat, 0.1);
  EXPECT_NEAR(h1->GetMean(), h2->GetMean(),
              h1->GetMeanError() + h2->GetMeanError());
  EXPECT_NEAR(h1->GetStdDev(), h2->GetStdDev(),
              h1->GetStdDevError() + h2->GetStdDevError());
}