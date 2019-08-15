#include "gtest/gtest.h"

#include "sct/centrality/multiplicity_model.h"
#include "sct/centrality/nbd_fit.h"
#include "sct/glauber/glauber_tree.h"
#include "sct/lib/logging.h"
#include "sct/utils/negative_binomial.h"
#include "sct/utils/random.h"

#include <memory>
#include <random>
#include <vector>

#include "boost/filesystem.hpp"

#include "stglauber/StCentralityMaker/StNbdFitMaker.h"
#include "stglauber/StCentralityMaker/StNegativeBinomial.h"
#include "stglauber/StGlauberTree/StGlauberTree.h"

#include "TH1D.h"
#include "TH2D.h"

TEST(compare_fitter, compare_fitter_diff_glauber) {

  std::string refmult_file = "testfiles/refmult.root";
  std::string sct_file = "testfiles/sct_tree.root";
  std::string stglauber_file = "testfiles/stglauber_npartncoll.root";

  if (!boost::filesystem::exists(refmult_file)) {
    LOG(ERROR) << "input file does not exist: " << refmult_file;
    return;
  }
  if (!boost::filesystem::exists(sct_file)) {
    LOG(ERROR) << "input file does not exist: " << sct_file;
    return;
  }
  if (!boost::filesystem::exists(stglauber_file)) {
    LOG(ERROR) << "input file does not exist: " << stglauber_file;
    return;
  }

  TFile refmult_in(refmult_file.c_str(), "READ");
  TFile sct_in(sct_file.c_str(), "READ");
  TFile stglauber_in(stglauber_file.c_str(), "READ");

  TH1D *refmult = (TH1D *)refmult_in.Get("refmult");
  TH2D *sct_npartncoll = (TH2D *)sct_in.Get("npartncoll_nominal");
  TH2D *stglauber_npartncoll = (TH2D *)stglauber_in.Get("hNcoll_Npart");

  double minmult = 100;
  double maxmult = 800;
  double npp = 2.00;
  double k = 0.80;
  double x = 0.14;
  double ppeff = 0.98;
  double auaueff = 0.84;
  double deff = ppeff - auaueff;
  double centmult = 540;
  double trigbias = 1.0;
  bool consteff = false;
  int nevents = 1e6;

  sct::NBDFit fitter(refmult, sct_npartncoll);
  fitter.minimumMultiplicityCut(minmult);
  fitter.useStGlauberChi2(true);
  fitter.useStGlauberNorm(true);
  fitter.setParameters(npp, k, x, ppeff, auaueff, centmult, trigbias, consteff);
  auto refit = fitter.fit(nevents);

  StNbdFitMaker maker;
  maker.SetParameters(npp, k, x, deff, trigbias, consteff);
  maker.SetMinimumMultiplicityCut(minmult);
  maker.ReadData(refmult_file.c_str(), stglauber_file.c_str(), "refmult");
  std::string tmpfilename = "/tmp/fit_test_stglauber_file.root";
  maker.Fit(nevents, tmpfilename.c_str());

  TFile stglauber_result(tmpfilename.c_str(), "READ");
  TH1D *stglauber_ref = (TH1D *)stglauber_result.Get("hRefMultSim");
  TH1D *stglauber_data = (TH1D *)stglauber_result.Get("refmult");

  EXPECT_LE(0.1, stglauber_ref->Chi2Test(refit->simu.get(), "UU NORM"));
}