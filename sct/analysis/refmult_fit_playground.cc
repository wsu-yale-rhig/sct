/* QA routine for a single fit to test the normalization routine, fit and chi2
 * to show how the different parts of the fitting process work.
 *
 */

#include "sct/centrality/centrality.h"
#include "sct/centrality/multiplicity_model.h"
#include "sct/centrality/nbd_fit.h"
#include "sct/lib/enumerations.h"
#include "sct/lib/flags.h"
#include "sct/lib/logging.h"
#include "sct/lib/string/string_cast.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/random.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "boost/filesystem.hpp"

#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// program settings
SCT_DEFINE_string(outDir, "tmp",
                  "Path to directory to store output ROOT files");
SCT_DEFINE_string(outFile, "example", "output file name (no extension)");
SCT_DEFINE_string(
    glauberFile, "npartncoll.root",
    "path to root file containing the npart x ncoll distribution");
SCT_DEFINE_string(glauberHistName, "npartncoll",
                  "name of npart x ncoll histogram");
SCT_DEFINE_string(dataFile, "refmult.root",
                  "path to root file containing data refmult distribution");
SCT_DEFINE_string(dataHistName, "refmult",
                  "name of reference multiplicty histogram");
SCT_DEFINE_int(events, 1e6, "number of events per fit");

// model settings
SCT_DEFINE_double(npp, 2.38, "Npp value for negative binomial");
SCT_DEFINE_double(k, 2.0, "k value for negative binomial");
SCT_DEFINE_double(x, 0.13, "x value for two component multiplicity");
SCT_DEFINE_double(ppEfficiency, 0.98, "pp efficiency");
SCT_DEFINE_double(AuAuEfficiency, 0.84, "0-5% central AuAu efficiency");
SCT_DEFINE_int(centMult, 540, "average 0-5% central multiplicity");
SCT_DEFINE_bool(constEff, false, "turn on to use only pp efficiency");
SCT_DEFINE_bool(useStGlauberChi2, false,
                "use StGlauber Chi2 calculation instead of ROOT");
SCT_DEFINE_bool(useStGlauberNorm, true,
                "use StGlauber Normalization instead of integral norm");
SCT_DEFINE_double(trigBias, 1.0, "trigger bias");
SCT_DEFINE_int(minMult, 100,
               "minimum multiplicity for chi2 comparisons in fit");

int main(int argc, char *argv[]) {
  // shut ROOT up :)
  gErrorIgnoreLevel = kWarning;

  // set help message
  std::string usage = "Test macro for fitting and chi2.";
  sct::SetUsageMessage(usage);

  sct::InitLogging(&argc, argv);
  sct::ParseCommandLineFlags(&argc, argv);

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_outDir);
  boost::filesystem::create_directories(dir);

  // load the input files data refmult histogram & glauber npart x ncoll
  TFile *glauber_file = new TFile(FLAGS_glauberFile.c_str(), "READ");
  TFile *data_file = new TFile(FLAGS_dataFile.c_str(), "READ");

  if (!glauber_file->IsOpen()) {
    LOG(ERROR) << "Glauber input file could not be opened: "
               << FLAGS_glauberFile << " not found or corrupt";
    return 1;
  }
  if (!data_file->IsOpen()) {
    LOG(ERROR) << "Data input file could not be opened: " << FLAGS_dataFile
               << " not found or corrupt";
    return 1;
  }

  // load data refmult histogram & glauber npart x ncoll
  TH2D *npartncoll = (TH2D *)glauber_file->Get(FLAGS_glauberHistName.c_str());
  TH1D *refmult = (TH1D *)data_file->Get(FLAGS_dataHistName.c_str());

  if (npartncoll == nullptr) {
    LOG(ERROR) << "NPart x NColl histogram: " << FLAGS_glauberHistName
               << " not found in file: " << FLAGS_glauberFile;
    return 1;
  }
  if (refmult == nullptr) {
    LOG(ERROR) << "Reference multiplicity histogram: " << FLAGS_dataHistName
               << " not found in file: " << FLAGS_dataFile;
    return 1;
  }

  // we will start the analysis routine - first, print out selected parameters
  LOG(INFO)
      << "Fitting refmult distribution with Glauber. Selected parameters:";
  LOG(INFO) << "Npp: " << FLAGS_npp << " K: " << FLAGS_k << " X: " << FLAGS_x;
  LOG(INFO) << "pp efficiency: " << FLAGS_ppEfficiency;
  LOG(INFO) << "AuAu efficiency: " << FLAGS_AuAuEfficiency;
  LOG(INFO) << "AuAu central multiplicity (for efficiency correction): "
            << FLAGS_centMult;
  LOG(INFO) << "trigger bias: " << FLAGS_trigBias;
  LOG(INFO) << "constant efficiency: " << (FLAGS_constEff ? "true" : "false");

  // create our fitting model
  LOG(INFO) << "First, fit with the sct::NBDFit tool";
  sct::NBDFit fitter(refmult, npartncoll);
  fitter.minimumMultiplicityCut(FLAGS_minMult);
  fitter.useStGlauberChi2(FLAGS_useStGlauberChi2);
  fitter.useStGlauberNorm(FLAGS_useStGlauberNorm);
  fitter.setParameters(FLAGS_npp, FLAGS_k, FLAGS_x, FLAGS_ppEfficiency,
                       FLAGS_AuAuEfficiency, FLAGS_centMult, FLAGS_trigBias,
                       FLAGS_constEff);
  auto refit = fitter.fit(FLAGS_events);
  LOG(INFO) << "Resulting chi2/ndf: " << refit->chi2 / refit->ndf;
  LOG(INFO) << "Histogram normalization in fitter: ";
  LOG(INFO) << "Integral for data: " << refit->data->Integral();
  LOG(INFO) << "Integral for glauber: " << refit->simu->Integral();
  LOG(INFO) << "Integral for data in norm range: "
            << refit->data->Integral(refit->data->GetXaxis()->FindBin(FLAGS_minMult),
                                 refit->data->GetNbinsX());
  LOG(INFO) << "Integral for glauber in norm range: "
            << refit->simu->Integral(
                   refit->simu->GetXaxis()->FindBin(FLAGS_minMult),
                   refit->simu->GetNbinsX());

  // now create the multiplicity model - do the same procedure as the fitter
  // by hand
  LOG(INFO) << "Perform Glauber -> refmult calculation by hand using "
               "sct::MultiplicityModel";
  sct::MultiplicityModel mult_model(
      FLAGS_npp, FLAGS_k, FLAGS_x, FLAGS_ppEfficiency, FLAGS_AuAuEfficiency,
      FLAGS_centMult, FLAGS_trigBias, FLAGS_constEff);

  TH1D *glauber_refmult = (TH1D *)refit->simu->Clone("glauber_ref");
  glauber_refmult->Reset();

  // loop over the glauber NPart x NColl distribution, fill the refmult
  // distribution
  for (int i = 0; i < FLAGS_events; ++i) {
    double npart = 0;
    double ncoll = 0;
    npartncoll->GetRandom2(npart, ncoll);
    if (npart < 2 || ncoll < 1)
      continue;
    double mult = mult_model.multiplicity(npart, ncoll);
    glauber_refmult->Fill(mult);
  }

  // Show that the fit is similar between the fitter and the multiplicity model
  auto simu_compare = fitter.chi2(refit->simu.get(), glauber_refmult);
  LOG(INFO) << "comparison of fitter and by-hand multiplicity model chi2/ndf: "
            << simu_compare.first / simu_compare.second;

  auto data_compare = fitter.chi2(refmult, glauber_refmult);
  LOG(INFO) << "comparison of data and by-hand multiplicity model chi2/ndf: "
            << data_compare.first / data_compare.second;

  TH1D *glauber_ref_scaled =
      (TH1D *)glauber_refmult->Clone("glauber_ref_scaled");

  // Now compare normalization
  LOG(INFO) << "Integral for data: " << refmult->Integral();
  LOG(INFO) << "Integral for glauber: " << glauber_refmult->Integral();
  LOG(INFO) << "Integral for data in norm range: "
            << refmult->Integral(refmult->GetXaxis()->FindBin(FLAGS_minMult),
                                 refmult->GetNbinsX());
  LOG(INFO) << "Integral for glauber in norm range: "
            << glauber_refmult->Integral(
                   glauber_refmult->GetXaxis()->FindBin(FLAGS_minMult),
                   glauber_refmult->GetNbinsX());

  double norm = fitter.norm(refmult, glauber_refmult);
  LOG(INFO) << "relative normalization: " << norm;
  glauber_refmult->Scale(norm);

  LOG(INFO) << "Post normalization comparison: ";
  LOG(INFO) << "Integral for data: " << refmult->Integral();
  LOG(INFO) << "Integral for glauber: " << glauber_refmult->Integral();
  LOG(INFO) << "Integral for data in norm range: "
            << refmult->Integral(refmult->GetXaxis()->FindBin(FLAGS_minMult),
                                 refmult->GetNbinsX());
  LOG(INFO) << "Integral for glauber in norm range: "
            << glauber_refmult->Integral(
                   glauber_refmult->GetXaxis()->FindBin(FLAGS_minMult),
                   glauber_refmult->GetNbinsX());
  
  // save results to file
  std::string out_name = FLAGS_outDir + "/" + FLAGS_outFile + ".root";
  TFile out_file(out_name.c_str(), "RECREATE");

  refmult->SetName("refmult");
  refmult->Write();

  npartncoll->SetName("npartncoll");
  npartncoll->Write();

  refit->simu->SetName("fitter_ref");
  refit->simu->Write();

  glauber_refmult->SetName("hand_ref");
  glauber_refmult->Write();

  out_file.Close();

      gflags::ShutDownCommandLineFlags();
  return 0;
}
