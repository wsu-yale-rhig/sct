/* Used to create a set of refmult vs tofmatch correlations and select cuts
 * based on the correlations to remove pileup events. In general the correlation
 * between refmult and tofmult should be linear and deviation from that linear
 * correlation can be taken as a signal of out-of-time pileup, when refmult can
 * be increased without a corresponding increase in the tofMult.
 */

#include "sct/lib/flags.h"
#include "sct/lib/logging.h"
#include "sct/lib/memory.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/print_helper.h"

#include <fstream>
#include <iostream>
#include <string>

#include "boost/filesystem.hpp"

#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// program settings
SCT_DEFINE_string(outDir, "tmp",
                  "Path to directory to store output ROOT file and QA");
SCT_DEFINE_string(outFile, "tofmult_mask.root", "output ROOT file name");
SCT_DEFINE_string(dataFile, "refmult.root",
                  "path to root file containing refmult tree");
SCT_DEFINE_string(treeName, "refMultTree", "refmult tree name");

SCT_DEFINE_string(refmultBranch, "refMult",
                  "name of reference multiplicty branch");
SCT_DEFINE_int(refmultMin, 0, "minimum refmult");
SCT_DEFINE_int(refmultMax, 800, "maximum refmult");
SCT_DEFINE_int(refmultProjectionWidth, 1,
               "number of refmult bins per projection (for fit)");

SCT_DEFINE_string(tofmatchBranch, "tofMatch",
                  "branch containing number of tracks matched to tof");
SCT_DEFINE_int(tofmatchMin, 0, "minimum tofMatch");
SCT_DEFINE_int(tofmatchMax, 800, "maximum tofMatch");

SCT_DEFINE_string(luminosityBranch, "lumi",
                  "name of branch containing luminosity information (zdc rate, "
                  "bbc rate, etc)");
SCT_DEFINE_double(lumiMin, 0, "minimum luminosity [kHz]");
SCT_DEFINE_double(lumiMax, 100, "maximum luminosity [kHz]");
SCT_DEFINE_double(lumiBins, 1,
                  "number of bins in lumi to use for tofMatch mask");

SCT_DEFINE_string(vzBranch, "vz", "name of vertex z branch");
SCT_DEFINE_double(vzMin, -30.0, "minimum Vz [cm]");
SCT_DEFINE_double(vzMax, 30.0, "maximum Vz [cm]");

SCT_DEFINE_string(vrBranch, "vr", "name of vertex r branch");
SCT_DEFINE_double(vrMax, 3.0, "maximum Vr [cm]");

SCT_DEFINE_string(dVzBranch, "vzvpdvz", "name of vz - vpd vz branch");
SCT_DEFINE_double(dVzMax, 3.0, "maximum dVz [cm]");

SCT_DEFINE_int(minEventsForFit, 100,
               "minimum number of events required in a projection for a fit to "
               "be attempted");
SCT_DEFINE_double(nSigma, 4.0, "controls the width of the mask");

bool AcceptEvent(double vz, double vr, double dVz, unsigned refmult,
                 double lumi, unsigned tofMatch);

// used to create the fit
struct fitParams {
  double norm_default = 0.0;
  double mean_default = 0.0;
  double width_default = 0.0;
  double refmult_bin_min = 0.0;
  double refmult_bin_max = 0.0;
  double lumi_bin_min = 0.0;
  double lumi_bin_max = 0.0;
};

// contains the fit information for a single tofMatch fit
struct fit {
  fitParams params;
  std::unique_ptr<TH1D> projection;
  std::unique_ptr<TF1> fit;
  std::string bin_name;
  bool success = false;

  double norm() {
    if (fit != nullptr)
      return fit->GetParameter(0);
    return 0;
  }

  double mean() {
    if (fit != nullptr)
      return fit->GetParameter(1);
    return 0;
  }

  double width() {
    if (fit != nullptr)
      return fit->GetParameter(2);
    return 0;
  }
};

void FitTofMult(std::vector<std::vector<fit>> &results, TH3D *hist);
fit FitTofMultSlice(fitParams pars, TH3D *hist);
void FillMask(std::vector<std::vector<fit>> &results, TH3D *mask);

int main(int argc, char *argv[]) {
  // set help message and initialize logging and command line flags
  std::string usage =
      "Creates a set of refmult x tofMatch pileup cuts by making a TH3D mask.";
  usage += "Input must have refmult, luminosity, tofMatch, vz, vr and dVz "
           "branches - output is a 3D refmult x tofMatch x luminosity cut.";
  sct::SetUsageMessage(usage);

  sct::InitLogging(&argc, argv);
  sct::ParseCommandLineFlags(&argc, argv);

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_outDir);
  boost::filesystem::create_directories(dir);

  // load the input file & load tree
  TFile *in_file = new TFile(FLAGS_dataFile.c_str(), "READ");
  if (!in_file || !in_file->IsOpen()) {
    LOG(FATAL) << "refmult file: " << FLAGS_dataFile << " could not be opened";
  }

  TTreeReader reader(FLAGS_treeName.c_str(), in_file);
  TTreeReaderValue<unsigned> refmult(reader, FLAGS_refmultBranch.c_str());
  TTreeReaderValue<double> vz(reader, FLAGS_vzBranch.c_str());
  TTreeReaderValue<double> lumi(reader, FLAGS_luminosityBranch.c_str());
  TTreeReaderValue<double> vr(reader, FLAGS_vrBranch.c_str());
  TTreeReaderValue<double> dVz(reader, FLAGS_dVzBranch.c_str());
  TTreeReaderValue<unsigned> tofMatch(reader, FLAGS_tofmatchBranch.c_str());

  // make an output file to store the corrected refmult distribution
  std::string out_file_name = FLAGS_outDir + "/" + FLAGS_outFile;
  TFile *out_file = new TFile(out_file_name.c_str(), "RECREATE");
  out_file->cd(0);
  // ROOT setup
  // -----------------------------------------
  // Histograms will calculate gaussian errors
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  // shut ROOT up :)
  gErrorIgnoreLevel = kWarning;

  // create histograms - one to be filled now which we will take projections of
  // for fitting, and one for the eventual mask

  unsigned refmult_bins = FLAGS_refmultMax - FLAGS_refmultMin;
  unsigned tofmatch_bins = FLAGS_tofmatchMax - FLAGS_tofmatchMin;

  TH3D *ref_tofmatch_lumi = new TH3D(
      "refmulttofmatchlumi", ";refmult;tofMatch;zdcX [kHz]", refmult_bins,
      FLAGS_refmultMin, FLAGS_refmultMax, tofmatch_bins, FLAGS_tofmatchMin,
      FLAGS_tofmatchMax, FLAGS_lumiBins, FLAGS_lumiMin, FLAGS_lumiMax);
  TH3D *tofmatch_mask = new TH3D(
      "tofmatch_mask", ";refmult;tofMatch;zdcX [kHz]", refmult_bins,
      FLAGS_refmultMin, FLAGS_refmultMax, tofmatch_bins, FLAGS_tofmatchMin,
      FLAGS_tofmatchMax, FLAGS_lumiBins, FLAGS_lumiMin, FLAGS_lumiMax);

  // loop over all events and fill the initial histogram
  while (reader.Next()) {
    double lumikhz = *lumi / 1000.0;

    if (!AcceptEvent(*vz, *vr, *dVz, *refmult, lumikhz, *tofMatch))
      continue;

    ref_tofmatch_lumi->Fill(*refmult, *tofMatch, lumikhz);
  }

  // perform the 1D fits that the mask will be created from. This is done by
  // taking projections for each luminosity bin and refmult bin and fitting the
  // tofMatch distribution with a gaussian.
  std::vector<std::vector<fit>> fits;
  FitTofMult(fits, ref_tofmatch_lumi);
  FillMask(fits, tofmatch_mask);

  out_file->cd();
  ref_tofmatch_lumi->Write();
  tofmatch_mask->Write();
  // for(auto& a : fits)
  //   for(auto& b :a)
  //     b.projection->Write();
  out_file->Close();

  gflags::ShutDownCommandLineFlags();
  return 0;
}

bool AcceptEvent(double vz, double vr, double dVz, unsigned refmult,
                 double lumi, unsigned tofMatch) {
  if (vz < FLAGS_vzMin || vz > FLAGS_vzMax)
    return false;
  if (vr > FLAGS_vrMax)
    return false;
  if (fabs(dVz) > FLAGS_dVzMax)
    return false;
  if (refmult < FLAGS_refmultMin || refmult > FLAGS_refmultMax)
    return false;
  if (lumi < FLAGS_lumiMin || lumi > FLAGS_lumiMax)
    return false;
  if (tofMatch < FLAGS_tofmatchMin || tofMatch > FLAGS_tofmatchMax)
    return false;
  return true;
}

void FitTofMult(std::vector<std::vector<fit>> &results, TH3D *hist) {

  if (hist->GetXaxis()->GetNbins() % FLAGS_refmultProjectionWidth != 0) {
    LOG(ERROR) << "Can not evenly divide " << hist->GetXaxis()->GetNbins()
               << " refmult bins by the requested width of "
               << FLAGS_refmultProjectionWidth;
    return;
  }
  unsigned refmult_bins =
      hist->GetXaxis()->GetNbins() / FLAGS_refmultProjectionWidth;

  results.clear();
  results.resize(refmult_bins);
  for (auto &bin : results)
    bin.resize(FLAGS_lumiBins);

  // loop over every refmult x luminosity bin and fit the projection
  for (int i = 0; i < refmult_bins; ++i) {
    int ref_bin_low = i * FLAGS_refmultProjectionWidth + 1;
    int ref_bin_high = (i + 1) * FLAGS_refmultProjectionWidth;
    for (int j = 0; j < FLAGS_lumiBins; ++j) {
      int lumi_bin_low = j + 1;
      int lumi_bin_high = j + 1;

      // create the default parameters for this projection
      fitParams pars;
      pars.norm_default = 0.01;
      pars.mean_default = hist->GetXaxis()->GetBinCenter(ref_bin_high);
      pars.width_default = pars.mean_default / 5.0;
      pars.refmult_bin_min = ref_bin_low;
      pars.refmult_bin_max = ref_bin_high;
      pars.lumi_bin_min = lumi_bin_low;
      pars.lumi_bin_max = lumi_bin_high;

      results[i][j] = FitTofMultSlice(pars, hist);
    }
  }
}

fit FitTofMultSlice(fitParams pars, TH3D *hist) {
  fit result;

  // create bin name
  std::string bin_name = sct::MakeString(
      "refmult_", hist->GetXaxis()->GetBinLowEdge(pars.refmult_bin_min), "_",
      hist->GetXaxis()->GetBinUpEdge(pars.refmult_bin_max), "_lumi_",
      hist->GetZaxis()->GetBinLowEdge(pars.lumi_bin_min), "_",
      hist->GetZaxis()->GetBinUpEdge(pars.lumi_bin_max));
  std::string projection_name = "projection_" + bin_name;
  std::string fit_name = "fit_" + bin_name;

  result.params = pars;
  result.bin_name = bin_name;
  result.projection = std::unique_ptr<TH1D>((TH1D *)hist->ProjectionY(
      projection_name.c_str(), pars.refmult_bin_min, pars.refmult_bin_max,
      pars.lumi_bin_min, pars.lumi_bin_max));
  result.projection->SetDirectory(0);
  result.fit = std::make_unique<TF1>(fit_name.c_str(), "gaus(0)",
                                     FLAGS_tofmatchMin, FLAGS_tofmatchMax);

  result.fit->SetParameter(0, pars.norm_default);
  result.fit->SetParameter(1, pars.mean_default);
  result.fit->SetParameter(2, pars.width_default);

  if (result.projection->Integral() < FLAGS_minEventsForFit)
    return result;

  result.projection->Scale(1.0 / result.projection->Integral());
  auto fit_result = result.projection->Fit(result.fit.get(), "EMS");
  result.success = true;

  return result;
}

void FillMask(std::vector<std::vector<fit>> &results, TH3D *mask) {
  mask->Reset();

  for (int i = 1; i <= mask->GetXaxis()->GetNbins(); ++i) {
    int refmult_bin = (i-1) / FLAGS_refmultProjectionWidth;
    for (int j = 1; j <= mask->GetYaxis()->GetNbins(); ++j) {
      for (int k = 1; k <= mask->GetZaxis()->GetNbins(); ++k) {
        int lumi_bin = k - 1;

        fit& result = results[refmult_bin][lumi_bin];

        if (!result.success)
          continue;
        
        if (result.width() > 30)
          continue;

        // get the widths from the gaussian
        double min = result.mean() - result.width() * FLAGS_nSigma;
        double max = result.mean() + result.width() * FLAGS_nSigma;
        double tofmatch_val = mask->GetYaxis()->GetBinLowEdge(j);
        bool val = tofmatch_val >= min && tofmatch_val <= max;

        mask->SetBinContent(i, j, k, val);
      }
    }
  }
}