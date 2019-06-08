/* With a full RefMultCorr definition its good to check a few things -
 * 1) the fraction of each events inside each centrality bin
 * 2) reweighting gives a flat ratio when compared to the glauber model
 *
 * apologies, this one requires quite a few command line arguments :)
 *
 */

#include "sct/centrality/refmultcorr_template.h"
#include "sct/lib/flags.h"
#include "sct/lib/logging.h"
#include "sct/lib/memory.h"
#include "sct/lib/string/string_cast.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/print_helper.h"

#include <fstream>
#include <iostream>
#include <string>

#include "boost/filesystem.hpp"

#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

// program settings
SCT_DEFINE_string(outDir, "tmp", "directory to save output to");
SCT_DEFINE_string(correctionFile, "refmultcorr_params.txt",
                  "results from pre_glauber_corrections");
SCT_DEFINE_string(glauberCentDefFile, "fit_results.txt",
                  "results from fit_glauber_to_data");
SCT_DEFINE_string(dataFile, "refmult.root",
                  "path to root file containing refmult tree");
SCT_DEFINE_string(treeName, "refMultTree", "refmult tree name");

SCT_DEFINE_string(refmultBranch, "refMult",
                  "name of reference multiplicty branch");
SCT_DEFINE_int(refmultMin, 0, "minimum refmult");
SCT_DEFINE_int(refmultMax, 800, "maximum refmult");

SCT_DEFINE_string(luminosityBranch, "lumi",
                  "name of branch containing luminosity information (zdc rate, "
                  "bbc rate, etc)");
SCT_DEFINE_double(lumiMin, 0, "minimum luminosity (kHz)");
SCT_DEFINE_double(lumiMax, 100, "maximum luminosity (kHz)");
SCT_DEFINE_double(lumiBins, 20, "number of bins to use for lumi correction");
SCT_DEFINE_double(lumiNorm, 0.0,
                  "normalization point - lumi correction = f(zdcNorm)/f(zdc)");

SCT_DEFINE_string(vzBranch, "vz", "name of vertex z branch");
SCT_DEFINE_double(vzMin, -30.0, "minimum Vz [cm]");
SCT_DEFINE_double(vzMax, 30.0, "maximum Vz [cm]");
SCT_DEFINE_double(vzBins, 20, "number of bins to use for Vz correction");
SCT_DEFINE_double(vzNorm, 0.0,
                  "normalization point - vz correction = f(vzNorm)/f(vz)");

SCT_DEFINE_string(vrBranch, "vr", "name of vertex r branch");
SCT_DEFINE_double(vrMax, 3.0, "maximum Vr");

SCT_DEFINE_string(dVzBranch, "vzvpdvz", "name of vz - vpd vz branch");
SCT_DEFINE_double(dVzMax, 3.0, "maximum dVz");

SCT_DEFINE_double(reweightingBound, 400.0,
                  "cutoff for reweighting - weight for any event with "
                  "refmultcorr above this threshold is 1");

SCT_DEFINE_string(glauberDistFile, "fit_results.root",
                  "file with the glauber refmult distribution histogram");
SCT_DEFINE_string(glauberDistName, "glauber",
                  "name of glauber refmult histogram");
SCT_DEFINE_string(glauberRatioName, "ratio_fit",
                  "name of glauber ratio histogram");

template <class T>
std::vector<T> GetParametersFromFile(const std::string &file,
                                     const std::string &key);
bool AcceptEvent(double vz, double vr, double dVz, unsigned refmult,
                 double lumi);

int main(int argc, char *argv[]) {
  // set help message and initialize logging and command line flags
  std::string usage =
      "Performs post-glauber crosschecks to see if the centrality definition "
      "is robust across centrality bins, luminosity bins, vz bins, etc";
  sct::SetUsageMessage(usage);

  sct::InitLogging(&argc, argv);
  sct::ParseCommandLineFlags(&argc, argv);

  // load the glauber file to get reference distributions
  TFile *glauber_file = new TFile(FLAGS_glauberDistFile.c_str(), "READ");
  if (!glauber_file->IsOpen()) {
    LOG(FATAL) << "Can't Open " << FLAGS_glauberDistFile;
  }
  TH1D *glauber_dist = (TH1D *)glauber_file->Get(FLAGS_glauberDistName.c_str());
  TH1D *glauber_ratio =
      (TH1D *)glauber_file->Get(FLAGS_glauberRatioName.c_str());

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

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_outDir);
  boost::filesystem::create_directories(dir);

  // parse the needed parameters for the glauber definition, vz, lumi
  // corrections, and weights
  std::vector<double> vz_par =
      GetParametersFromFile<double>(FLAGS_correctionFile, "vz parameters: ");
  std::vector<double> lumi_par =
      GetParametersFromFile<double>(FLAGS_correctionFile, "lumi parameters: ");
  std::vector<unsigned> cent_bounds = GetParametersFromFile<unsigned>(
      FLAGS_glauberCentDefFile, "nominal cent: ");
  std::vector<double> glauber_weights =
      GetParametersFromFile<double>(FLAGS_glauberCentDefFile, "weights: ");

  LOG(INFO) << "vz parameters: " << vz_par;
  if (vz_par.size() != 7)
    LOG(ERROR) << "incorrect number of vz parameters: expect 2, got "
               << vz_par.size();
  LOG(INFO) << "lumi parameters: " << lumi_par;
  if (lumi_par.size() != 2)
    LOG(ERROR) << "incorrect number of vz parameters: expect 7, got "
               << lumi_par.size();
  LOG(INFO) << "16 bin centrality def: " << cent_bounds;
  if (cent_bounds.size() != 16)
    LOG(ERROR) << "incorrect number of vz parameters: expect 16, got "
               << cent_bounds.size();
  LOG(INFO) << "glauber weights: " << glauber_weights;
  if (glauber_weights.size() != 7)
    LOG(ERROR) << "incorrect number of vz parameters: expect 7, got "
               << glauber_weights.size();

  // create our centrality definition
  sct::RefMultCorrTemplate centrality;
  centrality.setZDCParameters(lumi_par);
  centrality.setZDCNormalizationPoint(FLAGS_vzNorm);
  centrality.setVzParameters(vz_par);
  centrality.setVzNormalizationPoint(FLAGS_vzNorm);
  centrality.setCentralityBounds16Bin(cent_bounds);
  centrality.setWeightParameters(glauber_weights, FLAGS_reweightingBound);

  if (centrality.status() != true)
    LOG(FATAL) << "centrality definition not loaded properly, exiting";

  // ROOT setup
  // Histograms will calculate gaussian errors
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TH3::SetDefaultSumw2();

  // shut ROOT up :)
  gErrorIgnoreLevel = kWarning;

  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(true);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);

  // create different option sets for printing
  sct::histogramOpts hOpts;
  hOpts.marker_size = 0.5;

  sct::canvasOpts cOpts;

  sct::canvasOpts cOptsNoLeg;
  cOptsNoLeg.do_legend = false;

  sct::canvasOpts cOptsLogy;
  cOptsLogy.log_y = true;

  sct::canvasOpts cOptsLowerLeg;
  cOptsLowerLeg.leg_upper_bound = 0.4;
  cOptsLowerLeg.leg_left_bound = 0.2;
  cOptsLowerLeg.leg_lower_bound = 0.2;
  cOptsLowerLeg.leg_right_bound = 0.5;

  sct::canvasOpts cOptsLowerLegLogy;
  cOptsLowerLegLogy.leg_upper_bound = 0.45;
  cOptsLowerLegLogy.leg_left_bound = 0.2;
  cOptsLowerLegLogy.leg_lower_bound = 0.2;
  cOptsLowerLegLogy.leg_right_bound = 0.5;
  cOptsLowerLegLogy.log_y = true;

  // create histograms
  TH1D *raw_refmult =
      new TH1D("raw_refmult", ";refmult", FLAGS_refmultMax - FLAGS_refmultMin,
               FLAGS_refmultMin, FLAGS_refmultMax);
  TH1D *corr_refmult =
      new TH1D("corr_refmult", ";refmult", FLAGS_refmultMax - FLAGS_refmultMin,
               FLAGS_refmultMin, FLAGS_refmultMax);
  TH1D *weighted_refmult = new TH1D("weighted_refmult", ";refmult",
                                    FLAGS_refmultMax - FLAGS_refmultMin,
                                    FLAGS_refmultMin, FLAGS_refmultMax);
  TH3D *events_per_bin = new TH3D(
      "eventsperbin", ";zdcX;v_{z};cent;", FLAGS_lumiBins, FLAGS_lumiMin,
      FLAGS_lumiMax, FLAGS_vzBins, FLAGS_vzMin, FLAGS_vzMax, 16, -0.5, 15.5);
  TH3D *events_per_bin_weighted =
      new TH3D("eventsperbinweighted", ";zdcX;v_{z};cent;", FLAGS_lumiBins,
               FLAGS_lumiMin, FLAGS_lumiMax, FLAGS_vzBins, FLAGS_vzMin,
               FLAGS_vzMax, 16, -0.5, 15.5);

  while (reader.Next()) {
    double lumikhz = *lumi / 1000.0;

    if (!AcceptEvent(*vz, *vr, *dVz, *refmult, lumikhz))
      continue;

    centrality.setEvent(*refmult, lumikhz, *vz);
    double refmultcorr = centrality.refMultCorr();
    double weight = centrality.weight();
    int cent = centrality.centrality16();

    if (weight < 0 || cent < 0 || cent > 15)
      continue;

    raw_refmult->Fill(*refmult);
    corr_refmult->Fill(refmultcorr);
    weighted_refmult->Fill(refmultcorr, weight);
    events_per_bin->Fill(lumikhz, *vz, cent);
    events_per_bin_weighted->Fill(lumikhz, *vz, cent, weight);
  }

  // normalize all the refmult distributions to the glauber distribution
  int norm_bin_low = raw_refmult->GetXaxis()->FindBin(FLAGS_reweightingBound);
  int norm_bin_low_glauber =
      glauber_dist->GetXaxis()->FindBin(FLAGS_reweightingBound);

  raw_refmult->Scale(1.0 / raw_refmult->Integral(norm_bin_low, -1));
  corr_refmult->Scale(1.0 / corr_refmult->Integral(norm_bin_low, -1));
  weighted_refmult->Scale(1.0 / weighted_refmult->Integral(norm_bin_low, -1));
  glauber_dist->Scale(1.0 / glauber_dist->Integral(norm_bin_low, -1));

  sct::Overlay1D(glauber_dist, corr_refmult, "glauber", "unweighted", hOpts,
                 cOptsLogy, FLAGS_outDir, "glauber_unweighted_comparison", "",
                 "refmult", "arbitrary", "");
  sct::Overlay1D(glauber_dist, weighted_refmult, "glauber", "weighted", hOpts,
                 cOptsLogy, FLAGS_outDir, "glauber_weighted_comparison", "",
                 "refmult", "arbitrary", "");

  // normalize the event counts per event
  events_per_bin->Scale(1.0 / events_per_bin->Integral());
  events_per_bin_weighted->Scale(1.0 / events_per_bin_weighted->Integral());

  // plot the ratios
  TH1D *unweighted_ratio = (TH1D *)glauber_dist->Clone("unweighted_ref_ratio");
  TH1D *weighted_ratio = (TH1D *)glauber_dist->Clone("weighted_ref_ratio");
  unweighted_ratio->Divide(corr_refmult);
  unweighted_ratio->GetYaxis()->SetRangeUser(0, 3);
  weighted_ratio->Divide(weighted_refmult);
  weighted_ratio->GetYaxis()->SetRangeUser(0, 3);

  TF1 *weighted_ratio_fit =
      new TF1("w_ratio_fit", "[0]", 10, FLAGS_reweightingBound);
  weighted_ratio_fit->SetParameter(0, 1);
  weighted_ratio->Fit(weighted_ratio_fit, "ME", "", 10, FLAGS_reweightingBound);
  sct::PrettyPrint1D(unweighted_ratio, hOpts, cOpts, "unweighted ratio",
                     FLAGS_outDir, "unweighted_glauber_ratio", "", "refmult",
                     "glauber/refmultcorr", "");
  sct::PrettyPrint1D(weighted_ratio, hOpts, cOptsNoLeg, "weighted ratio",
                     FLAGS_outDir, "weighted_glauber_ratio", "", "refmult",
                     "glauber/refmultcorr", "");

  TH1D *events = events_per_bin->ProjectionZ();
  TH1D *events_weighted = events_per_bin_weighted->ProjectionZ();
  events->Scale(1.0 / events->Integral() * 8 / 10);
  events_weighted->Scale(1.0 / events_weighted->Integral() * 8 / 10);
  events->GetYaxis()->SetRangeUser(0.02, 0.08);
  events_weighted->GetYaxis()->SetRangeUser(0.02, 0.08);
  sct::Overlay1D(events, events_weighted, "unweighted", "weighted", hOpts,
                 cOptsLowerLeg, FLAGS_outDir, "centrality_bin_event_fraction",
                 "", "centrality bin [5%]", "fraction of events", "");

  // plot weighted event fraction as a function of zdcX
  std::vector<TH1D *> zdcx_event_fraction;
  std::vector<std::string> zdcx_event_fraction_name;
  for (int i = 1; i <= events_per_bin_weighted->GetXaxis()->GetNbins(); ++i) {
    std::string name = sct::MakeString(events_per_bin_weighted->GetName(), i);
    TH1D *tmp = events_per_bin_weighted->ProjectionZ(name.c_str(), i, i, 1, -1);
    tmp->Scale(1.0 / tmp->Integral() * 8 / 10);
    tmp->GetYaxis()->SetRangeUser(0.02, 0.08);
    std::string leg_name =
        sct::MakeString(std::setprecision(3),
                        events_per_bin_weighted->GetXaxis()->GetBinLowEdge(i),
                        " < zdcX [kHz] < ",
                        events_per_bin_weighted->GetXaxis()->GetBinUpEdge(i));
    zdcx_event_fraction.push_back(tmp);
    zdcx_event_fraction_name.push_back(leg_name);

    sct::Overlay1D(zdcx_event_fraction, zdcx_event_fraction_name, hOpts,
                 cOptsLowerLeg, FLAGS_outDir, "centrality_bin_event_fraction_zdcx",
                 "", "centrality bin [5%]", "fraction of events", "");
  }

  gflags::ShutDownCommandLineFlags();
  return 0;
}

template <class T>
std::vector<T> GetParametersFromFile(const std::string &file,
                                     const std::string &key) {
  std::ifstream input_file;
  input_file.open(file);
  std::string line;
  while (std::getline(input_file, line)) {
    if (sct::Consume(line, key)) {
      return sct::ParseArgStringToVec<T>(line, " ", false);
    }
  }
  return std::vector<T>();
}

bool AcceptEvent(double vz, double vr, double dVz, unsigned refmult,
                 double lumi) {
  if (vz < FLAGS_vzMin || vz > FLAGS_vzMax)
    return false;
  if (vr > FLAGS_vrMax)
    return false;
  if (dVz > FLAGS_dVzMax)
    return false;
  if (refmult < FLAGS_refmultMin || refmult > FLAGS_refmultMax)
    return false;
  if (lumi < FLAGS_lumiMin || lumi > FLAGS_lumiMax)
    return false;
  return true;
}