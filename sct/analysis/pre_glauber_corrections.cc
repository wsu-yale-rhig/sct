/* Refmult can depend on detector conditions - for instance,
 * as luminosity increases, average tracking efficiency decreases.
 * Before comparing a refmult distribution to glauber to extract
 * a centrality definition, the refmult distribution should be corrected
 * for these effects. This routine corrects for luminosity dependence
 * and Vz dependence using the following functional forms:
 *
 * correction to <refmult> as a function of luminosity
 * fit <refmult> vs luminosity with a first order polynomial, then:
 * corrected refmult = refmult * f(zdcNorm) / f(refmult)
 *
 * correction to refmult as a function of Vz: the high multiplicity tail
 * of the refmult distribution is fit with an error function in multiple bins:
 * vz fit = [0] + [0]*erf([1]*(x-[2]))
 * parameter 2 is then plotted as a function of Vz, and fit with a 6th order
 * polynomial.
 * corrected refmult = refmult * f(vzNorm) / f(vz)
 *
 * The routines produce an output tree with fully corrected refmult, as well as
 * some QA histograms for the corrections.
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
SCT_DEFINE_string(outFile, "corrected_refmult.root", "output ROOT file name");
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
SCT_DEFINE_double(lumiBins, 60,
                  "number of bins to use for luminosity correction");
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

// applies event cuts as defined in program flags for each stage of the
// corrections
bool AcceptEvent(double vz, double vr, double dVz, unsigned refmult,
                 double lumi);

// used internally to bin and project TH3s along different axes
template <class T>
using result_container =
    std::pair<std::vector<std::shared_ptr<T>>, std::vector<std::string>>;

result_container<TProfile>
ProfileRefmultInBinsOfVz(TH3D *hist, std::string name_prefix, int nsplits = 4);

result_container<TH1D>
ProjectRefMultInBinsOfVz(TH3D *hist, std::string name_prefix, int nsplits = 4);

result_container<TH1D>
ProjectRefMultInBinsOfZDC(TH3D *hist, std::string name_prefix, int nsplits = 4);

// used to print out QA plots at each step of the corrections
void RefmultQA(TH3D *h, std::string name, sct::histogramOpts hOpts_hist,
               sct::canvasOpts cOpts_hist, sct::histogramOpts hOpts_prof,
               sct::canvasOpts cOpts_prof, int bins = 4, bool fit_flat = true);

int main(int argc, char *argv[]) {
  // set help message and initialize logging and command line flags
  std::string usage =
      "Performs luminosity and Vz corrections to a refmult distribution tree. ";
  usage +=
      "Input must have refmult, luminosity & vz branches - output is a second "
      "tree with corrected refmult.";
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

  // make an output file to store the corrected refmult distribution
  std::string out_file_name = FLAGS_outDir + "/" + FLAGS_outFile;
  TFile *out_file = new TFile(out_file_name.c_str(), "RECREATE");

  // fit and histogram parameters
  // -----------------------------------------
  // user ranges for projections/fits
  double avg_refmult_y_min = 120;
  double avg_refmult_y_max = 220;
  int avg_refmult_rebinning_lumi = 4;

  double vz_tail_min = 450;
  double vz_tail_max = 650;

  double vz_pol_min = FLAGS_vzMin;
  double vz_pol_max = FLAGS_vzMax;

  double h_param_y_min = 500;
  double h_param_y_max = 600;

  // Fit function definitions and default parameter values
  // -----------------------------------------
  // avg refmult as a function of luminosity (linear function) defaults
  std::string lumi_function = "[0]+[1]*x";
  double lumi_func_p0_default = 200.0;
  double lumi_func_p1_default = -0.1;

  // refmult tail distribution as a function of vz (error function) defaults
  std::string vz_function = "[0] + [0] * TMath::Erf([1]*(x-[2]))";
  double vz_func_p0_default = 0.0003;
  double vz_func_p1_default = -0.014;
  double vz_func_p2_default = 520;
  // double vz_func_p0_default = 0.0003;
  // double vz_func_p1_default = -0.02;
  // double vz_func_p2_default = 520;

  // erf param as a function of vz (pol6) defaults
  std::string vz_pol_function = "pol6(0)";
  double vz_pol_func_p0_default = 550;
  double vz_pol_func_p1_default = 0.1;
  double vz_pol_func_p2_default = 0.1;
  double vz_pol_func_p3_default = 0.1;
  double vz_pol_func_p4_default = 0.1;
  double vz_pol_func_p5_default = 0.1;
  double vz_pol_func_p6_default = 0.1;

  // ROOT setup
  // -----------------------------------------
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

  // UNCORRECTED REFMULT - GENERATE LUMINOSITY CORRECTION
  // --------------------------------------------------------------------------

  // event loop
  // first calculate average refmult as a function of luminosity
  TH3D *uncorr_lumi =
      new TH3D("uncorr_lumi", ";V_{z};luminosity [kHz];refmult", FLAGS_vzBins,
               FLAGS_vzMin, FLAGS_vzMax, FLAGS_lumiBins, FLAGS_lumiMin,
               FLAGS_lumiMax, (FLAGS_refmultMax - FLAGS_refmultMin),
               FLAGS_refmultMin, FLAGS_refmultMax);

  while (reader.Next()) {
    double lumikhz = *lumi / 1000.0;

    if (!AcceptEvent(*vz, *vr, *dVz, *refmult, lumikhz))
      continue;

    uncorr_lumi->Fill(*vz, lumikhz, *refmult);
  }

  LOG(INFO) << "Fitting average refmult as a function of luminosity";
  TProfile *uncorr_lumi_1d = (TProfile *)((TH2D *)uncorr_lumi->Project3D("ZY"))
                                 ->ProfileX("uncorr_profile");
  uncorr_lumi_1d->Rebin(avg_refmult_rebinning_lumi);

  // now fit with a first order polynomial for our luminosity correction
  TF1 *luminosity_correction =
      new TF1("luminosity_correction", lumi_function.c_str(), FLAGS_lumiMin,
              FLAGS_lumiMax);
  luminosity_correction->SetParameters(lumi_func_p0_default,
                                       lumi_func_p1_default);
  uncorr_lumi_1d->Fit(luminosity_correction);

  // print uncorrected average refmult
  uncorr_lumi_1d->GetYaxis()->SetRangeUser(avg_refmult_y_min,
                                           avg_refmult_y_max);
  sct::PrettyPrint1D(uncorr_lumi_1d, hOpts, cOptsNoLeg, "", FLAGS_outDir,
                     "uncorrected_avg_refmult_vs_lumi", "", "ZDC Rate [kHz]",
                     "<refmult>", "");
  RefmultQA(uncorr_lumi, "uncorrected", hOpts, cOptsLowerLegLogy, hOpts, cOpts,
            4, false);

  // LUMINOSITY CORRECTED - GENERATE VZ CORRECTION
  // --------------------------------------------------------------------------
  // loop over the data again, this time scale each entry by the luminosity
  // correction
  LOG(INFO) << "generating luminosity corrected refmult distribution";
  reader.Restart();
  TH3D *corr_lumi =
      new TH3D("corr_lumi", ";V_{z};luminosity [kHz];refmult", FLAGS_vzBins,
               FLAGS_vzMin, FLAGS_vzMax, FLAGS_lumiBins, FLAGS_lumiMin,
               FLAGS_lumiMax, (FLAGS_refmultMax - FLAGS_refmultMin),
               FLAGS_refmultMin, FLAGS_refmultMax);

  while (reader.Next()) {
    double lumikhz = *lumi / 1000.0;

    if (!AcceptEvent(*vz, *vr, *dVz, *refmult, lumikhz))
      continue;

    double norm = luminosity_correction->Eval(FLAGS_lumiNorm);
    double zdc_corr = luminosity_correction->Eval(lumikhz);

    double corr_refmult = static_cast<double>(*refmult) * norm / zdc_corr;

    corr_lumi->Fill(*vz, lumikhz, corr_refmult);
  }

  // now generate a 1D profile, <refmult> vs luminosity, and fit to see if its
  // approximately flat
  TProfile *corr_lumi_1d = (TProfile *)((TH2D *)corr_lumi->Project3D("ZY"))
                               ->ProfileX("corr_profile");
  corr_lumi_1d->Rebin(avg_refmult_rebinning_lumi);
  corr_lumi_1d->GetYaxis()->SetRangeUser(avg_refmult_y_min, avg_refmult_y_max);
  TF1 *corr_ref_lumi_fit =
      new TF1("lumi_corrected_ref_lumi_fit", lumi_function.c_str(),
              FLAGS_lumiMin, FLAGS_lumiMax);
  corr_lumi_1d->Fit(corr_ref_lumi_fit);

  sct::PrettyPrint1D(corr_lumi_1d, hOpts, cOptsNoLeg, "", FLAGS_outDir,
                     "lumi_corrected_avg_refmult_vs_lumi", "", "ZDC Rate [kHz]",
                     "<refmult>", "");
  RefmultQA(corr_lumi, "lumi_corrected", hOpts, cOptsLowerLegLogy, hOpts,
            cOpts);

  // now, do vz correction
  LOG(INFO) << "Performing Vz correction: extracting fit parameter h";
  std::vector<TH1D *> vz_uncorr_binned_refmult;
  std::vector<TF1 *> vz_uncorr_binned_refmult_fit;
  TH1D *uncorr_vz_fit_h_param =
      new TH1D("uncorr_vz_fit_h_param", ";v_{Z};h", FLAGS_vzBins, FLAGS_vzMin,
               FLAGS_vzMax);

  // we will print out each tail fit
  boost::filesystem::path tail_fit_dir(FLAGS_outDir);
  tail_fit_dir /= "lumi_corrected_tail_fit";
  boost::filesystem::create_directories(tail_fit_dir);
  std::string tail_fit_dir_path = tail_fit_dir.string();

  for (int i = 1; i <= FLAGS_vzBins; ++i) {
    vz_uncorr_binned_refmult.push_back(corr_lumi->ProjectionZ(
        sct::MakeString("refmult_vz_bin_", i).c_str(), i, i, 1, -1));
    vz_uncorr_binned_refmult[i - 1]->Scale(
        1.0 / vz_uncorr_binned_refmult[i - 1]->Integral());
    vz_uncorr_binned_refmult_fit.push_back(
        new TF1(sct::MakeString("refmult_vz_bin_fit_", i).c_str(),
                vz_function.c_str(), vz_tail_min, vz_tail_max));
    vz_uncorr_binned_refmult_fit[i - 1]->SetParameters(
        vz_func_p0_default, vz_func_p1_default, vz_func_p2_default);
    vz_uncorr_binned_refmult[i - 1]->Fit(vz_uncorr_binned_refmult_fit[i - 1],
                                         "ME", "", vz_tail_min, vz_tail_max);
    uncorr_vz_fit_h_param->SetBinContent(
        i, vz_uncorr_binned_refmult_fit[i - 1]->GetParameter(2));
    uncorr_vz_fit_h_param->SetBinError(
        i, vz_uncorr_binned_refmult_fit[i - 1]->GetParError(2));

    vz_uncorr_binned_refmult[i - 1]->GetXaxis()->SetRangeUser(vz_tail_min,
                                                              vz_tail_max);
    std::string hist_title = sct::MakeString(
        std::setprecision(3), corr_lumi->GetXaxis()->GetBinLowEdge(i),
        " < v_{z} < ", corr_lumi->GetXaxis()->GetBinUpEdge(i));
    gStyle->SetOptTitle(true);
    sct::PrettyPrint1D(vz_uncorr_binned_refmult[i - 1], hOpts, cOptsNoLeg, "",
                       tail_fit_dir_path,
                       sct::MakeString("lumi_corrected_tail_fit_h_", i),
                       hist_title, "refMult", "fit paramter h", "");
    gStyle->SetOptTitle(false);
  }
  LOG(INFO) << "Fitting h parameter distribution with sixth order polynomial";
  TF1 *vz_correction = new TF1("uncorr_vz_fit_param_fit",
                               vz_pol_function.c_str(), vz_pol_min, vz_pol_max);
  vz_correction->SetParameters(vz_pol_func_p0_default, vz_pol_func_p1_default,
                               vz_pol_func_p2_default, vz_pol_func_p3_default,
                               vz_pol_func_p4_default, vz_pol_func_p5_default,
                               vz_pol_func_p6_default);

  uncorr_vz_fit_h_param->Fit(vz_correction, "MEF");
  uncorr_vz_fit_h_param->GetYaxis()->SetRangeUser(h_param_y_min, h_param_y_max);
  sct::PrettyPrint1D(uncorr_vz_fit_h_param, hOpts, cOptsNoLeg, "", FLAGS_outDir,
                     "lumi_corrected_fit_h", "", "v_{z} [cm]", "fit paramter h",
                     "");

  // LUMINOSITY & VZ CORRECTED - GENERATE FULLY CORRECTED REFMULT
  // --------------------------------------------------------------------------

  // create luminosity & vz corrected refmult
  LOG(INFO) << "generating luminosity and Vz corrected refmult distribution";
  reader.Restart();
  TH3D *corr_lumi_vz =
      new TH3D("corr_lumi_vz", ";V_{z};luminosity [kHz];refmult", FLAGS_vzBins,
               FLAGS_vzMin, FLAGS_vzMax, FLAGS_lumiBins, FLAGS_lumiMin,
               FLAGS_lumiMax, (FLAGS_refmultMax - FLAGS_refmultMin),
               FLAGS_refmultMin, FLAGS_refmultMax);

  while (reader.Next()) {
    double lumikhz = *lumi / 1000.0;

    if (!AcceptEvent(*vz, *vr, *dVz, *refmult, lumikhz))
      continue;

    double zdc_norm = luminosity_correction->Eval(FLAGS_lumiNorm);
    double zdc_corr = luminosity_correction->Eval(lumikhz);

    double vz_norm = vz_correction->Eval(FLAGS_vzNorm);
    double vz_corr = vz_correction->Eval(*vz);

    double corr_refmult = static_cast<double>(*refmult) * (zdc_norm * vz_norm) /
                          (zdc_corr * vz_corr);

    corr_lumi_vz->Fill(*vz, lumikhz, corr_refmult);
  }

  RefmultQA(corr_lumi_vz, "lumi_vz_corrected", hOpts, cOptsLowerLegLogy, hOpts,
            cOpts, 4);

  // check vz h distribution
  LOG(INFO) << "Checking Vz correction: extracting fit parameter h";
  std::vector<TH1D *> vz_corr_binned_refmult;
  std::vector<TF1 *> vz_corr_binned_refmult_fit;
  TH1D *corr_vz_fit_h_param = new TH1D("corr_vz_fit_h_param", ";v_{Z};h",
                                       FLAGS_vzBins, FLAGS_vzMin, FLAGS_vzMax);

  // we will print out each corrected tail fit
  boost::filesystem::path corrected_tail_fit_dir(FLAGS_outDir);
  corrected_tail_fit_dir /= "lumi_vz_corrected_tail_fit";
  boost::filesystem::create_directories(corrected_tail_fit_dir);
  std::string corrected_tail_fit_dir_path = corrected_tail_fit_dir.string();

  for (int i = 1; i <= FLAGS_vzBins; ++i) {
    vz_corr_binned_refmult.push_back(corr_lumi_vz->ProjectionZ(
        sct::MakeString("refmult_vz_bin2_", i).c_str(), i, i, 1, -1));
    vz_corr_binned_refmult[i - 1]->Scale(
        1.0 / vz_corr_binned_refmult[i - 1]->Integral());
    vz_corr_binned_refmult_fit.push_back(
        new TF1(sct::MakeString("refmult_vz_bin_fit2_", i).c_str(),
                vz_function.c_str(), vz_tail_min, vz_tail_max));
    vz_corr_binned_refmult_fit[i - 1]->SetParameters(
        vz_func_p0_default, vz_func_p1_default, vz_func_p2_default);
    vz_corr_binned_refmult[i - 1]->Fit(vz_corr_binned_refmult_fit[i - 1], "ME",
                                       "", vz_tail_min, vz_tail_max);
    corr_vz_fit_h_param->SetBinContent(
        i, vz_corr_binned_refmult_fit[i - 1]->GetParameter(2));
    corr_vz_fit_h_param->SetBinError(
        i, vz_corr_binned_refmult_fit[i - 1]->GetParError(2));

    vz_corr_binned_refmult[i - 1]->GetXaxis()->SetRangeUser(vz_tail_min,
                                                            vz_tail_max);
    std::string hist_title = sct::MakeString(
        std::setprecision(3), corr_lumi_vz->GetXaxis()->GetBinLowEdge(i),
        " < v_{z} < ", corr_lumi_vz->GetXaxis()->GetBinUpEdge(i));
    gStyle->SetOptTitle(true);
    sct::PrettyPrint1D(vz_corr_binned_refmult[i - 1], hOpts, cOptsNoLeg, "",
                       corrected_tail_fit_dir_path,
                       sct::MakeString("lumi_vz_corrected_tail_fit_h_", i),
                       hist_title, "refMult", "fit paramter h", "");
    gStyle->SetOptTitle(false);
  }

  // fit the corrected h parameter with a constant
  TF1 *lumi_vz_corrected_h_fit =
      new TF1("full_corr_h", "[0]", FLAGS_vzMin, FLAGS_vzMax);
  lumi_vz_corrected_h_fit->SetParameter(0, vz_pol_func_p0_default);
  corr_vz_fit_h_param->Fit(lumi_vz_corrected_h_fit);
  sct::PrettyPrint1D(corr_vz_fit_h_param, hOpts, cOptsNoLeg, "", FLAGS_outDir,
                     "lumi_vz_corrected_fit_h", "", "v_{z} [cm]",
                     "fit paramter h", "");

  // overlay pre and post correction
  TH1D *uncorr_refmult = uncorr_lumi->ProjectionZ();
  TH1D *full_corr_refmult = corr_lumi_vz->ProjectionZ();

  full_corr_refmult->Scale(1.0 / full_corr_refmult->Integral());
  uncorr_refmult->Scale(1.0 / uncorr_refmult->Integral());
  PrintWithRatio(uncorr_refmult, full_corr_refmult, "uncorrected refmult",
                 "corrected refmult", hOpts, cOptsLowerLegLogy, FLAGS_outDir,
                 "corrected_refmult", "", "refmult", "fraction");

  // write uncorrected and corrected refmult to file
  full_corr_refmult->SetName("refmult");
  full_corr_refmult->Write();
  uncorr_refmult->SetName("refmult_uncorrected");
  uncorr_refmult->Write();
  out_file->Close();

  // write parameters to file
  std::ofstream param_file(FLAGS_outDir + "/" + "refmultcorr_params.txt");
  param_file << "lumi parameters: ";
  for (int i = 0; i < luminosity_correction->GetNpar(); ++i) {
    param_file << luminosity_correction->GetParameter(i);
    param_file << (i == luminosity_correction->GetNpar() - 1 ? "\n" : " ");
  }
  param_file << "vz parameters: ";
  for (int i = 0; i < vz_correction->GetNpar(); ++i) {
    param_file << vz_correction->GetParameter(i);
    param_file << (i == vz_correction->GetNpar() - 1 ? "\n" : " ");
  }

  param_file.close();

  gflags::ShutDownCommandLineFlags();
  return 0;
}

bool AcceptEvent(double vz, double vr, double dVz, unsigned refmult,
                 double lumi) {
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
  return true;
}

result_container<TProfile>
ProfileRefmultInBinsOfVz(TH3D *hist, std::string name_prefix, int nsplits) {
  std::vector<std::shared_ptr<TProfile>> lumi_1d_vz;
  std::vector<std::string> lumi_1d_vz_string;
  for (int i = 0; i < nsplits; ++i) {
    int low_bin = i * FLAGS_vzBins / nsplits + 1;
    int high_bin = (i + 1) * FLAGS_vzBins / nsplits;
    hist->GetXaxis()->SetRange(low_bin, high_bin);

    std::string name = name_prefix + "_lumi_avg_refmult_" + std::to_string(i);
    lumi_1d_vz.push_back(std::shared_ptr<TProfile>(
        (TProfile *)((TH2D *)hist->Project3D("ZY"))
            ->ProfileX(sct::MakeString(name, "_profile_", i).c_str())));
    lumi_1d_vz[i]->SetName(name.c_str());

    // create string with vz range
    double dvz = (FLAGS_vzMax - FLAGS_vzMin) / FLAGS_vzBins;
    double vz_low = FLAGS_vzMin + dvz * (low_bin - 1);
    double vz_high = FLAGS_vzMin + dvz * high_bin;
    lumi_1d_vz_string.push_back(
        sct::MakeString(vz_low, " < v_{z} < ", vz_high));
  }
  return {lumi_1d_vz, lumi_1d_vz_string};
}

result_container<TH1D>
ProjectRefMultInBinsOfVz(TH3D *hist, std::string name_prefix, int nsplits) {
  std::vector<std::shared_ptr<TH1D>> refmult_1d_vz;
  std::vector<std::string> refmult_1d_vz_string;
  for (int i = 0; i < nsplits; ++i) {
    int low_bin = i * FLAGS_vzBins / nsplits + 1;
    int high_bin = (i + 1) * FLAGS_vzBins / nsplits;

    std::string name = name_prefix + "_refmult_vz_bin_" + std::to_string(i);
    refmult_1d_vz.push_back(std::shared_ptr<TH1D>(
        hist->ProjectionZ(sct::MakeString(name, "_projection_", i).c_str(),
                          low_bin, high_bin, 0, -1)));
    refmult_1d_vz[i]->SetName(name.c_str());

    // create string with vz range
    double dvz = (FLAGS_vzMax - FLAGS_vzMin) / FLAGS_vzBins;
    double vz_low = FLAGS_vzMin + dvz * (low_bin - 1);
    double vz_high = FLAGS_vzMin + dvz * high_bin;
    refmult_1d_vz_string.push_back(
        sct::MakeString(vz_low, " < v_{z} < ", vz_high));
  }
  return {refmult_1d_vz, refmult_1d_vz_string};
}

result_container<TH1D>
ProjectRefMultInBinsOfZDC(TH3D *hist, std::string name_prefix, int nsplits) {
  std::vector<std::shared_ptr<TH1D>> refmult_1d_lumi;
  std::vector<std::string> refmult_1d_lumi_string;
  for (int i = 0; i < nsplits; ++i) {
    int low_bin = i * FLAGS_lumiBins / nsplits + 1;
    int high_bin = (i + 1) * FLAGS_lumiBins / nsplits;

    std::string name = name_prefix + "_refmult_lumi_bin_" + std::to_string(i);
    refmult_1d_lumi.push_back(std::shared_ptr<TH1D>(
        hist->ProjectionZ(sct::MakeString(name, "_projection_", i).c_str(), 0,
                          -1, low_bin, high_bin)));
    refmult_1d_lumi[i]->SetName(name.c_str());

    // create string with vz range
    double dlumi = (FLAGS_lumiMax - FLAGS_lumiMin) / FLAGS_lumiBins;
    double lumi_low = FLAGS_lumiMin + dlumi * (low_bin - 1);
    double lumi_high = FLAGS_lumiMin + dlumi * high_bin;
    refmult_1d_lumi_string.push_back(
        sct::MakeString(lumi_low, " < zdcX[kHz] < ", lumi_high));
  }
  return {refmult_1d_lumi, refmult_1d_lumi_string};
}

void RefmultQA(TH3D *h, std::string name, sct::histogramOpts hOpts_hist,
               sct::canvasOpts cOpts_hist, sct::histogramOpts hOpts_prof,
               sct::canvasOpts cOpts_prof, int bins, bool fit_flat) {
  // print the uncorrected average refmult in bins of vz & bins of zdc

  result_container<TH1D> refmult_lumi_bins =
      ProjectRefMultInBinsOfZDC(h, name, bins);

  for (auto h : refmult_lumi_bins.first)
    h->Scale(1.0 / h->Integral());
  Overlay1D(refmult_lumi_bins.first, refmult_lumi_bins.second, hOpts_hist,
            cOpts_hist, FLAGS_outDir, name + "_ref_lumi_bin", "", "refmult",
            "fraction", "");

  result_container<TH1D> refmult_vz_bins =
      ProjectRefMultInBinsOfVz(h, name, bins);
  for (auto h : refmult_vz_bins.first)
    h->Scale(1.0 / h->Integral());
  Overlay1D(refmult_vz_bins.first, refmult_vz_bins.second, hOpts_hist,
            cOpts_hist, FLAGS_outDir, name + "_ref_vz_bin", "", "refmult",
            "fraction", "");

  // print out average refmult as a function of luminosity in bins of vz
  result_container<TProfile> lumi_1d_vz =
      ProfileRefmultInBinsOfVz(h, name, bins);
  for (int i = 0; i < lumi_1d_vz.first.size(); ++i) {
    std::string fit_function = "";
    TF1 *fit = nullptr;
    if (fit_flat) {
      fit = new TF1(
          sct::MakeString(lumi_1d_vz.first[i]->GetName(), "_fit").c_str(),
          "[0]", 0, 100);
      fit->SetParameter(0, 160);
    } else {
      fit = new TF1(
          sct::MakeString(lumi_1d_vz.first[i]->GetName(), "_fit").c_str(),
          "[0]+[1]*x", 0, 100);
      fit->SetParameter(0, 160);
      fit->SetParameter(1, -0.01);
    }

    fit->SetLineColor(hOpts_prof.colors[i % hOpts_prof.colors.size()]);
    TFitResultPtr res = lumi_1d_vz.first[i]->Fit(fit, "MES");
    double chi2 = res.Get()->Chi2();
    unsigned ndf = res.Get()->Ndf();
    lumi_1d_vz.second[i] +=
        sct::MakeString(" chi2 = ", std::setprecision(4), chi2, "/", ndf);
  }
  lumi_1d_vz.first[0]->GetYaxis()->SetRangeUser(120, 220);
  gStyle->SetOptFit(false);
  cOpts_prof.leg_left_bound = 0.4;
  Overlay1D(lumi_1d_vz.first, lumi_1d_vz.second, hOpts_prof, cOpts_prof,
            FLAGS_outDir, name + "_avg_ref_vs_lumi_vz", "", "ZDC Rate [kHz]",
            "<refmult>", "");
  gStyle->SetOptFit(true);
}
