// sct/analysis/test_refmultcorr.cc

#include "sct/core/base.hh"
#include "sct/core/flags.hh"
#include "sct/core/logging.hh"
#include "sct/utils/print_helper.hh"
#include "sct/centrality/centrality_def.hh"

#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <boost/filesystem.hpp>

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TStyle.h"

// program settings
SCT_DEFINE_string(outDir, "tmp", "Path to directory to store output");
SCT_DEFINE_string(outFile, "refmultcorr_check", "output file name (no extension)");
SCT_DEFINE_string(dataFile, "refmult.root", "file containing raw refmult tree");
SCT_DEFINE_string(treeName, "refMultTree", "name of refmult tree");
SCT_DEFINE_string(paramFile, "refmult_params.txt",
                  "name of file containing output from pre_glauber_corrections");
SCT_DEFINE_string(centFile, "refmult_weights.txt",
                  "name of file containing output from fit_glauber_to_data.hh");

// if you want to compare the differences between two centrality definitions, you can use
// these
SCT_DEFINE_string(compareParamFile, "", "second centrality definition parameter file for comparison");
SCT_DEFINE_string(compareWeightFile, "", "second centrality definition weight file for comparison");

SCT_DEFINE_string(luminosityBranch, "lumi",
                  "name of branch containing luminosity information (zdc rate, bbc rate, etc)");
SCT_DEFINE_int(lumiBins, 3, "number of bins in luminosity to compare");
SCT_DEFINE_double(lumiMin, 0, "minimum luminosity (zdc rate, in kHz)");
SCT_DEFINE_double(lumiMax, 100, "maximum luminosity (zdc rate, in kHz)");
SCT_DEFINE_double(lumiNorm, 0.0, "zdc normalization point");
SCT_DEFINE_double(compareLumiNorm, 50000.0, "zdc normalization point for comparison definition");

SCT_DEFINE_string(refmultBranch, "refMult", "name of reference multiplicty branch");
SCT_DEFINE_int(refmultMin, 0, "minimum refmult");
SCT_DEFINE_int(refmultMax, 800, "maximum refmult");

SCT_DEFINE_string(vzBranch, "vz", "name of vertex z branch");
SCT_DEFINE_int(vzBins, 5, "number of bins in Vz to compare");
SCT_DEFINE_double(vzMin, -30.0, "minimum vz (in centimeters");
SCT_DEFINE_double(vzMax, 30.0, "maximum vz (in centimeters");
SCT_DEFINE_double(vzNorm, 0.0, "vz normalization point");
SCT_DEFINE_double(compareVzNorm, 0.0, "vz normalization point");

SCT_DEFINE_string(vrBranch, "vr", "name of vertex r branch");
SCT_DEFINE_double(vrMax, 3.0, "maximum Vr");

SCT_DEFINE_string(dVzBranch, "vzvpdvz", "name of vz - vpd vz branch");
SCT_DEFINE_double(dVzMax, 3.0, "maximum dVz");

SCT_DEFINE_string(runIDBranch, "runId", "name of runid branch");
SCT_DEFINE_string(eventIDBranch, "eventId", "name of eventid branch");

void ParseFileForUInt(std::string filename, std::string id_line, std::vector<unsigned>& container) {
  std::ifstream file(filename);
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      if (line.find(id_line) != std::string::npos) {
        if (std::getline(file, line)) {
          std::vector<std::string> strings;
          sct::SplitString(line, strings);
          for (auto entry : strings) {
            container.push_back(atoi(entry.c_str()));
          }
        }
        else {
          LOG(ERROR) << "Expected unsigned int parameters, found EOF, exiting";
          return;
        }
      }
    }
  }
  else {
    LOG(ERROR) << "Could not open file: " << filename << " to find: " << id_line;
  }
}

void ParseFileForDouble(std::string filename, std::string id_line, std::vector<double>& container) {
  std::ifstream file(filename);
  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      if (line.find(id_line) != std::string::npos) {
        if (std::getline(file, line)) {
          std::vector<std::string> strings;
          sct::SplitString(line, strings);
          for (auto entry : strings) {
            container.push_back(atof(entry.c_str()));
          }
        }
        else {
          LOG(ERROR) << "Expected double parameters, found EOF, exiting";
          return;
        }
      }
    }
  }
  else {
    LOG(ERROR) << "Could not open file: " << filename << " to find: " << id_line;
  }
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

int main(int argc, char* argv[]) {
  
  // set help message
  std::string usage = "test refmult corrections and weights accuracy across ";
  usage += "multiple bins of luminosity & Vz";
  sct::SetUsageMessage(usage);
  
  sct::InitLogging(&argc, argv);
  sct::ParseCommandLineFlags(&argc, argv);
  
  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_outDir);
  boost::filesystem::create_directories(dir);
  
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
  
  // read in output from parameter file, and store vz and luminosity parameters
  std::vector<double> lumi_params;
  std::vector<double> vz_params;
  
  ParseFileForDouble(FLAGS_paramFile, "luminosity parameters", lumi_params);
  ParseFileForDouble(FLAGS_paramFile, "vz parameters", vz_params);
  
  if (vz_params.size() == 0 || lumi_params.size() == 0) {
    LOG(ERROR) << "Couldn't load either vz or luminosity correction parameters, exiting";
    return 1;
  }
  LOG(INFO) << "luminosity weights: " << lumi_params;
  LOG(INFO) << "vz weights: " << vz_params;
  
  // now parse the centrality definition file & the weights
  std::vector<unsigned> cent_bin_nominal;
  std::vector<unsigned> cent_bin_p5;
  std::vector<unsigned> cent_bin_m5;
  std::vector<double> weights;
  
  ParseFileForUInt(FLAGS_centFile, "nominal centrality", cent_bin_nominal);
  ParseFileForUInt(FLAGS_centFile, "+5 xsec", cent_bin_p5);
  ParseFileForUInt(FLAGS_centFile, "-5 xsec", cent_bin_m5);
  ParseFileForDouble(FLAGS_centFile, "glauber weighting", weights);
  
  if (cent_bin_nominal.empty() || cent_bin_m5.empty() || cent_bin_p5.empty() ||
      weights.empty()) {
    LOG(ERROR) << "Couldn't load centrality bins or glauber weights correctly, exiting";
    return 1;
  }
  
  LOG(INFO) << "centrality bins: " << cent_bin_nominal;
  LOG(INFO) << "centrality bins (xsec -5%): " << cent_bin_m5;
  LOG(INFO) << "centrality bins (xsec +5%): " << cent_bin_p5;
  LOG(INFO) << "glauber weights: " << weights;
  
  // load our CentralityDef
  sct::CentralityDef refmultCorr;
  refmultCorr.setVzRange(FLAGS_vzMin, FLAGS_vzMax);
  refmultCorr.setZDCRange(FLAGS_lumiMin * 1000.0, FLAGS_lumiMax * 1000.0);
  refmultCorr.setWeightParameters(weights);
  refmultCorr.setVzParameters(vz_params);
  refmultCorr.setZDCParameters(lumi_params);
  refmultCorr.setCentralityBounds16Bin(cent_bin_nominal);
  refmultCorr.setZDCNormalizationPoint(FLAGS_lumiNorm);
  refmultCorr.setVzNormalizationPoint(FLAGS_vzNorm);
  
  // load our refmult tree
  TFile* tree_file = new TFile(FLAGS_dataFile.c_str(), "READ");
  
  if (!tree_file || !tree_file->IsOpen()) {
    LOG(ERROR) << "refmult file: " << FLAGS_dataFile << " could not be opened";
    return 1;
  }
  
  TTreeReader reader(FLAGS_treeName.c_str(), tree_file);
  TTreeReaderValue<unsigned> refmult(reader, FLAGS_refmultBranch.c_str());
  TTreeReaderValue<unsigned> runid(reader, FLAGS_runIDBranch.c_str());
  TTreeReaderValue<unsigned> eventid(reader, FLAGS_eventIDBranch.c_str());
  TTreeReaderValue<double> vz(reader, FLAGS_vzBranch.c_str());
  TTreeReaderValue<double> lumi(reader, FLAGS_luminosityBranch.c_str());
  TTreeReaderValue<double> vr(reader, FLAGS_vrBranch.c_str());
  TTreeReaderValue<double> dVz(reader, FLAGS_dVzBranch.c_str());
  
  
  // create our histograms
  TH1D* hRefmult = new TH1D("refmult", "", 800, 0, 800);
  TH1D* hCorrRef = new TH1D("refmultcorr", "", 800, 0, 800);
  TH1D* hCorrRefWeighted = new TH1D("refmultcorrweighted", "", 800, 0, 800);
  TH1D* centrality = new TH1D("cent", "", 16, 0, 16);
  TH1D* centrality_weighted = new TH1D("weighted_cent", "", 16, 0, 16);
  
  TH3D* binned_cent = new TH3D("cent3d", ";lumi [kHz];V_{z};centBin",
                               FLAGS_lumiBins,FLAGS_lumiMin, FLAGS_lumiMax,
                               FLAGS_vzBins, FLAGS_vzMin, FLAGS_vzMax,
                               16, 0, 16);
  TH3D* binned_ref = new TH3D("ref3d", ";lumi [kHz];V_{z};centBin",
                              FLAGS_lumiBins,FLAGS_lumiMin, FLAGS_lumiMax,
                              FLAGS_vzBins, FLAGS_vzMin, FLAGS_vzMax,
                              (FLAGS_refmultMax-FLAGS_refmultMin),
                              FLAGS_refmultMin, FLAGS_refmultMax);
  
  // needed if we compare to another centrality definition
  std::vector<std::unordered_set<std::pair<int, int>, sct::PairHash>> default_ids(16, std::unordered_set<std::pair<int, int>, sct::PairHash>());

  while(reader.Next()) {
    
    double lumiKhz = *lumi / 1000.0;
    if (!AcceptEvent(*vz, *vr, *dVz, *refmult, lumiKhz))
      continue;
    
    refmultCorr.setEvent(*runid, *refmult, *lumi, *vz);
    double corrected_refmult =refmultCorr.refMultCorr();
    int bin16 = refmultCorr.centrality16();
    double weight = refmultCorr.weight();
    
    if (bin16 < 0)
      continue;
    
    hRefmult->Fill(*refmult);
    hCorrRef->Fill(corrected_refmult);
    hCorrRefWeighted->Fill(corrected_refmult, weight);
    centrality->Fill(bin16);
    centrality_weighted->Fill(bin16, weight);
    
    binned_cent->Fill(lumiKhz, *vz, bin16, weight);
    binned_ref->Fill(lumiKhz, *vz, corrected_refmult, weight);
    
    // now input them into comparison bins
    default_ids[bin16].insert({*runid, *eventid});
  }
  
  // now print the centrality contents weighted and unweighted
  centrality->Scale(0.8 / centrality->Integral());
  centrality_weighted->Scale(0.8 / centrality_weighted->Integral());
  
  sct::histogramOpts hOpts;
  sct::canvasOpts    cOpts;
  
  centrality->GetYaxis()->SetRangeUser(0.0, 0.1);
  sct::Overlay1D(centrality, centrality_weighted, "unweighted 5% bins", "weighted 5% bins",
                 hOpts, cOpts, FLAGS_outDir, "centrality_distribution", "", "centrality bin",
                 "fraction of events", "");
  
  // now look at the centralities as a function of luminosity/vz
  
  // luminosity first
  TH2D* lumi_cent = (TH2D*) binned_cent->Project3D("XZ");
  std::vector<TH1D*> lumi_cent_projections;
  std::vector<std::string> lumi_cent_projection_names;
  for (int i = 1; i <= lumi_cent->GetNbinsY(); ++i) {
    std::string name = "lumi_cent_" + std::to_string(i);
    lumi_cent_projections.push_back(lumi_cent->ProjectionX(name.c_str(), i, i));
    lumi_cent_projections[i-1]->Scale(0.8 / lumi_cent_projections[i-1]->Integral());
    double dlumi = (FLAGS_lumiMax -FLAGS_lumiMin) / FLAGS_lumiBins;
    double bin = i-1;
    
    std::stringstream bin_name;
    bin_name << std::setprecision(2) << FLAGS_lumiMin + dlumi * bin << " < zdc Rate [kHz] < " << FLAGS_lumiMin + dlumi * (bin+1);
    lumi_cent_projection_names.push_back(bin_name.str());
  }
  
  // same procedure for vz
  TH2D* vz_cent = (TH2D*) binned_cent->Project3D("YZ");
  std::vector<TH1D*> vz_cent_projections;
  std::vector<std::string> vz_cent_projection_names;
  for (int i = 1; i <= vz_cent->GetNbinsY(); ++i) {
    std::string name = "vz_cent_" + std::to_string(i);
    vz_cent_projections.push_back(vz_cent->ProjectionX(name.c_str(), i, i));
    vz_cent_projections[i-1]->Scale(0.8 / vz_cent_projections[i-1]->Integral());
    double dvz = (FLAGS_vzMax -FLAGS_vzMin) / FLAGS_vzBins;
    double bin = i-1;
    
    std::stringstream bin_name;
    bin_name << std::setprecision(2) << FLAGS_vzMin + dvz * bin << " < vz [cm] < " << FLAGS_vzMin + dvz * (bin+1);
    vz_cent_projection_names.push_back(bin_name.str());
  }
  
  // print the two
  lumi_cent_projections[0]->GetYaxis()->SetRangeUser(0.0, 0.1);
  vz_cent_projections[0]->GetYaxis()->SetRangeUser(0.0, 0.1);
  Overlay1D(lumi_cent_projections, lumi_cent_projection_names, hOpts, cOpts, FLAGS_outDir, "lumi_cent_bins",
            "", "centrality bin", "fraction", "");
  Overlay1D(vz_cent_projections, vz_cent_projection_names, hOpts, cOpts, FLAGS_outDir, "vz_cent_bins",
            "", "centrality bin", "fraction", "");

  // create our output root file
  std::string out_name = FLAGS_outDir + "/" + FLAGS_outFile + ".root";
  TFile out(out_name.c_str(), "RECREATE");
  
  // if requested, we can compare the similarity of two centrality definitions
  if (!FLAGS_compareParamFile.empty() && !FLAGS_compareWeightFile.empty()) {
    
    // load new centrality definitions and weights
    std::vector<double> alt_lumi_params;
    std::vector<double> alt_vz_params;
    
    ParseFileForDouble(FLAGS_compareParamFile, "luminosity parameters", alt_lumi_params);
    ParseFileForDouble(FLAGS_compareParamFile, "vz parameters", alt_vz_params);
    
    if (alt_vz_params.size() == 0 || alt_lumi_params.size() == 0) {
      LOG(ERROR) << "Couldn't load either alternate vz or luminosity correction parameters, exiting";
      return 1;
    }
    LOG(INFO) << "alternate luminosity weights: " << alt_lumi_params;
    LOG(INFO) << "alternate vz weights: " << alt_vz_params;
    
    // now parse the centrality definition file & the weights
    std::vector<unsigned> alt_cent_bin_nominal;
    std::vector<unsigned> alt_cent_bin_p5;
    std::vector<unsigned> alt_cent_bin_m5;
    std::vector<double> alt_weights;
    
    ParseFileForUInt(FLAGS_compareWeightFile, "nominal centrality", alt_cent_bin_nominal);
    ParseFileForUInt(FLAGS_compareWeightFile, "+5 xsec", alt_cent_bin_p5);
    ParseFileForUInt(FLAGS_compareWeightFile, "-5 xsec", alt_cent_bin_m5);
    ParseFileForDouble(FLAGS_compareWeightFile, "glauber weighting", alt_weights);
    
    if (alt_cent_bin_nominal.empty() || alt_cent_bin_m5.empty() || alt_cent_bin_p5.empty() ||
        alt_weights.empty()) {
      LOG(ERROR) << "Couldn't load alternate centrality bins or glauber weights correctly, exiting";
      return 1;
    }
    
    LOG(INFO) << "alternalte centrality bins: " << alt_cent_bin_nominal;
    LOG(INFO) << "alternate centrality bins (xsec -5%): " << alt_cent_bin_m5;
    LOG(INFO) << "alternate centrality bins (xsec +5%): " << alt_cent_bin_p5;
    LOG(INFO) << "alternate glauber weights: " << alt_weights;
    
    sct::CentralityDef altRefmultCorr;
    altRefmultCorr.setVzRange(FLAGS_vzMin, FLAGS_vzMax);
    altRefmultCorr.setZDCRange(FLAGS_lumiMin * 1000.0, FLAGS_lumiMax * 1000.0);
    altRefmultCorr.setWeightParameters(alt_weights);
    altRefmultCorr.setVzParameters(alt_vz_params);
    altRefmultCorr.setZDCParameters(alt_lumi_params);
    altRefmultCorr.setCentralityBounds16Bin(alt_cent_bin_nominal);
    altRefmultCorr.setZDCNormalizationPoint(FLAGS_compareLumiNorm);
    altRefmultCorr.setVzNormalizationPoint(FLAGS_compareVzNorm);
    
    TH3D* alt_binned_cent = new TH3D("alt_cent3d", ";lumi [kHz];V_{z};centBin",
                                 FLAGS_lumiBins,FLAGS_lumiMin, FLAGS_lumiMax,
                                 FLAGS_vzBins, FLAGS_vzMin, FLAGS_vzMax,
                                 16, 0, 16);
    TH3D* alt_binned_ref = new TH3D("alt_ref3d", ";lumi [kHz];V_{z};centBin",
                                FLAGS_lumiBins,FLAGS_lumiMin, FLAGS_lumiMax,
                                FLAGS_vzBins, FLAGS_vzMin, FLAGS_vzMax,
                                (FLAGS_refmultMax-FLAGS_refmultMin),
                                FLAGS_refmultMin, FLAGS_refmultMax);
    
    // create our histograms
    TH1D* alt_hRefmult = new TH1D("alt_refmult", "", 800, 0, 800);
    TH1D* alt_hCorrRef = new TH1D("alt_refmultcorr", "", 800, 0, 800);
    TH1D* alt_hCorrRefWeighted = new TH1D("alt_refmultcorrweighted", "", 800, 0, 800);
    TH1D* alt_centrality = new TH1D("alt_cent", "", 16, 0, 16);
    TH1D* alt_centrality_weighted = new TH1D("alt_weighted_cent", "", 16, 0, 16);
    TH2D* bin_shift = new TH2D("bin_shift", ";default centrality; alternate centrality", 16, 0, 16, 16, 0, 16);
    // and we need to store run/event ids
    std::vector<std::unordered_set<std::pair<int, int>, sct::PairHash>> compare_ids(16, std::unordered_set<std::pair<int, int>, sct::PairHash>());
    
    reader.Restart();
    
    while(reader.Next()) {
      
      double lumiKhz = *lumi / 1000.0;
      if (!AcceptEvent(*vz, *vr, *dVz, *refmult, lumiKhz))
        continue;
      
      refmultCorr.setEvent(*runid, *refmult, *lumi, *vz);
      altRefmultCorr.setEvent(*runid, *refmult, *lumi, *vz);
      double corrected_refmult = refmultCorr.refMultCorr();
      double alt_corrected_refmult = altRefmultCorr.refMultCorr();
      int bin16 = refmultCorr.centrality16();
      int alt_bin16 = altRefmultCorr.centrality16();
      double weight = refmultCorr.weight();
      double alt_weight = altRefmultCorr.weight();
      
      if (alt_bin16 < 0)
        continue;
      
      bin_shift->Fill(bin16, alt_bin16);
      
      alt_hRefmult->Fill(*refmult);
      alt_hCorrRef->Fill(corrected_refmult);
      alt_hCorrRefWeighted->Fill(corrected_refmult, weight);
      alt_centrality->Fill(bin16);
      alt_centrality_weighted->Fill(bin16, weight);
      
      alt_binned_cent->Fill(lumiKhz, *vz, bin16, weight);
      alt_binned_ref->Fill(lumiKhz, *vz, corrected_refmult, weight);
      
      compare_ids[alt_bin16].insert({*runid, *eventid});
    }
    
   
    // now print the centrality contents weighted and unweighted
    alt_centrality->Scale(0.8 / alt_centrality->Integral());
    alt_centrality_weighted->Scale(0.8 / alt_centrality_weighted->Integral());
    
    centrality->GetYaxis()->SetRangeUser(0.0, 0.1);
    centrality_weighted->GetYaxis()->SetRangeUser(0.0, 0.1);
    alt_centrality->GetYaxis()->SetRangeUser(0.0, 0.1);
    
    Overlay1D(centrality, alt_centrality, "centrality", "comparison centrality", hOpts, cOpts,
              FLAGS_outDir, "centrality_comparison_unweighted", "", "centrality bin", "fraction");
    Overlay1D(centrality_weighted, alt_centrality_weighted, "centrality", "comparison centrality",
              hOpts, cOpts, FLAGS_outDir, "centrality_comparison_weighted", "", "centrality bin",
              "fraction");
    Overlay1D(alt_centrality, alt_centrality_weighted, "unweighted comparison", "weighted comparison",
              hOpts, cOpts, FLAGS_outDir, "centrality_comparison", "", "centrality bin", "fration");
    Print2DSimple(bin_shift, hOpts, cOpts, FLAGS_outDir, "bin_shift", "", "default centrality",
                  "comparison centrality");
    
    // now do comparisons of # of events that have shifted from one bin to another
    TH1D* cent_shift_fraction = new TH1D("centshift", ";centrality;fraction", 16, 0, 16);
    for (int i = 0; i < 16; ++i) {
      LOG(INFO) << "comparing centrality bin: " << i;
      auto& defaults = default_ids[i];
      auto& compares = compare_ids[i];
      
      double nEvents = defaults.size();
      double nEvents_compare = compares.size();
      double match = 0.0;
      
      // matching!!
      for (auto id : defaults) {
        if (compares.find(id) != compares.end())
          match++;
      }
      
      double stable_fraction = 2.0 * match / (nEvents + nEvents_compare);
      cent_shift_fraction->SetBinContent(i, stable_fraction);
      cent_shift_fraction->SetBinError(i, sqrt(1.0/nEvents + 1.0/(nEvents + nEvents_compare)));
      sct::canvasOpts cOptsLowerLeg;
      cOptsLowerLeg.leg_upper_bound = 0.4;
      cOptsLowerLeg.leg_left_bound = 0.2;
      cOptsLowerLeg.leg_lower_bound = 0.2;
      cOptsLowerLeg.leg_right_bound = 0.5;
      cent_shift_fraction->GetYaxis()->SetRangeUser(0.5, 1.0);
      sct::PrettyPrint1D(cent_shift_fraction, hOpts, cOptsLowerLeg, "", FLAGS_outDir, "stable_event_fraction",
                         "", "centrality_bin", "fraction", "fraction of stable events");
    }
    
    // write to root file
    alt_hRefmult->Write();
    alt_hCorrRef->Write();
    alt_hCorrRefWeighted->Write();
    alt_centrality->Write();
    alt_centrality_weighted->Write();
    alt_binned_ref->Write();
    alt_binned_cent->Write();
    cent_shift_fraction->Write();
  }
  
  // write to disk and close
  hRefmult->Write();
  hCorrRef->Write();
  hCorrRefWeighted->Write();
  centrality->Write();
  centrality_weighted->Write();
  binned_ref->Write();
  binned_cent->Write();
  out.Close();
  

  gflags::ShutDownCommandLineFlags();
  return 0;
}
