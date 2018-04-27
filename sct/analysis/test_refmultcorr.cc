// sct/analysis/test_refmultcorr.cc

#include "sct/core/base.hh"
#include "sct/core/flags.hh"
#include "sct/core/logging.hh"
#include "sct/utils/print_helper.hh"
#include "sct/centrality/centrality_def.hh"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/filesystem.hpp>

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"

// program settings
SCT_DEFINE_string(outDir, "tmp", "Path to directory to store output");
SCT_DEFINE_string(outFile, "refmultcorr_check", "output file name (no extension)");
SCT_DEFINE_string(dataFile, "refmult.root", "file containing raw refmult tree");
SCT_DEFINE_string(treeName, "refMultTree", "name of refmult tree");
SCT_DEFINE_string(paramFile, "refmult_params.txt",
                  "name of file containing output from pre_glauber_corrections");
SCT_DEFINE_string(centFile, "refmult_weights.txt",
                  "name of file containing output from fit_glauber_to_data.hh");

SCT_DEFINE_string(luminosityBranch, "lumi",
                  "name of branch containing luminosity information (zdc rate, bbc rate, etc)");
SCT_DEFINE_int(lumiBins, 3, "number of bins in luminosity to compare");
SCT_DEFINE_double(lumiMin, 0, "minimum luminosity (zdc rate, in kHz)");
SCT_DEFINE_double(lumiMax, 100, "maximum luminosity (zdc rate, in kHz)");
SCT_DEFINE_double(lumiNorm, 0.0, "zdc normalization point");

SCT_DEFINE_string(refmultBranch, "refMult", "name of reference multiplicty branch");
SCT_DEFINE_int(refmultMin, 0, "minimum refmult");
SCT_DEFINE_int(refmultMax, 800, "maximum refmult");

SCT_DEFINE_string(vzBranch, "vz", "name of vertex z branch");
SCT_DEFINE_int(vzBins, 5, "number of bins in Vz to compare");
SCT_DEFINE_double(vzMin, -30.0, "minimum vz (in centimeters");
SCT_DEFINE_double(vzMax, 30.0, "maximum vz (in centimeters");
SCT_DEFINE_double(vzNorm, 0.0, "vz normalization point");

SCT_DEFINE_string(vrBranch, "vr", "name of vertex r branch");
SCT_DEFINE_double(vrMax, 3.0, "maximum Vr");

SCT_DEFINE_string(dVzBranch, "vzvpdvz", "name of vz - vpd vz branch");
SCT_DEFINE_double(dVzMax, 3.0, "maximum dVz");

SCT_DEFINE_string(runIDBranch, "runId", "name of runid branch");

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
  TH3D* binned_ref = new TH3D("cent3d", ";lumi [kHz];V_{z};centBin",
                              FLAGS_lumiBins,FLAGS_lumiMin, FLAGS_lumiMax,
                              FLAGS_vzBins, FLAGS_vzMin, FLAGS_vzMax,
                              (FLAGS_refmultMax-FLAGS_refmultMin),
                              FLAGS_refmultMin, FLAGS_refmultMax);
  
  while(reader.Next()) {
    
    double lumiKhz = *lumi / 1000.0;
    if (!AcceptEvent(*vz, *vr, *dVz, *refmult, lumiKhz))
      continue;
    
    refmultCorr.setEvent(*runid, *refmult, *lumi, *vz);
    double corrected_refmult =refmultCorr.refMultCorr();
    unsigned bin16 = refmultCorr.centrality16();
    double weight = refmultCorr.weight();
    
    hRefmult->Fill(*refmult);
    hCorrRef->Fill(corrected_refmult);
    hCorrRefWeighted->Fill(corrected_refmult, weight);
    centrality->Fill(bin16);
    centrality_weighted->Fill(bin16, weight);
    
    binned_cent->Fill(lumiKhz, *vz, bin16);
    binned_ref->Fill(lumiKhz, *vz, corrected_refmult);
  }
  
  // write our output
  std::string out_name = FLAGS_outDir + "/" + FLAGS_outFile + ".root";
  TFile out(out_name.c_str(), "RECREATE");
  hRefmult->Write();
  hCorrRef->Write();
  hCorrRefWeighted->Write();
  centrality->Write();
  centrality_weighted->Write();
  binned_ref->Write();
  binned_cent->Write();
  out.Close();
  
  // print the 
  
  
  gflags::ShutDownCommandLineFlags();
  return 0;
}
