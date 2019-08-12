/* Scans through a 3D grid of parameters (Npp, k, x) and builds up  a refmult
 * distribution by sampling the npart vs ncoll distribution from the glauber
 * results N times for every set of parameters. Compares the generated refmult
 * distribution to a data refmult distribution using a chi2 similarity measure.
 * Reports the best fit parameters.
 *
 */

#include "sct/centrality/centrality.h"
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
SCT_DEFINE_string(outFile, "centrality_results",
                  "output file name (no extension)");
SCT_DEFINE_string(inputFile, "glauber.root",
                  "path to root file containing the refmult distributions");
SCT_DEFINE_string(glauberHistName, "glauber",
                  "name of glauber refmult distribution");
SCT_DEFINE_string(dataHistName, "refmult",
                  "name of reference multiplicty histogram");

// if you want to create a reweighted refmult distribution, you can pass in a
// tree with refmult info,
SCT_DEFINE_string(refmultTreeFile, "", "file containing refmult tree");
SCT_DEFINE_string(refmultTreeName, "refMultTree", "name of refmult tree");
SCT_DEFINE_string(refmultBranchName, "refMult", "name of refmult branch");
SCT_DEFINE_string(preGlauberCorrFile, "",
                  "name of file containing lumi & vz correction parameters");

SCT_DEFINE_string(refmultBranch, "refMult",
                  "name of reference multiplicty branch");
SCT_DEFINE_int(refmultMin, 0, "minimum refmult");
SCT_DEFINE_int(refmultMax, 800, "maximum refmult");

SCT_DEFINE_string(luminosityBranch, "lumi",
                  "name of branch containing luminosity information (zdc rate, "
                  "bbc rate, etc)");
SCT_DEFINE_double(lumiMin, 0, "minimum luminosity [kHz]");
SCT_DEFINE_double(lumiMax, 100, "maximum luminosity [kHz]");
SCT_DEFINE_double(
    lumiNorm, 0.0,
    "normalization point - lumi correction = f(zdcNorm)/f(zdc) [kHz]");

SCT_DEFINE_string(vzBranch, "vz", "name of vertex z branch");
SCT_DEFINE_double(vzMin, -30.0, "minimum Vz [cm]");
SCT_DEFINE_double(vzMax, 30.0, "maximum Vz [cm]");
SCT_DEFINE_double(vzNorm, 0.0,
                  "normalization point - vz correction = f(vzNorm)/f(vz) [cm]");

SCT_DEFINE_string(vrBranch, "vr", "name of vertex r branch");
SCT_DEFINE_double(vrMax, 3.0, "maximum Vr [cm]");

SCT_DEFINE_string(dVzBranch, "vzvpdvz", "name of vz - vpd vz branch");
SCT_DEFINE_double(dVzMax, 3.0, "maximum dVz [cm]");

// used for histogram normalization when remaking the scaled glauber
// distribution
SCT_DEFINE_bool(useStGlauberNorm, true,
                "use StGlauber Normalization instead of integral norm");
SCT_DEFINE_int(minMult, 100, "minimum multiplicity for normalization");
SCT_DEFINE_int(fitCutoffLow, 0, "minimum value for centrality reweighting fit");
SCT_DEFINE_int(fitCutoffHigh, 400, "maximum value for centrality reweighting");

double LumiScaling(double zdcX, double zdc_norm_point,
                   std::vector<double> &pars) {
  if (pars.size() != 2) {
    LOG(ERROR) << "expected 2 parameters for luminosity scaling, returning 1";
    return 1.0;
  }
  double zdc_scaling = pars[0] + pars[1] * zdcX;
  double zdc_norm = pars[0] + pars[1] * zdc_norm_point;
  return zdc_norm / zdc_scaling;
}

double VzScaling(double vz, double vz_norm_point, std::vector<double> &pars) {
  if (pars.size() != 7) {
    LOG(ERROR) << "expected 7 parameters for luminosity scaling, returning 1";
    return 1.0;
  }
  double vz_scaling = 0.0;
  double vz_norm = 0.0;
  for (int i = 0; i < pars.size(); ++i) {
    vz_scaling += pars[i] * pow(vz, i);
    vz_norm += pars[i] * pow(vz_norm_point, i);
  }
  if ((vz_norm / vz_scaling) <= 0)
    return 1.0;
  return vz_norm / vz_scaling;
}

double GlauberScaling(double refmult, double cutoff, std::vector<double> pars) {
  if (refmult > cutoff)
    return 1;
  double weight = pars[0] + pars[1] / (pars[2] * refmult + pars[3]) +
                  pars[4] * (pars[2] * refmult + pars[3]);
  weight += pars[5] / pow(pars[2] * refmult + pars[3], 2) +
            pars[6] * pow(pars[2] * refmult + pars[3], 2);
  return 1.0 / weight;
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

int main(int argc, char *argv[]) {
  // shut ROOT up :)
  gErrorIgnoreLevel = kWarning;

  // set help message
  std::string usage =
      "Calculates centrality from a given glauber refmult distribution";
  sct::SetUsageMessage(usage);

  sct::InitLogging(&argc, argv);
  sct::ParseCommandLineFlags(&argc, argv);

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_outDir);
  boost::filesystem::create_directories(dir);

  // load the input files data refmult histogram & glauber npart x ncoll
  TFile *glauber_file = new TFile(FLAGS_inputFile.c_str(), "READ");

  if (!glauber_file->IsOpen()) {
    LOG(ERROR) << "Glauber input file could not be opened: " << FLAGS_inputFile
               << " not found or corrupt";
    return 1;
  }

  // load data refmult histogram & glauber refmult
  TH1D *glauber = (TH1D *)glauber_file->Get(FLAGS_glauberHistName.c_str());
  TH1D *refmult = (TH1D *)glauber_file->Get(FLAGS_dataHistName.c_str());

  if (glauber == nullptr) {
    LOG(ERROR) << "Glauber reference multiplicity histogram: "
               << FLAGS_glauberHistName
               << " not found in file: " << FLAGS_inputFile;
    return 1;
  }
  if (refmult == nullptr) {
    LOG(ERROR) << "Reference multiplicity histogram: " << FLAGS_dataHistName
               << " not found in file: " << FLAGS_inputFile;
    return 1;
  }

  // create output file
  std::string output_name = FLAGS_outDir + "/" + FLAGS_outFile + ".root";
  TFile out(output_name.c_str(), "RECREATE");

  // create our fitting model for histogram normalization
  sct::NBDFit fitter;
  fitter.minimumMultiplicityCut(FLAGS_minMult);
  fitter.useStGlauberNorm(FLAGS_useStGlauberNorm);

  // now calculate centrality definition for the best fit
  sct::Centrality cent;
  cent.setDataRefmult(refmult);
  cent.setSimuRefmult(glauber);
  std::vector<unsigned> cent_bounds = cent.centralityBins(sct::XSecMod::None);
  std::vector<unsigned> cent_bounds_p5 =
      cent.centralityBins(sct::XSecMod::Plus5);
  std::vector<unsigned> cent_bounds_m5 =
      cent.centralityBins(sct::XSecMod::Minus5);
  LOG(INFO) << "finished calculating centrality";
  // get weights
  auto weights = cent.weights(FLAGS_fitCutoffLow, FLAGS_fitCutoffHigh);

  // write the ratio to file
  weights.second->SetName("ratio_fit");
  weights.second->Write();

  // if we have a refmult tree file, we can create a reweighted refmult
  // distribution and take the ratio with the glauber distribution
  if (!FLAGS_refmultTreeFile.empty()) {
    LOG(INFO) << "will be calculating corrected ratio";
    // first, read in vz/luminosity corrections if we use them
    std::vector<double> vz_pars;
    std::vector<double> lumi_pars;
    if (!FLAGS_preGlauberCorrFile.empty()) {
      LOG(INFO) << "using luminosity & vz corrections";
      std::ifstream glauber_in;
      glauber_in.open(FLAGS_preGlauberCorrFile);
      if (glauber_in.is_open()) {
        std::string tmp_string;
        while (std::getline(glauber_in, tmp_string)) {
          auto result = sct::ParseStrToVec<double>(tmp_string);
          if (result.size() == 2)
            lumi_pars = result;
          else if (result.size() == 7)
            vz_pars = result;
        }
      }
    }

    // setup TTreeReader
    LOG(INFO) << "reading in refmult tree";
    TFile refmult_tree_file(FLAGS_refmultTreeFile.c_str(), "read");
    TTreeReader reader(FLAGS_refmultTreeName.c_str(), &refmult_tree_file);
    TTreeReaderValue<unsigned> ref(reader, FLAGS_refmultBranch.c_str());
    TTreeReaderValue<double> vz(reader, FLAGS_vzBranch.c_str());
    TTreeReaderValue<double> zdcx(reader, FLAGS_luminosityBranch.c_str());
    TTreeReaderValue<double> vr(reader, FLAGS_vrBranch.c_str());
    TTreeReaderValue<double> dVz(reader, FLAGS_dVzBranch.c_str());

    // create our histogram
    TH1D *refmult_scaled = new TH1D("weighted_refmult", "", 800, 0, 800);
    refmult_scaled->SetDirectory(0);

    // and qa histograms
    TH1D *hVz =
        new TH1D("vz", ";v_{z}", 100, FLAGS_vzMin - 10, FLAGS_vzMax + 10);
    hVz->SetDirectory(0);
    TH1D *hdVz =
        new TH1D("dvz", ";dv_{z}", 100, -FLAGS_dVzMax - 5, FLAGS_dVzMax + 5);
    hdVz->SetDirectory(0);
    TH1D *hVr =
        new TH1D("vr", ";v_{r}", 100, -FLAGS_vrMax - 1, FLAGS_vrMax + 1);
    hVr->SetDirectory(0);
    TH1D *hzdcx = new TH1D("zdcx", ";ZDCx [kHz]", 100, FLAGS_lumiMin - 10,
                           FLAGS_lumiMax + 10);
    hzdcx->SetDirectory(0);

    while (reader.Next()) {
      double refmultcorr = *ref;

      if (!AcceptEvent(*vz, *vr, *dVz, *ref, *zdcx / 1000.0))
        continue;

      hVz->Fill(*vz);
      hdVz->Fill(*dVz);
      hVr->Fill(*vr);
      hzdcx->Fill(*zdcx / 1000.0);

      if (lumi_pars.size() && vz_pars.size()) {
        double lumi_scaling =
            LumiScaling(*zdcx / 1000.0, FLAGS_lumiNorm, lumi_pars);
        double vz_scaling = VzScaling(*vz, FLAGS_vzNorm, vz_pars);
        refmultcorr *= lumi_scaling;
        refmultcorr *= vz_scaling;
      }

      double weight =
          GlauberScaling(refmultcorr, FLAGS_fitCutoffHigh, weights.first);
      refmult_scaled->Fill(refmultcorr, 1.0 / weight);
    }

    TH1D *scaled_ratio = new TH1D(*refmult_scaled);
    scaled_ratio->SetDirectory(0);
    scaled_ratio->SetName("corrected_ratio");
    scaled_ratio->Scale(1.0 / scaled_ratio->Integral(
                                  scaled_ratio->FindBin(FLAGS_minMult), 800));

    // normalize it
    double scale_factor = fitter.norm(scaled_ratio, glauber);
    scaled_ratio->Scale(1.0 / scale_factor);
    scaled_ratio->Divide(glauber);

    // write to file
    LOG(INFO) << "done with corrected spectra";
    refmult_tree_file.Close();
    out.cd();
    scaled_ratio->Write();
    refmult_scaled->Write();

    LOG(INFO) << "writing qa to file";
    hVz->Write();
    hdVz->Write();
    hVr->Write();
    hzdcx->Write();

    // write glauber and refmult
    glauber->Write();
    refmult->Write();
  }

  LOG(INFO) << "writing centrality definition to file";
  std::ofstream cent_file;
  cent_file.open(FLAGS_outDir + "/" + FLAGS_outFile + ".txt",
                 std::ios::out);

  cent_file << "nominal cent: ";
  for (auto i : cent_bounds)
    cent_file << i << " ";
  cent_file << "\n";
  cent_file << "+5% xsec cent: ";
  for (auto i : cent_bounds_p5)
    cent_file << i << " ";
  cent_file << "\n";
  cent_file << "-5% xsec cent: ";
  for (auto i : cent_bounds_m5)
    cent_file << i << " ";
  cent_file << "\n";
  cent_file << "weights: ";
  for (auto par : weights.first)
    cent_file << par << " ";
  cent_file << "\n";

  cent_file.close();

  gflags::ShutDownCommandLineFlags();
  return 0;
}
