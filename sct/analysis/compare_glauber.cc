/* Used to compare the new Glauber MC results to the
 * StFastGlauberMC results for QA of the new library.
 * (it will be labeled in the code where either the
 * old result should be considered incorrect, or where
 * methods have changed so that the results should not
 * be comparable.)
 *
 * Uses TTreeReader instead of sct::tree because the
 * output trees of the two libraries are slightly different
 */

#include "sct/lib/flags.h"
#include "sct/lib/logging.h"
#include "sct/lib/enumerations.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/print_helper.h"

#include <boost/filesystem.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"

// program settings
SCT_DEFINE_string(outDir, "tmp",
                      "Path to directory to store output files");
SCT_DEFINE_string(outPrefix, "", "output file name prefix");
SCT_DEFINE_string(newInput, "glauber_new.root", "new glauber ROOT file");
SCT_DEFINE_string(newTreeName, "event", "new glauber TTree name");
SCT_DEFINE_string(oldInput, "glauber_old.root", "old glauber ROOT file");
SCT_DEFINE_string(oldTreeName, "tree", "old glauber TTree name");

int main(int argc, char* argv[]) {
  
  constexpr unsigned npart_idx = static_cast<unsigned>(sct::GlauberWeight::NPart);
  constexpr unsigned ncoll_idx = static_cast<unsigned>(sct::GlauberWeight::NColl);
  constexpr unsigned spec_idx = static_cast<unsigned>(sct::GlauberWeight::Spectators);
  
  // set help message
  std::string usage = "Compares new MCGlauber results to StFastGlauberMC";
  sct::SetUsageMessage(usage);
  
  sct::InitLogging(&argc, argv);
  sct::ParseCommandLineFlags(&argc, argv);
  
  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_outDir);
  boost::filesystem::create_directories(dir);
 
  // load in the new glauber tree
  TFile *new_file = new TFile(FLAGS_newInput.c_str(), "READ");
  TTreeReader new_reader(FLAGS_newTreeName.c_str(), new_file);
  TTreeReaderValue<double> new_b(new_reader, "b");
  TTreeReaderValue<unsigned> new_npart(new_reader, "npart");
  TTreeReaderValue<unsigned> new_ncoll(new_reader, "ncoll");
  TTreeReaderValue<unsigned> new_nspec(new_reader, "nspec");
  TTreeReaderArray<double> new_sumx(new_reader, "sumx");
  TTreeReaderArray<double> new_sumy(new_reader, "sumy");
  TTreeReaderArray<double> new_sumx2(new_reader, "sumx2");
  TTreeReaderArray<double> new_sumy2(new_reader, "sumy2");
  TTreeReaderArray<double> new_sumxy(new_reader, "sumxy");
  TTreeReaderArray<double> new_eccrp2(new_reader, "eccrp2");
  TTreeReaderArray<double> new_eccpp2(new_reader, "eccpp2");
  TTreeReaderArray<double> new_eccpp3(new_reader, "eccpp3");
  TTreeReaderArray<double> new_eccpp4(new_reader, "eccpp4");
  TTreeReaderArray<double> new_pp2(new_reader, "pp2");
  TTreeReaderArray<double> new_pp3(new_reader, "pp3");
  TTreeReaderArray<double> new_pp4(new_reader, "pp4");
  
  TH2D* h_gen_rcostheta = (TH2D*) new_file->Get("rcosthetaA");
  TH1D* h_new_ws = h_gen_rcostheta->ProjectionX();
  
  // load in the old glauber tree
  TFile *old_file = new TFile(FLAGS_oldInput.c_str(), "READ");
  TTreeReader old_reader(FLAGS_oldTreeName.c_str(), old_file);
  TTreeReaderValue<double> old_b(old_reader, "b");
  TTreeReaderValue<unsigned> old_npart(old_reader, "npart");
  TTreeReaderValue<unsigned> old_ncoll(old_reader, "ncoll");
  TTreeReaderArray<double> old_sumx(old_reader, "sumx");
  TTreeReaderArray<double> old_sumy(old_reader, "sumy");
  TTreeReaderArray<double> old_sumx2(old_reader, "sumx2");
  TTreeReaderArray<double> old_sumy2(old_reader, "sumy2");
  TTreeReaderArray<double> old_sumxy(old_reader, "sumxy");
  TTreeReaderArray<double> old_eccrp2(old_reader, "eccrp2");
  TTreeReaderArray<double> old_eccpp2(old_reader, "eccpp2");
  TTreeReaderArray<double> old_eccpp3(old_reader, "eccpp3");
  TTreeReaderArray<double> old_eccpp4(old_reader, "eccpp4");
  TTreeReaderArray<double> old_pp2(old_reader, "pp2");
  TTreeReaderArray<double> old_pp3(old_reader, "pp3");
  TTreeReaderArray<double> old_pp4(old_reader, "pp4");
  
  // woods saxon
  TH1D* h_old_ws = (TH1D*) old_file->Get("hWoodsSaxon_0");
  
  // event statistics
  TH1D* h_new_b = new TH1D("new_b", ";impact parameter[fm]", 100, 0, 20);
  TH1D* h_old_b = new TH1D("old_b", ";impact parameter[fm]", 100, 0, 20);
  TH2D* h_new_npnc = new TH2D("new_npartncoll", ";npart;ncoll", 100, 0, 400, 100, 0, 1500);
  TH1D* h_new_nspec = new TH1D("new_nspec", ";nspec", 100, 0, 400);
  TH2D* h_new_bnpart = new TH2D("new_bnpart", ";b[fm];npart", 100, 0, 20, 100, 0, 400);
  TH2D* h_new_bncoll = new TH2D("new_bncoll", ";b[fm];ncoll", 100, 0, 20, 100, 0, 1500);
  TH2D* h_old_npnc = new TH2D("old_npartncoll", ";npart;ncoll", 100, 0, 400, 100, 0, 1500);
  TH2D* h_old_bnpart = new TH2D("old_bnpart", ";b[fm];npart", 100, 0, 20, 100, 0, 400);
  TH2D* h_old_bncoll = new TH2D("old_bncoll", ";b[fm];ncoll", 100, 0, 20, 100, 0, 1500);
  
  // participant/spectator kinematics
  TH1D* h_new_part_avg_x = new TH1D("new_part_avgx", "", 100, -10, 10);
  TH1D* h_new_part_avg_y = new TH1D("new_part_avgy", "", 100, -10, 10);
  TH1D* h_new_part_avg_x2 = new TH1D("new_part_avgx2", "", 100, 0, 10);
  TH1D* h_new_part_avg_y2 = new TH1D("new_part_avgy2", "", 100, 0, 10);
  TH1D* h_new_part_avg_xy = new TH1D("new_part_avgxy", "", 100, -10, 10);
  TH2D* h_new_part_avg_xy_2D = new TH2D("new_part_avgxy2D", "", 100, -10, 10, 100, -10, 10);
  
  TH1D* h_new_coll_avg_x = new TH1D("new_coll_avgx", "", 100, -10, 10);
  TH1D* h_new_coll_avg_y = new TH1D("new_coll_avgy", "", 100, -10, 10);
  TH1D* h_new_coll_avg_x2 = new TH1D("new_coll_avgx2", "", 100, 0, 10);
  TH1D* h_new_coll_avg_y2 = new TH1D("new_coll_avgy2", "", 100, 0, 10);
  TH1D* h_new_coll_avg_xy = new TH1D("new_coll_avgxy", "", 100, -10, 10);
  TH2D* h_new_coll_avg_xy_2D = new TH2D("new_coll_avgxy2D", "", 100, -10, 10, 100, -10, 10);
  
  TH1D* h_new_spec_avg_x = new TH1D("new_spec_avgx", "", 100, -10, 10);
  TH1D* h_new_spec_avg_y = new TH1D("new_spec_avgy", "", 100, -10, 10);
  TH1D* h_new_spec_avg_x2 = new TH1D("new_spec_avgx2", "", 100, 0, 10);
  TH1D* h_new_spec_avg_y2 = new TH1D("new_spec_avgy2", "", 100, 0, 10);
  TH1D* h_new_spec_avg_xy = new TH1D("new_spec_avgxy", "", 100, -10, 10);
  TH2D* h_new_spec_avg_xy_2D = new TH2D("new_spec_avgxy2D", "", 100, -10, 10, 100, -10, 10);
  
  TH1D* h_old_part_avg_x = new TH1D("old_part_avgx", "", 100, -10, 10);
  TH1D* h_old_part_avg_y = new TH1D("old_part_avgy", "", 100, -10, 10);
  TH1D* h_old_part_avg_x2 = new TH1D("old_part_avgx2", "", 100, 0, 10);
  TH1D* h_old_part_avg_y2 = new TH1D("old_part_avgy2", "", 100, 0, 10);
  TH1D* h_old_part_avg_xy = new TH1D("old_part_avgxy", "", 100, -10, 10);
  TH2D* h_old_part_avg_xy_2D = new TH2D("old_part_avgxy2D", "", 100, -10, 10, 100, -10, 10);
  
  TH1D* h_old_coll_avg_x = new TH1D("old_coll_avgx", "", 100, -10, 10);
  TH1D* h_old_coll_avg_y = new TH1D("old_coll_avgy", "", 100, -10, 10);
  TH1D* h_old_coll_avg_x2 = new TH1D("old_coll_avgx2", "", 100, 0, 10);
  TH1D* h_old_coll_avg_y2 = new TH1D("old_coll_avgy2", "", 100, 0, 10);
  TH1D* h_old_coll_avg_xy = new TH1D("old_coll_avgxy", "", 100, -10, 10);
  TH2D* h_old_coll_avg_xy_2D = new TH2D("old_coll_avgxy2D", "", 100, -10, 10, 100, -10, 10);
  
  TH1D* h_old_spec_avg_x = new TH1D("old_spec_avgx", "", 100, -10, 10);
  TH1D* h_old_spec_avg_y = new TH1D("old_spec_avgy", "", 100, -10, 10);
  TH1D* h_old_spec_avg_x2 = new TH1D("old_spec_avgx2", "", 100, 0, 10);
  TH1D* h_old_spec_avg_y2 = new TH1D("old_spec_avgy2", "", 100, 0, 10);
  TH1D* h_old_spec_avg_xy = new TH1D("old_spec_avgxy", "", 100, -10, 10);
  TH2D* h_old_spec_avg_xy_2D = new TH2D("old_spec_avgxy2D", "", 100, -10, 10, 100, -10, 10);
  
  // reaction plane eccentricity
  TH1D* h_new_part_rpecc2 = new TH1D("new_part_rpecc2", ";eccentricity", 100, 0, 1);
  TH1D* h_new_coll_rpecc2 = new TH1D("new_coll_rpecc2", ";eccentricity", 100, 0, 1);
  TH1D* h_new_spec_rpecc2 = new TH1D("new_spec_rpecc2", ";eccentricity", 100, 0, 1);
  TH1D* h_old_part_rpecc2 = new TH1D("old_part_rpecc2", ";eccentricity", 100, 0, 1);
  TH1D* h_old_coll_rpecc2 = new TH1D("old_coll_rpecc2", ";eccentricity", 100, 0, 1);
  TH1D* h_old_spec_rpecc2 = new TH1D("old_spec_rpecc2", ";eccentricity", 100, 0, 1);
  
  // participant plane
  TH1D* h_new_part_pp2 = new TH1D("new_part_pp2", ";pp2", 100, -sct::pi, sct::pi);
  TH1D* h_new_part_pp3 = new TH1D("new_part_pp3", ";pp3", 100, -sct::pi, sct::pi);
  TH1D* h_new_part_pp4 = new TH1D("new_part_pp4", ";pp4", 100, -sct::pi, sct::pi);
  TH1D* h_new_coll_pp2 = new TH1D("new_coll_pp2", ";pp2", 100, -sct::pi, sct::pi);
  TH1D* h_new_coll_pp3 = new TH1D("new_coll_pp3", ";pp3", 100, -sct::pi, sct::pi);
  TH1D* h_new_coll_pp4 = new TH1D("new_coll_pp4", ";pp4", 100, -sct::pi, sct::pi);
  TH1D* h_new_spec_pp2 = new TH1D("new_spec_pp2", ";pp2", 100, -sct::pi, sct::pi);
  TH1D* h_new_spec_pp3 = new TH1D("new_spec_pp3", ";pp3", 100, -sct::pi, sct::pi);
  TH1D* h_new_spec_pp4 = new TH1D("new_spec_pp4", ";pp4", 100, -sct::pi, sct::pi);
  
  TH1D* h_old_part_pp2 = new TH1D("old_part_pp2", ";pp2", 100, -sct::pi, sct::pi);
  TH1D* h_old_part_pp3 = new TH1D("old_part_pp3", ";pp3", 100, -sct::pi, sct::pi);
  TH1D* h_old_part_pp4 = new TH1D("old_part_pp4", ";pp4", 100, -sct::pi, sct::pi);
  TH1D* h_old_coll_pp2 = new TH1D("old_coll_pp2", ";pp2", 100, -sct::pi, sct::pi);
  TH1D* h_old_coll_pp3 = new TH1D("old_coll_pp3", ";pp3", 100, -sct::pi, sct::pi);
  TH1D* h_old_coll_pp4 = new TH1D("old_coll_pp4", ";pp4", 100, -sct::pi, sct::pi);
  TH1D* h_old_spec_pp2 = new TH1D("old_spec_pp2", ";pp2", 100, -sct::pi, sct::pi);
  TH1D* h_old_spec_pp3 = new TH1D("old_spec_pp3", ";pp3", 100, -sct::pi, sct::pi);
  TH1D* h_old_spec_pp4 = new TH1D("old_spec_pp4", ";pp4", 100, -sct::pi, sct::pi);
  
  // eccentricity as a function of npart
  TProfile* h_new_part_ecc2 = new TProfile("new_npart_ecc2", "", 100, 0, 400);
  TProfile* h_new_part_ecc3 = new TProfile("new_npart_ecc3", "", 100, 0, 400);
  TProfile* h_new_part_ecc4 = new TProfile("new_npart_ecc4", "", 100, 0, 400);
  
  // loop over old data first
  while(new_reader.Next()) {
    h_new_b->Fill(*new_b);
    h_new_npnc->Fill(*new_npart, *new_ncoll);
    h_new_nspec->Fill(*new_nspec);
    h_new_bnpart->Fill(*new_b, *new_npart);
    h_new_bncoll->Fill(*new_b, *new_ncoll);
    
    h_new_part_avg_x->Fill(new_sumx[npart_idx]);
    h_new_part_avg_y->Fill(new_sumy[npart_idx]);
    h_new_part_avg_x2->Fill(sqrt(new_sumx2[npart_idx]));
    h_new_part_avg_y2->Fill(sqrt(new_sumy2[npart_idx]));
    h_new_part_avg_xy->Fill(new_sumxy[npart_idx]);
    h_new_part_avg_xy_2D->Fill(new_sumx[npart_idx], new_sumy[npart_idx]);
    h_new_coll_avg_x->Fill(new_sumx[ncoll_idx]);
    h_new_coll_avg_y->Fill(new_sumy[ncoll_idx]);
    h_new_coll_avg_x2->Fill(sqrt(new_sumx2[ncoll_idx]));
    h_new_coll_avg_y2->Fill(sqrt(new_sumy2[ncoll_idx]));
    h_new_coll_avg_xy->Fill(new_sumxy[ncoll_idx]);
    h_new_coll_avg_xy_2D->Fill(new_sumx[ncoll_idx], new_sumy[ncoll_idx]);
    h_new_spec_avg_x->Fill(new_sumx[spec_idx]);
    h_new_spec_avg_y->Fill(new_sumy[spec_idx]);
    h_new_spec_avg_x2->Fill(sqrt(new_sumx2[spec_idx]));
    h_new_spec_avg_y2->Fill(sqrt(new_sumy2[spec_idx]));
    h_new_spec_avg_xy->Fill(new_sumxy[spec_idx]);
    h_new_spec_avg_xy_2D->Fill(new_sumx[spec_idx], new_sumy[spec_idx]);
    
    if (*new_npart <= 2)
      continue;
    
    h_new_part_rpecc2->Fill(new_eccrp2[npart_idx]);
    h_new_coll_rpecc2->Fill(new_eccrp2[ncoll_idx]);
    h_new_spec_rpecc2->Fill(new_eccrp2[spec_idx]);
    
    h_new_part_pp2->Fill(new_pp2[npart_idx]);
    h_new_part_pp3->Fill(new_pp3[npart_idx]);
    h_new_part_pp4->Fill(new_pp4[npart_idx]);
    h_new_coll_pp2->Fill(new_pp2[ncoll_idx]);
    h_new_coll_pp3->Fill(new_pp3[ncoll_idx]);
    h_new_coll_pp4->Fill(new_pp4[ncoll_idx]);
    h_new_spec_pp2->Fill(new_pp2[spec_idx]);
    h_new_spec_pp3->Fill(new_pp3[spec_idx]);
    h_new_spec_pp4->Fill(new_pp4[spec_idx]);
    
    if (new_eccpp2[npart_idx] == 1.0 || new_eccpp3[npart_idx] == 1.0)
      continue;
    
    h_new_part_ecc2->Fill(*new_npart, new_eccpp2[npart_idx]);
    h_new_part_ecc3->Fill(*new_npart, new_eccpp3[npart_idx]);
    h_new_part_ecc4->Fill(*new_npart, new_eccpp4[npart_idx]);
  }
  
  // loop over old data
  while (old_reader.Next()) {
    h_old_b->Fill(*old_b);
    h_old_npnc->Fill(*old_npart, *old_ncoll);
    h_old_bnpart->Fill(*old_b, *old_npart);
    h_old_bncoll->Fill(*old_b, *old_ncoll);
    
    h_old_part_avg_x->Fill(old_sumx[0]);
    h_old_part_avg_y->Fill(old_sumy[0]);
    h_old_part_avg_x2->Fill(sqrt(old_sumx2[0]));
    h_old_part_avg_y2->Fill(sqrt(old_sumy2[0]));
    h_old_part_avg_xy->Fill(old_sumxy[0]);
    h_old_part_avg_xy_2D->Fill(old_sumx[0], old_sumy[0]);
    h_old_coll_avg_x->Fill(old_sumx[1]);
    h_old_coll_avg_y->Fill(old_sumy[1]);
    h_old_coll_avg_x2->Fill(sqrt(old_sumx2[1]));
    h_old_coll_avg_y2->Fill(sqrt(old_sumy2[1]));
    h_old_coll_avg_xy->Fill(old_sumxy[1]);
    h_old_coll_avg_xy_2D->Fill(old_sumx[1], old_sumy[1]);
    h_old_spec_avg_x->Fill(old_sumx[3]);
    h_old_spec_avg_y->Fill(old_sumy[3]);
    h_old_spec_avg_x2->Fill(sqrt(old_sumx2[3]));
    h_old_spec_avg_y2->Fill(sqrt(old_sumy2[3]));
    h_old_spec_avg_xy->Fill(old_sumxy[3]);
    h_old_spec_avg_xy_2D->Fill(old_sumx[3], old_sumy[3]);
    
    if (*old_npart <= 2)
      continue;
    
    h_old_part_rpecc2->Fill(old_eccrp2[0]);
    h_old_coll_rpecc2->Fill(old_eccrp2[1]);
    h_old_spec_rpecc2->Fill(old_eccrp2[3]);
    
    h_old_part_pp2->Fill(old_pp2[0]);
    h_old_part_pp3->Fill(old_pp3[0]);
    h_old_part_pp4->Fill(old_pp4[0]);
    h_old_coll_pp2->Fill(old_pp2[1]);
    h_old_coll_pp3->Fill(old_pp3[1]);
    h_old_coll_pp4->Fill(old_pp4[1]);
    h_old_spec_pp2->Fill(old_pp2[3]);
    h_old_spec_pp3->Fill(old_pp3[3]);
    h_old_spec_pp4->Fill(old_pp4[3]);
  }
  
  h_new_ws->Scale(1.0 / h_new_ws->Integral());
  h_old_ws->Scale(1.0 / h_old_ws->Integral());
  
  h_new_b->Scale(1.0 / h_new_b->Integral());
  h_old_b->Scale(1.0 / h_old_b->Integral());
  h_new_npnc->Scale(1.0 / h_new_npnc->Integral());
  h_old_npnc->Scale(1.0 / h_old_npnc->Integral());
  h_new_nspec->Scale(1.0 / h_new_nspec->Integral());
  h_new_bnpart->Scale(1.0 / h_new_bnpart->Integral());
  h_new_bncoll->Scale(1.0 / h_new_bncoll->Integral());
  h_old_bnpart->Scale(1.0 / h_old_bnpart->Integral());
  h_old_bncoll->Scale(1.0 / h_old_bncoll->Integral());
  
  h_new_part_avg_x->Scale(1.0 / h_new_part_avg_x->Integral());
  h_new_part_avg_y->Scale(1.0 / h_new_part_avg_y->Integral());
  h_new_part_avg_x2->Scale(1.0 / h_new_part_avg_x2->Integral());
  h_new_part_avg_y2->Scale(1.0 / h_new_part_avg_y2->Integral());
  h_new_part_avg_xy->Scale(1.0 / h_new_part_avg_xy->Integral());
  h_new_part_avg_xy_2D->Scale(1.0 / h_new_part_avg_xy_2D->Integral());
  h_new_coll_avg_x->Scale(1.0 / h_new_coll_avg_x->Integral());
  h_new_coll_avg_y->Scale(1.0 / h_new_coll_avg_y->Integral());
  h_new_coll_avg_x2->Scale(1.0 / h_new_coll_avg_x2->Integral());
  h_new_coll_avg_y2->Scale(1.0 / h_new_coll_avg_y2->Integral());
  h_new_coll_avg_xy->Scale(1.0 / h_new_coll_avg_xy->Integral());
  h_new_coll_avg_xy_2D->Scale(1.0 / h_new_coll_avg_xy_2D->Integral());
  h_new_spec_avg_x->Scale(1.0 / h_new_spec_avg_x->Integral());
  h_new_spec_avg_y->Scale(1.0 / h_new_spec_avg_y->Integral());
  h_new_spec_avg_x2->Scale(1.0 / h_new_spec_avg_x2->Integral());
  h_new_spec_avg_y2->Scale(1.0 / h_new_spec_avg_y2->Integral());
  h_new_spec_avg_xy->Scale(1.0 / h_new_spec_avg_xy->Integral());
  h_new_spec_avg_xy_2D->Scale(1.0 / h_new_spec_avg_xy_2D->Integral());

  h_old_part_avg_x->Scale(1.0 / h_old_part_avg_x->Integral());
  h_old_part_avg_y->Scale(1.0 / h_old_part_avg_y->Integral());
  h_old_part_avg_x2->Scale(1.0 / h_old_part_avg_x2->Integral());
  h_old_part_avg_y2->Scale(1.0 / h_old_part_avg_y2->Integral());
  h_old_part_avg_xy->Scale(1.0 / h_old_part_avg_xy->Integral());
  h_old_part_avg_xy_2D->Scale(1.0 / h_old_part_avg_xy_2D->Integral());
  h_old_coll_avg_x->Scale(1.0 / h_old_coll_avg_x->Integral());
  h_old_coll_avg_y->Scale(1.0 / h_old_coll_avg_y->Integral());
  h_old_coll_avg_x2->Scale(1.0 / h_old_coll_avg_x2->Integral());
  h_old_coll_avg_y2->Scale(1.0 / h_old_coll_avg_y2->Integral());
  h_old_coll_avg_xy->Scale(1.0 / h_old_coll_avg_xy->Integral());
  h_old_coll_avg_xy_2D->Scale(1.0 / h_old_coll_avg_xy_2D->Integral());
  h_old_spec_avg_x->Scale(1.0 / h_old_spec_avg_x->Integral());
  h_old_spec_avg_y->Scale(1.0 / h_old_spec_avg_y->Integral());
  h_old_spec_avg_x2->Scale(1.0 / h_old_spec_avg_x2->Integral());
  h_old_spec_avg_y2->Scale(1.0 / h_old_spec_avg_y2->Integral());
  h_old_spec_avg_xy->Scale(1.0 / h_old_spec_avg_xy->Integral());
  h_old_spec_avg_xy_2D->Scale(1.0 / h_old_spec_avg_xy_2D->Integral());
  
  h_new_part_rpecc2->Scale(1.0 / h_new_part_rpecc2->Integral());
  h_new_coll_rpecc2->Scale(1.0 / h_new_coll_rpecc2->Integral());
  h_new_spec_rpecc2->Scale(1.0 / h_new_spec_rpecc2->Integral());
  
  h_old_part_rpecc2->Scale(1.0 / h_old_part_rpecc2->Integral());
  h_old_coll_rpecc2->Scale(1.0 / h_old_coll_rpecc2->Integral());
  h_old_spec_rpecc2->Scale(1.0 / h_old_spec_rpecc2->Integral());
  
  h_new_part_pp2->Scale(1.0 / h_new_part_pp2->Integral());
  h_new_part_pp3->Scale(1.0 / h_new_part_pp3->Integral());
  h_new_part_pp4->Scale(1.0 / h_new_part_pp4->Integral());
  h_new_coll_pp2->Scale(1.0 / h_new_coll_pp2->Integral());
  h_new_coll_pp3->Scale(1.0 / h_new_coll_pp3->Integral());
  h_new_coll_pp4->Scale(1.0 / h_new_coll_pp4->Integral());
  h_new_spec_pp2->Scale(1.0 / h_new_spec_pp2->Integral());
  h_new_spec_pp3->Scale(1.0 / h_new_spec_pp3->Integral());
  h_new_spec_pp4->Scale(1.0 / h_new_spec_pp4->Integral());
  
  h_old_part_pp2->Scale(1.0 / h_old_part_pp2->Integral());
  h_old_part_pp3->Scale(1.0 / h_old_part_pp3->Integral());
  h_old_part_pp4->Scale(1.0 / h_old_part_pp4->Integral());
  h_old_coll_pp2->Scale(1.0 / h_old_coll_pp2->Integral());
  h_old_coll_pp3->Scale(1.0 / h_old_coll_pp3->Integral());
  h_old_coll_pp4->Scale(1.0 / h_old_coll_pp4->Integral());
  h_old_spec_pp2->Scale(1.0 / h_old_spec_pp2->Integral());
  h_old_spec_pp3->Scale(1.0 / h_old_spec_pp3->Integral());
  h_old_spec_pp4->Scale(1.0 / h_old_spec_pp4->Integral());
  
  // print comparisons
  
  // set drawing preferences for histograms and graphs
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(false);
  gStyle->SetLegendBorderSize(0);
  
  // create histogram options
  sct::histogramOpts hOpts;
  hOpts.marker_size = 1.0;
  sct::canvasOpts cOpts;
  sct::canvasOpts cOptsNoLeg;
  cOptsNoLeg.do_legend = false;
  sct::canvasOpts cOptsLogY;
  cOptsLogY.log_y = true;
  
  // woods-saxon distribution
  sct::PrintWithRatio(h_new_ws, h_old_ws, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "woodssaxon"), "",
                          "radius[fm]", "fraction");
  
  // impact parameter
  
  sct::PrintWithRatio(h_new_b, h_old_b, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "b"), "",
                          "impact parameter[fm]", "fraction");
  
  // now get the npart and ncoll distributions
  TH1D* h_new_npart = h_new_npnc->ProjectionX();
  h_new_npart->Scale(1.0 / h_new_npart->Integral());
  TH1D* h_old_npart = h_old_npnc->ProjectionX();
  h_old_npart->Scale(1.0 / h_old_npart->Integral());
  TH1D* h_new_ncoll = h_new_npnc->ProjectionY();
  h_new_ncoll->Scale(1.0 / h_new_ncoll->Integral());
  TH1D* h_old_ncoll = h_old_npnc->ProjectionY();
  h_old_ncoll->Scale(1.0 / h_old_ncoll->Integral());
  
  sct::PrintWithRatio(h_new_npart, h_old_npart, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOptsLogY, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "npart"), "",
                          "participants", "probability");
  sct::PrintWithRatio(h_new_ncoll, h_old_ncoll, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOptsLogY, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "ncoll"), "",
                          "binary collisions", "probability");
  
  // average ncoll and average npart as a function of b
  TH1D* h_new_b_avg_part = h_new_bnpart->ProfileX();
  TH1D* h_old_b_avg_part = h_old_bnpart->ProfileX();
  TH1D* h_new_b_avg_coll = h_new_bncoll->ProfileX();
  TH1D* h_old_b_avg_coll = h_old_bncoll->ProfileX();
  
  sct::PrintWithRatio(h_new_b_avg_part, h_old_b_avg_part, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "avgnpart"), "",
                          "impact parameter", "<participants>");
  sct::PrintWithRatio(h_new_b_avg_coll, h_old_b_avg_coll, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "avgncoll"), "",
                          "impact parameter", "<binary collisions>");
  
  // average x
  sct::PrintWithRatio(h_new_part_avg_x, h_old_part_avg_x, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "partavgx"), "",
                          "<x>", "fraction", "participant <x>");
  sct::PrintWithRatio(h_new_coll_avg_x, h_old_coll_avg_x, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "collavgx"), "",
                          "<x>", "fraction", "binary collision <x>");
  sct::PrintWithRatio(h_new_spec_avg_x, h_old_spec_avg_x, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "specavgx"), "",
                          "<x>", "fraction", "spectator <x>");
  
  // average y
  sct::PrintWithRatio(h_new_part_avg_y, h_old_part_avg_y, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "partavgy"), "",
                          "<y>", "fraction", "participant <y>");
  sct::PrintWithRatio(h_new_coll_avg_y, h_old_coll_avg_y, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "collavgy"), "",
                          "<y>", "fraction", "binary collision <y>");
  sct::PrintWithRatio(h_new_spec_avg_y, h_old_spec_avg_y, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "specavgy"), "",
                          "<y>", "fraction", "spectator <y>");
  
  // average x^2
  sct::PrintWithRatio(h_new_part_avg_x2, h_old_part_avg_x2, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "partavgx2"), "",
                          "<|x|>", "fraction", "participant <|x|>");
  sct::PrintWithRatio(h_new_coll_avg_x2, h_old_coll_avg_x2, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "collavgx2"), "",
                          "<|x|>", "fraction", "binary collision <|x|>");
  sct::PrintWithRatio(h_new_spec_avg_x2, h_old_spec_avg_x2, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "specavgx2"), "",
                          "<|x|>", "fraction", "spectator <|x|>");
  
  // average y^2
  sct::PrintWithRatio(h_new_part_avg_y2, h_old_part_avg_y2, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "partavgy2"), "",
                          "<|y|>", "fraction", "participant <|y|>");
  sct::PrintWithRatio(h_new_coll_avg_y2, h_old_coll_avg_y2, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "collavgy2"), "",
                          "<|y|>", "fraction", "binary collision <|y|>");
  sct::PrintWithRatio(h_new_spec_avg_y2, h_old_spec_avg_y2, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "specavgy2"), "",
                          "<|y|>", "fraction", "spectator <|y|>");
  
  // average xy
  sct::PrintWithRatio(h_new_part_avg_xy, h_old_part_avg_xy, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "partavgxy"), "",
                          "<xy>", "fraction", "participant <xy>");
  sct::PrintWithRatio(h_new_coll_avg_xy, h_old_coll_avg_xy, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "collavgxy"), "",
                          "<xy>", "fraction", "binary collision <xy>");
  sct::PrintWithRatio(h_new_spec_avg_xy, h_old_spec_avg_xy, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "specavgxy"), "",
                          "<xy>", "fraction", "spectator <xy>");
  
  // reaction plane eccentricity
  sct::PrintWithRatio(h_new_part_rpecc2, h_old_part_rpecc2, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "partrpecc2"), "",
                          "#epsilon_{RP}", "fraction", "participant #epsilon_{RP}");
  sct::PrintWithRatio(h_new_coll_rpecc2, h_old_coll_rpecc2, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "collrpecc2"), "",
                          "#epsilon_{RP}", "fraction", "binary collision #epsilon_{RP}");
  sct::PrintWithRatio(h_new_spec_rpecc2, h_old_spec_rpecc2, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "specrpecc2"), "",
                          "#epsilon_{RP}", "fraction", "spectator #epsilon_{RP}");
  
  // participant plane second order
  sct::PrintWithRatio(h_new_part_pp2, h_old_part_pp2, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "partpp2"), "",
                          "#Psi_{pp}^{2}", "fraction", "participant plane");
  sct::PrintWithRatio(h_new_coll_pp2, h_old_coll_pp2, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "collpp2"), "",
                          "#Psi_{pp}^{2}", "fraction", "collision weighted participant plane");
  sct::PrintWithRatio(h_new_spec_pp2, h_old_spec_pp2, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "specpp2"), "",
                          "#Psi_{pp}^{2}", "fraction", "spectator participant plane");
  
  // participant plane third order
  sct::PrintWithRatio(h_new_part_pp3, h_old_part_pp3, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "partpp3"), "",
                          "#Psi_{pp}^{3}", "fraction", "participant plane");
  sct::PrintWithRatio(h_new_coll_pp3, h_old_coll_pp3, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "collpp3"), "",
                          "#Psi_{pp}^{3}", "fraction", "collision weighted participant plane");
  sct::PrintWithRatio(h_new_spec_pp3, h_old_spec_pp3, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "specpp3"), "",
                          "#Psi_{pp}^{3}", "fraction", "spectator participant plane");
  
  // participant plane fourth order
  sct::PrintWithRatio(h_new_part_pp4, h_old_part_pp4, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "partpp4"), "",
                          "#Psi_{pp}^{4}", "fraction", "participant plane");
  sct::PrintWithRatio(h_new_coll_pp4, h_old_coll_pp4, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "collpp4"), "",
                          "#Psi_{pp}^{4}", "fraction", "collision weighted participant plane");
  sct::PrintWithRatio(h_new_spec_pp4, h_old_spec_pp4, "New Glauber library",
                          "StFastGlauberMC", hOpts, cOpts, FLAGS_outDir,
                          sct::MakeString(FLAGS_outPrefix, "specpp4"), "",
                          "#Psi_{pp}^{4}", "fraction", "spectator participant plane");
  
  // eccentricity
  sct::Overlay1D(h_new_part_ecc2, h_new_part_ecc3, "eccentricity", "triangularity",
                     hOpts, cOpts, FLAGS_outDir, "ppecc", "", "participants", "<#epsilon_{n}>");
  
  
  gflags::ShutDownCommandLineFlags();
  return 0;
}
