/* Creates a 2D plot of an image for a given impact parameter, collision
 * energy and species.
 */

#include "sct/glauber/mc_glauber.h"
#include "sct/glauber/nucleon.h"
#include "sct/glauber/nucleus.h"
#include "sct/lib/assert.h"
#include "sct/lib/enumerations.h"
#include "sct/lib/flags.h"
#include "sct/lib/logging.h"
#include "sct/lib/memory.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/nucleus_info.h"
#include "sct/utils/print_helper.h"
#include "sct/utils/random.h"

#include <algorithm>
#include <string>
#include <vector>
#include <iomanip>

#include "boost/filesystem.hpp"

#include "TCanvas.h"
#include "TEllipse.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TText.h"
#include "TVector2.h"

SCT_DEFINE_string(outDir, "tmp", "Path to directory to store output results");
SCT_DEFINE_int(events, 1, "number of events");
SCT_DEFINE_string(extension, ".png", ".png or .pdf for output files");
SCT_DEFINE_string(modification, "nominal",
                  "Can choose [nominal, large, small, largexsec, smallxsec, "
                  "gaussian]. In general, nominal should be used");
SCT_DEFINE_double(b_min, 0.0, "minimum impact parameter [fm]");
SCT_DEFINE_double(b_max, 20.0, "maximum impact parameter [fm]");
SCT_DEFINE_bool(
    do_centrality, false,
    "generates 1000 events to calculate an impact parameter based centrality "
    "before generating events for centrality estimation");
SCT_DEFINE_bool(single_impact_parameter, false,
                "generate all events with b=b_min");

// collision settings
SCT_DEFINE_string(speciesA, "Au197", "nucleus species for nucleus A");
SCT_DEFINE_string(speciesB, "Au197", "nucleus species for nucleus B");
SCT_DEFINE_bool(deformationA, false,
                "turn on nucleus deformation for nucleus A");
SCT_DEFINE_bool(deformationB, false,
                "turn on nucleus deformation for nucleus B");
SCT_DEFINE_int(energy, 200, "collision energy");
SCT_DEFINE_bool(random_seed, true,
                "if true, uses a random seed for the thread RNGs - otherwise "
                "sets a default value");

sct::GlauberSpecies LookupSpecies(std::string name);
std::string FormattedSpeciesName(std::string name);
sct::CollisionEnergy LookupCollisionEnergy(int e);
sct::GlauberMod LookupGlauberMod(std::string name);
std::vector<double> FindCentBounds(TH1D *impact_parameter);
void PrintEvent(sct::MCGlauber &gen, std::string outFile, std::string nameA,
                std::string nameB, std::string mod, std::string energy,
                double nucleonxsec, std::string cent, TH2D *axes);

int main(int argc, char *argv[]) {

  // set help message
  std::string usage = "Create a visualization of glauber events.";
  sct::SetUsageMessage(usage);

  sct::InitLogging(&argc, argv);
  sct::ParseCommandLineFlags(&argc, argv);

  if (FLAGS_extension != ".png" && FLAGS_extension != ".pdf") {
    LOG(ERROR) << "file extension not valid: should select .png or .pdf";
    return 0;
  }

  // set root printing styles
  gStyle->SetOptStat(false);
  gStyle->SetOptFit(true);
  gStyle->SetOptTitle(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetHatchesSpacing(1.0);
  gStyle->SetHatchesLineWidth(2);

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_outDir);
  boost::filesystem::create_directories(dir);

  // since the flag is a string we parse it to a sct::GlauberSpecies
  sct::GlauberSpecies speciesA = LookupSpecies(FLAGS_speciesA);
  sct::GlauberSpecies speciesB = LookupSpecies(FLAGS_speciesB);

  sct::CollisionEnergy energy = LookupCollisionEnergy(FLAGS_energy);

  sct::GlauberMod mod = LookupGlauberMod(FLAGS_modification);

  // RNG seed
  if (FLAGS_random_seed) {
    std::random_device rng;
    sct::Random::instance().seed(rng());
  } else {
    sct::Random::instance().seed(sct::Counter::instance().counter());
  }

  // create the generator
  sct::MCGlauber generator(speciesA, speciesB, energy, mod, FLAGS_deformationA,
                           FLAGS_deformationB);

  std::vector<double> cent_bounds;
  if (FLAGS_do_centrality) {
    LOG(INFO) << "Centrality estimation requested: generating impact parameter "
                 "distribution";
    generator.run(5000);
    cent_bounds = FindCentBounds(generator.acceptedImpactParameter());
  }
  std::reverse(cent_bounds.begin(), cent_bounds.end());

  // set impact parameter
  generator.setImpactParameterRange(FLAGS_b_min, FLAGS_b_max);
  if (FLAGS_single_impact_parameter) {
    generator.setImpactParameterRange(FLAGS_b_min, FLAGS_b_min + 0.0001);
  }

  // create histogram to print the coordinate axes
  sct::histogramOpts hopts;
  hopts.title_offset_y = 0.7;
  hopts.center_title_x = true;
  double x_limit = 18;
  double y_limit = 18;
  TH2D *h_axes = new TH2D("axes", ";x [fm]; y [fm]", 10, -x_limit, x_limit, 10,
                          -y_limit, y_limit);
  hopts.SetHistogram(h_axes);

  for (int i = 0; i < FLAGS_events; ++i) {
    if (i % 10 == 0)
      LOG(INFO) << "Generating event: " << i;
    sct::EventStatus ret = generator.generate();
    if (ret == sct::EventStatus::Error) {
      LOG(ERROR) << "Error during event generation";
      continue;
    }
    if (ret == sct::EventStatus::Miss) {
      i--;
      continue;
    }
    sct::Nucleus *nucleusA = generator.nucleusA();
    sct::Nucleus *nucleusB = generator.nucleusB();

    boost::filesystem::path image_path(FLAGS_outDir);
    std::string cent_string;
    if (FLAGS_do_centrality) {
      int cent_bin = -1;
      double b = generator.results()->B();
      for (int c = 0; c < cent_bounds.size(); ++c) {
        if (b < cent_bounds[c]) {
          cent_bin = c;
          cent_string = sct::MakeString(sct::Centrality16BinLowerBound[c], "-",
                                        sct::Centrality16BinUpperBound[c], "%");
          break;
        }
        cent_string = "80-100%";
      }
      image_path /=
          sct::MakeString("event_", i, "_centrality_", cent_bin, ".png");
    } else
      image_path /=
          sct::MakeString("event_", i, "_b_", generator.results()->B(), ".png");
    image_path.replace_extension(FLAGS_extension);

    std::string image_filename = image_path.string();
    PrintEvent(generator, image_filename, FLAGS_speciesA, FLAGS_speciesB,
               FLAGS_modification, sct::MakeString(FLAGS_energy),
               generator.lookupXSec(energy), cent_string, h_axes);
  }

  gflags::ShutDownCommandLineFlags();
  return 0;
}

sct::GlauberSpecies LookupSpecies(std::string name) {
  std::transform(name.begin(), name.end(), name.begin(), ::tolower);
  for (auto &sp : sct::speciesString) {
    auto &key = sp.first;
    auto &string_key = sp.second;
    if (string_key.compare(name) == 0)
      return key;
  }
  SCT_THROW("unknown species: ", name, " can not run");
}

std::string FormattedSpeciesName(std::string name) {
  std::size_t found = name.find_first_of("1234567890");
  std::string number = "";
  while (found != std::string::npos) {
    number += name[found];
    name.erase(found, 1);
    found = name.find_first_of("1234567890");
  }
  std::string ret = sct::MakeString("^{", number, "}", name);
  return ret;
}

sct::CollisionEnergy LookupCollisionEnergy(int e) {
  sct::CollisionEnergy energy;
  switch (e) {
  case 2760:
    energy = sct::CollisionEnergy::E2760;
    break;
  case 200:
    energy = sct::CollisionEnergy::E200;
    break;
  case 62:
    energy = sct::CollisionEnergy::E62;
    break;
  case 39:
    energy = sct::CollisionEnergy::E39;
    break;
  case 27:
    energy = sct::CollisionEnergy::E27;
    break;
  case 19:
    energy = sct::CollisionEnergy::E19;
    break;
  case 14:
    energy = sct::CollisionEnergy::E14;
    break;
  case 11:
    energy = sct::CollisionEnergy::E11;
    break;
  case 7:
    energy = sct::CollisionEnergy::E7;
    break;
  default:
    SCT_THROW("Requested collision energy not known to SCT: ", e);
  }
  return energy;
}

sct::GlauberMod LookupGlauberMod(std::string name) {
  std::transform(name.begin(), name.end(), name.begin(), ::tolower);
  for (auto &sp : sct::glauberModToString) {
    auto &key = sp.first;
    auto &string_key = sp.second;
    if (string_key.compare(name) == 0)
      return key;
  }
  SCT_THROW("unknown glauber modification: ", name, " can not run");
}

// helper function to try and avoid rounding errors, when doing
// the fractional calculations in the centrality definitions
template <typename T> T Round(T t, int digits) {
  if (t == 0.0) // otherwise it will return 'nan' due to the log10() of zero
    return 0.0;

  double factor = pow(10.0, digits - ceil(log10(fabs(t))));
  return round(t * factor) / factor;
}

std::vector<double> FindCentBounds(TH1D *h) {
  std::vector<double> bounds;

  // copy our centrality bounds from sct/lib/enumeration.h
  std::vector<double> cent_min = sct::Centrality16BinLowerBound;
  std::vector<double> cent_max = sct::Centrality16BinUpperBound;

  // get the number of cuts
  unsigned n_cent_bins = cent_min.size();
  bounds.resize(n_cent_bins);

  std::vector<double> cent_cuts;
  for (int i = 0; i < n_cent_bins; ++i)
    cent_cuts.push_back(cent_max[i] / 100.0);
  std::reverse(cent_cuts.begin(), cent_cuts.end());

  unsigned cent_bin = 0;
  unsigned nbins_hist = h->GetNbinsX();
  double norm = h->Integral();

  for (int i = nbins_hist; i >= 1; --i) {
    // define the bins to integrate over
    unsigned bin_low = i;
    unsigned bin_high = nbins_hist;

    // get the integral
    double integral = h->Integral(bin_low, bin_high);

    // and the fraction of the total integral
    double fraction = integral / norm;

    // we need to round the fraction to avoid some rounding errors in the
    // floating point arithmetic, which of course involves more errors...
    fraction = Round(fraction, 12);

    for (int bin = cent_bin; bin < n_cent_bins; ++bin) {
      if (1.0 - fraction < cent_cuts[bin]) {
        bounds[cent_bin++] = h->GetXaxis()->GetBinCenter(bin_low);
      }
    }
  }
  return bounds;
}

double GetParticipantPlane(sct::Nucleus *nucleus1, sct::Nucleus *nucleus2,
                           int order = 2) {

  double mQx = 0.0;
  double mQy = 0.0;

  for (auto &nucleus : {nucleus1, nucleus2})
    for (int i = 0; i < nucleus->size(); ++i) {
      sct::Nucleon &n = (*nucleus)[i];
      if (n.nColl() == 0)
        continue;
      mQx += cos(n.phi() * order);
      mQy += sin(n.phi() * order);
    }

  constexpr double pi = 3.14159265359;
  TVector2 mQ(mQx, mQy);
  double psi = mQ.Phi() / order;
  return psi * 180.0 / pi;
}

void PrintEvent(sct::MCGlauber &gen, std::string outFile, std::string nameA,
                std::string nameB, std::string mod, std::string energy,
                double nucleonxsec, std::string cent, TH2D *axes) {

  // setup the canvas margins and coordinates
  sct::canvasOpts copts;
  double margins = 0.13;
  copts.left_margin = margins;
  copts.right_margin = margins;
  copts.upper_margin = margins;
  copts.bottom_margin = margins;
  TCanvas c("c", "c", 720, 720);
  copts.SetMargins(&c);
  c.SetFixedAspectRatio();
  double width = gPad->GetWw() * gPad->GetAbsWNDC();
  double height = gPad->GetWh() * gPad->GetAbsHNDC();
  double xmin = -1;
  double xmax = 1;
  // double ymin = - ((xmax - xmin) * height / width) / 2.0;
  // double ymax = ((xmax - xmin) * height / width) / 2.0;
  double ymin = -1;
  double ymax = 1;
  c.Range(xmin, ymin, xmax, ymax);

  // plot the axes first
  axes->Draw();

  // draw line along impact parameter
  TLine impact_par(-18.0, 0.0, 18.0, 0.0);
  impact_par.Draw();

  // now plot the nuclei as circles
  // first, find b in the canvas coordinates
  double b = gen.results()->B();
  // get nucleus radius
  double radA = sct::NucleusInfo::instance().radius(LookupSpecies(nameA));
  double radB = sct::NucleusInfo::instance().radius(LookupSpecies(nameB));

  TEllipse n1(-b / 2.0, 0.0, radA, radA);
  n1.SetFillStyle(0);
  n1.SetLineWidth(3.0);
  TEllipse n2(b / 2.0, 0.0, radB, radB);
  n2.SetFillStyle(0);
  n2.SetLineWidth(3.0);
  n1.Draw();
  n2.Draw();

  // print relevant event info
  std::string part_string =
      sct::MakeString("N_{part}: ", gen.results()->nPart());
  std::string bin_string =
      sct::MakeString("N_{binary}: ", gen.results()->nColl());
  std::string b_string = sct::MakeString("b: ", std::setprecision(3), gen.results()->B(), " fm");
  std::string system_name = sct::MakeString(FormattedSpeciesName(nameA), " + ", FormattedSpeciesName(nameB));
  std::string system_energy =
      sct::MakeString("#sqrt{s_{NN}} = ", energy, " GeV");
  std::string nucleon_xsec =
      sct::MakeString("NN xsec = ", nucleonxsec, " fm^{2}");
  // centrality if present
  std::string cent_string = "centrality: " + cent;

  // set TPaveTexts
  double x_lower = -16;
  double x_upper = -8;
  double y_lower = 11;
  double y_upper = 17;
  if (!cent.empty()) {
    x_lower = -15.5;
    x_upper = -7.5;
    y_lower = 10.5;
  }
  TPaveText pt_left(x_lower, y_lower, x_upper, y_upper);
  pt_left.SetBorderSize(0);
  pt_left.SetFillStyle(0);
  pt_left.AddText(part_string.c_str())->SetTextSize(0.03);
  pt_left.AddText(bin_string.c_str())->SetTextSize(0.03);
  pt_left.AddText(b_string.c_str())->SetTextSize(0.03);
  if (!cent.empty()) {
    pt_left.AddText(cent_string.c_str())->SetTextSize(0.03);
  }

  TPaveText pt_right(8, 11, 16, 17);
  pt_right.SetBorderSize(0);
  pt_right.SetFillStyle(0);
  pt_right.AddText(system_name.c_str())->SetTextSize(0.03);
  pt_right.AddText(system_energy.c_str())->SetTextSize(0.03);
  pt_right.AddText(nucleon_xsec.c_str())->SetTextSize(0.03);

  pt_left.Draw();
  pt_right.Draw();

  // radii for nucleus from NN xsec
  constexpr double pi = 3.14159265359;
  double nucleon_radius = sqrt(nucleonxsec / pi) / 2.0;

  // draw nucleons in first nucleus
  for (int i = 0; i < gen.nucleusA()->size(); ++i) {
    sct::Nucleon &n = (*gen.nucleusA())[i];
    double x = n.x();
    double y = n.y();
    TEllipse t(x, y, nucleon_radius, nucleon_radius);
    t.SetFillStyle(0);
    t.SetLineColor(kRed);
    if (n.nColl() == 0)
      t.SetLineColorAlpha(kRed, 0.2);
    t.DrawClone();
  }

  // draw nucleons in second nucleus
  for (int i = 0; i < gen.nucleusB()->size(); ++i) {
    sct::Nucleon &n = (*gen.nucleusB())[i];
    double x = n.x();
    double y = n.y();
    TEllipse t(x, y, nucleon_radius, nucleon_radius);
    t.SetFillStyle(0);
    t.SetLineColor(kBlue);
    if (n.nColl() == 0)
      t.SetLineColorAlpha(kBlue, 0.2);
    t.DrawClone();
  }

  // draw a "legend" for spectators and participants
  double legend_radius = nucleon_radius * 1.5;
  TEllipse n1_spec(-15, -15, legend_radius, legend_radius);
  n1_spec.SetFillStyle(0);
  n1_spec.SetLineColorAlpha(kRed, 0.2);
  TEllipse n1_part(-15, -13, legend_radius, legend_radius);
  n1_part.SetFillStyle(0);
  n1_part.SetLineColorAlpha(kRed, 1.0);

  TEllipse n2_spec(-12, -15, legend_radius, legend_radius);
  n2_spec.SetFillStyle(0);
  n2_spec.SetLineColorAlpha(kBlue, 0.2);
  TEllipse n2_part(-12, -13, legend_radius, legend_radius);
  n2_part.SetFillStyle(0);
  n2_part.SetLineColorAlpha(kBlue, 1.0);

  TEllipse n_leg(-13.5, -11, legend_radius, legend_radius);
  n_leg.SetFillStyle(0);
  n_leg.SetLineColorAlpha(kBlack, 1.0);

  n1_spec.Draw();
  n1_part.Draw();
  n2_spec.Draw();
  n2_part.Draw();
  n_leg.Draw();

  TPaveText leg1(-9.5, -14, -4, -12);
  leg1.SetBorderSize(0);
  leg1.SetFillStyle(0);
  leg1.AddText("participants")->SetTextSize(0.025);
  leg1.Draw();

  TPaveText leg2(-9.5, -16, -4, -14);
  leg2.SetBorderSize(0);
  leg2.SetFillStyle(0);
  leg2.AddText("spectators")->SetTextSize(0.025);
  leg2.Draw();

  TPaveText leg3(-9.5, -12, -4, -10);
  leg3.SetBorderSize(0);
  leg3.SetFillStyle(0);
  leg3.AddText("nuclei")->SetTextSize(0.025);
  leg3.Draw();

  // add participant plane
  double part_plane_2 = GetParticipantPlane(gen.nucleusA(), gen.nucleusB(), 2);

  TVector2 pp2_norm(15.0, 0.0);
  pp2_norm = pp2_norm.Rotate(part_plane_2 * pi / 180.0);
  pp2_norm = pp2_norm.Rotate(pi / 2.0);

  // create TLine for psi2
  double pp2_x = pp2_norm.X();
  double pp2_y = pp2_norm.Y();

  TLine pp2_line(-pp2_x, -pp2_y, pp2_x, pp2_y);

  pp2_line.SetLineColor(kMagenta + 1);
  pp2_line.SetLineWidth(2);
  pp2_line.Draw();

  // create TLines for psi
  double part_plane_3 = GetParticipantPlane(gen.nucleusA(), gen.nucleusB(), 3);
  TVector2 pp3_norm(10.0, 0.0);
  pp3_norm = pp3_norm.Rotate(part_plane_3 * pi / 180.0);
  pp3_norm = pp3_norm.Rotate(pi / 3.0);

  double pp3_x1 = pp3_norm.X();
  double pp3_y1 = pp3_norm.Y();

  pp3_norm = pp3_norm.Rotate(2.0 * pi / 3.0);
  double pp3_x2 = pp3_norm.X();
  double pp3_y2 = pp3_norm.Y();

  pp3_norm = pp3_norm.Rotate(2.0 * pi / 3.0);
  double pp3_x3 = pp3_norm.X();
  double pp3_y3 = pp3_norm.Y();

  TLine pp3_line1(0.0, 0.0, pp3_x1, pp3_y1);
  TLine pp3_line2(0.0, 0.0, pp3_x2, pp3_y2);
  TLine pp3_line3(0.0, 0.0, pp3_x3, pp3_y3);

  pp3_line1.SetLineColor(kCyan + 1);
  pp3_line1.SetLineWidth(2);
  pp3_line1.Draw();
  pp3_line2.SetLineColor(kCyan + 1);
  pp3_line2.SetLineWidth(2);
  pp3_line2.Draw();
  pp3_line3.SetLineColor(kCyan + 1);
  pp3_line3.SetLineWidth(2);
  pp3_line3.Draw();

  // final legend
  TLine pp2_leg(3.0, -12, 5.0, -12);
  TLine pp3_leg(3.0, -15, 5.0, -15);

  pp2_leg.SetLineColor(kMagenta + 1);
  pp2_leg.SetLineWidth(2);
  pp2_leg.Draw();

  pp3_leg.SetLineColor(kCyan + 1);
  pp3_leg.SetLineWidth(2);
  pp3_leg.Draw();

  TPaveText pp2_leg_text(8, -13, 14, -11);
  pp2_leg_text.SetBorderSize(0);
  pp2_leg_text.SetFillStyle(0);
  pp2_leg_text.AddText("#Psi_{2} participant plane")->SetTextSize(0.025);
  pp2_leg_text.Draw();

  TPaveText pp3_leg_text(8, -16, 14, -14);
  pp3_leg_text.SetBorderSize(0);
  pp3_leg_text.SetFillStyle(0);
  pp3_leg_text.AddText("#Psi_{3} participant plane")->SetTextSize(0.025);
  pp3_leg_text.Draw();

  // save to disk
  c.SaveAs(outFile.c_str());

  return;
}