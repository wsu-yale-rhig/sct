/* Runs the Glauber Monte-Carlo model
 *
 * Can run systematic variations of the glauber model
 * as well - used to produce systematics for centrality
 * definitions.
 *
 */

#include "sct/glauber/glauber_tree.h"
#include "sct/glauber/mc_glauber.h"
#include "sct/lib/enumerations.h"
#include "sct/lib/flags.h"
#include "sct/lib/logging.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/nucleus_info.h"
#include "sct/utils/random.h"

#include <algorithm>
#include <string>
#include <thread>

#include "boost/filesystem.hpp"

#include "TFile.h"
#include "TH2.h"
#include "TROOT.h"
#include "TThread.h"
#include "TTree.h"

// program settings
SCT_DEFINE_string(outDir, "glauber",
                  "Path to directory to store output ROOT files");
SCT_DEFINE_int(events, 1e6, "number of events");
SCT_DEFINE_bool(systematic, true,
                "runs with systematic variations of the glauber parameters");
SCT_DEFINE_string(modification, "nominal",
                  "only applies when systematic = false. Can choose "
                  "[nominal, large, small, largexsec, smallxsec, gaussian]");

// EXPERIMENTAL - currently, ROOT seems to have a problem with
// this, even though I turned on thread safety...
SCT_DEFINE_bool(multithread, false, "allow parallelization of glauber jobs");

// collision settings
SCT_DEFINE_string(species, "au197", "nucleus species");
SCT_DEFINE_bool(deformation, false, "turn on nucleus deformation");
SCT_DEFINE_int(energy, 200, "collision energy");

void PrintSettings() {
  LOG(INFO) << "Running MC Glauber:";
  LOG(INFO) << "nucleus: " << FLAGS_species;
  LOG(INFO) << "collision energy: " << FLAGS_energy << " GeV";
  LOG(INFO) << "deformation: " << (FLAGS_deformation ? "On" : "Off");
  LOG(INFO) << "number of events: " << FLAGS_events;
}

// job to run a single parameter set
void RunGlauber(sct::GlauberSpecies species, sct::CollisionEnergy energy,
                sct::GlauberMod mod, std::string outDir, bool deformation,
                int nEvents) {
  std::string modstring = sct::glauberModToString[mod];

  // create and run the generator
  sct::MCGlauber generator(species, species, energy, mod, deformation,
                           deformation);
  generator.run(nEvents);
  sct::GlauberTree* result = generator.results();

  int binx = sct::NucleusInfo::instance().massNumber(species) * 2 + 1;
  int biny = sct::NucleusInfo::instance().massNumber(species) * 8;

  TH2D* ncollnpart = new TH2D(sct::MakeString("npartncoll_", modstring).c_str(),
                              ";nPart;nColl", binx, 0, binx, biny, 0, biny);
  for (int i = 0; i < result->getEntries(); ++i) {
    result->getEntry(i);
    ncollnpart->Fill(result->nPart(), result->nColl());
  }

  std::string outName = sct::MakeString(
      outDir, "/glauber_", modstring, "_", FLAGS_species, "_", FLAGS_energy,
      "gev_deformed_", (FLAGS_deformation ? "true" : "false"), ".root");
  // create output file
  TFile file(outName.c_str(), "RECREATE");

  // get woods-saxon distributions
  TH2D* WSA = generator.nuclearPDFA();
  WSA->SetName("rcosthetaA");
  WSA->Write();
  TH2D* WSB = generator.nuclearPDFB();
  WSB->SetName("rcosthetaB");
  WSB->Write();
  TH1D* genIP = generator.generatedImpactParameter();
  genIP->SetName("generatedIP");
  genIP->Write();
  TH1D* acceptIP = generator.acceptedImpactParameter();
  acceptIP->SetName("acceptedIP");
  acceptIP->Write();
  result->write();
  ncollnpart->Write();

  file.Close("nodelete");
  delete ncollnpart;
}

int main(int argc, char* argv[]) {
  // otherwise ROOT will commit suicide when we try to launch
  // multiple threads
  ROOT::EnableThreadSafety();

  // set help message
  std::string usage = "Runs an MC glauber simulation for the specified nuclei ";
  usage += "at the specified energy. can perform systematic variations of the ";
  usage += "glauber & nucleus parameters.";
  sct::SetUsageMessage(usage);

  sct::InitLogging(&argc, argv);
  sct::ParseCommandLineFlags(&argc, argv);

  // read nucleus type from command line arguments
  sct::GlauberSpecies species;
  if (FLAGS_species == "Au197" || FLAGS_species == "au197")
    species = sct::GlauberSpecies::Au197;
  else if (FLAGS_species == "Sm154" || FLAGS_species == "sm154")
    species = sct::GlauberSpecies::Sm154;
  else if (FLAGS_species == "U238" || FLAGS_species == "u238")
    species = sct::GlauberSpecies::U238;
  else if (FLAGS_species == "Pb208" || FLAGS_species == "pb208")
    species = sct::GlauberSpecies::Pb208;
  else if (FLAGS_species == "Cu63" || FLAGS_species == "cu63")
    species = sct::GlauberSpecies::Cu63;
  else {
    LOG(ERROR) << "Requested unknown nucleus: exiting";
    return 1;
  }

  // read energy from command line arguments
  sct::CollisionEnergy energy;
  switch (FLAGS_energy) {
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
      LOG(ERROR) << "requested unknown energy: exiting";
      return 1;
  }

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_outDir);
  boost::filesystem::create_directories(dir);

  // list of the glauber modifications
  std::vector<sct::GlauberMod> modifiers{
      sct::GlauberMod::Nominal,   sct::GlauberMod::Large,
      sct::GlauberMod::Small,     sct::GlauberMod::LargeXSec,
      sct::GlauberMod::SmallXSec, sct::GlauberMod::Gauss};

  PrintSettings();

  // if we are not running systematics, we only have one setting to run
  if (FLAGS_systematic == false) {
    std::string modstring = FLAGS_modification;
    std::transform(modstring.begin(), modstring.end(), modstring.begin(),
                   ::tolower);
    sct::GlauberMod mod = sct::stringToGlauberMod[modstring];

    RunGlauber(species, energy, mod, FLAGS_outDir, FLAGS_deformation,
               FLAGS_events);
  } else {
    // if we can multithread we'll run all the systematics at once
    if (FLAGS_multithread) {
      std::vector<std::thread> workers;
      for (auto mod : modifiers) {
        LOG(INFO) << "Running Glauber with modification: "
                  << sct::glauberModToString[mod];
        workers.push_back(std::thread(RunGlauber, species, energy, mod,
                                      FLAGS_outDir, FLAGS_deformation,
                                      FLAGS_events));
      }
      for (int i = 0; i < workers.size(); ++i) {
        workers[i].join();
        LOG(INFO) << "Running Glauber with modification: "
                  << sct::glauberModToString[modifiers[i]] << " done";
      }
    }

    // otherwise, run sequentially
    else {
      for (auto mod : modifiers) {
        LOG(INFO) << "Running Glauber with modification: "
                  << sct::glauberModToString[mod];
        RunGlauber(species, energy, mod, FLAGS_outDir, FLAGS_deformation,
                   FLAGS_events);
        LOG(INFO) << "done";
      }
    }
  }

  // exit out successfully :)
  LOG(INFO) << "Glauber Monte-Carlo successful";
  gflags::ShutDownCommandLineFlags();
  return 0;
}
