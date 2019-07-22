/* Given a set of multiplicity model parameters, will estimate
 * the errors on the glauber weights by comparing the results
 * using differing initial parameters to the glauber model -
 * increasing or decreasing the NN cross section, varying the
 * nucleus radius, etc. Generates systematic errors.
 *
 * The multiplicity model parameters should be taken from the best
 * fit reported from fit_glauber_to_data. The refmult distribution
 * should be the same distribution used in fit_glauber_to_data,
 * and the glauber directory should contain the nominal as well as
 * variations of the glauber model, from run_glauber.
 *
 */

#include "sct/lib/enumerations.h"
#include "sct/lib/flags.h"
#include "sct/lib/logging.h"
#include "sct/lib/string/string_utils.h"
#include "sct/systematics/glauber_systematics.h"
#include "sct/utils/histogram_info.h"

#include <set>
#include <string>

#include "boost/filesystem.hpp"
#include "boost/range/iterator_range.hpp"

// program settings
SCT_DEFINE_string(glauberDir, "glauber",
                  "directory containing results from run_glauber");
SCT_DEFINE_string(outDir, "tmp", "Path to directory to store output and QA");
SCT_DEFINE_string(outputFile, "glauber_systematics.root", "output file name");

SCT_DEFINE_string(species, "au197", "nucleus species");
SCT_DEFINE_bool(deformation, false, "turn on nucleus deformation");
SCT_DEFINE_int(energy, 200, "collision energy");
SCT_DEFINE_bool(doNppVariations, true, "turn on Npp variations (+-5%)");

// centrality definition - 5% bins, 0-80%. Can set a hard-coded default to be
// used in the code - look for std::vector<double> default_centrality
SCT_DEFINE_string(centDef, "",
                  "centrality bin lower refmult cut, 5% bins from 0-80%, "
                  "separated by commas");

SCT_DEFINE_double(npp, 2.0,
                  "multiplicity model parameter Npp (taken from fit)");
SCT_DEFINE_double(k, 2.0, "multiplicity model parameter K (taken from fit)");
SCT_DEFINE_double(x, 0.35, "multiplicity model parameter X (taken from fit)");
SCT_DEFINE_double(
    ppEff, 0.98, "multiplicity model parameter pp Efficiency (taken from fit)");
SCT_DEFINE_double(
    aaEff, 0.84, "multiplicity model parameter AA Efficiency (taken from fit)");
SCT_DEFINE_double(aaMult, 540,
                  "multiplicity model parameter AA Average 0-5% multiplicity "
                  "(taken from fit)");
SCT_DEFINE_double(
    trigEff, 1.0,
    "multiplicity model parameter trigger efficiency (taken from fit)");
SCT_DEFINE_bool(constEff, false,
                "multiplicity model parameter set efficiency to a constant "
                "(taken from fit)");

// function to check if a file should be used by doing some simple substring
// matching
bool AcceptFile(std::string file, sct::GlauberSpecies species,
                sct::CollisionEnergy energy) {
  // check to make sure it is a glauber root file
  if (file.find("glauber") == 0 && file.find(".root") == file.size() - 5) {
    // check if it has the proper deformation
    if ((FLAGS_deformation == true &&
         file.find("deformed_true") == std::string::npos) ||
        (FLAGS_deformation == false &&
         file.find("deformed_false") == std::string::npos))
      return false;

    // now check to make sure it is for the proper nucleus
    if (file.find(sct::speciesString[species]) != std::string::npos &&
        file.find(std::to_string(static_cast<int>(energy))) !=
            std::string::npos) {
      return true;
    }
  }
  return false;
}

int main(int argc, char* argv[]) {
  // set help message
  std::string usage = "Runs systematic evaluation of centrality definition.";
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
    LOG(ERROR) << "Requested unknown nucleus: " << FLAGS_species << ": exiting";
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
      LOG(ERROR) << "requested unknown energy: " << FLAGS_energy << ": exiting";
      return 1;
  }

  // list of the glauber modifications
  std::set<sct::GlauberMod> modifiers{
      sct::GlauberMod::Nominal,   sct::GlauberMod::Large,
      sct::GlauberMod::Small,     sct::GlauberMod::LargeXSec,
      sct::GlauberMod::SmallXSec, sct::GlauberMod::Gauss};
  std::set<std::string> mod_string{"nominal",   "large",     "small",
                                   "largexsec", "smallxsec", "gauss"};

  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path out_dir(FLAGS_outDir);
  boost::filesystem::create_directories(out_dir);

  // and create filename for output root file
  out_dir /= FLAGS_outputFile;
  std::string out_filename = out_dir.string();

  // now find all the glauber output files in the specified directory
  boost::filesystem::path glauber_path(FLAGS_glauberDir);
  std::set<std::string> glauber_input;
  if (boost::filesystem::is_directory(glauber_path)) {
    for (auto& entry : boost::make_iterator_range(
             boost::filesystem::directory_iterator(glauber_path), {})) {
      std::string file = entry.path().filename().string();
      std::transform(file.begin(), file.end(), file.begin(), ::tolower);

      if (AcceptFile(file, species, energy))
        glauber_input.insert(entry.path().string());
    }
  }

  // check the size
  if (glauber_input.size() != modifiers.size()) {
    LOG(ERROR)
        << "Incorrect number of input glauber files found in input directiory: "
        << FLAGS_glauberDir;
    LOG(ERROR) << "Found(" << glauber_input.size()
               << ") modified glauber runs: " << glauber_input;
    LOG(ERROR) << "expected files for(" << modifiers.size()
               << ") modifiers: " << mod_string;
    LOG(ERROR) << "exiting";
    return 1;
  }

  LOG(INFO) << "Found (" << glauber_input.size()
            << ") systematic modifications: " << glauber_input;

  // load centrality definition
  std::vector<double> default_centrality{460, 390, 329, 275, 228, 187, 152, 122,
                                         96,  75,  57,  42,  31,  22,  15,  10};
  std::vector<double> centrality;

  if (FLAGS_centDef.empty())
    centrality = default_centrality;
  else
    centrality = sct::ParseArgStringToVec<double>(FLAGS_centDef);

  int centrality_bins =
      sct::HistogramInfo::instance().bins(sct::GlauberObservable::Centrality);
  if (centrality.size() != centrality_bins) {
    LOG(ERROR) << "centrality definition not valid: requires "
               << centrality_bins << " bins, but got " << centrality.size();
    return 1;
  }

  // create the analysis class and intialize it with the correct mutliplicity
  // model and input files
  sct::GlauberSystematics systematics(out_filename);
  systematics.setMultiplicityModel(FLAGS_npp, FLAGS_k, FLAGS_x, FLAGS_ppEff,
                                   FLAGS_aaEff, FLAGS_aaMult, FLAGS_trigEff,
                                   FLAGS_constEff);
  systematics.setCentralityBins(centrality);
  systematics.doNppVariations(FLAGS_doNppVariations);
  systematics.load(glauber_input);

  // run the analysis and write to disk
  systematics.run();
  systematics.write();

  // now, we use the results from the analysis to generate the systematics
  // on the centrality definition
  TFile systematics_file(out_filename.c_str(), "READ");

  // now calculate results one at a time, by looping over each observable
  std::set<sct::GlauberObservable> observable_list{
      sct::GlauberObservable::Npart,        sct::GlauberObservable::Ncoll,
      sct::GlauberObservable::Multiplicity, sct::GlauberObservable::B,
      sct::GlauberObservable::RPArea,       sct::GlauberObservable::PPArea,
      sct::GlauberObservable::PP2Ecc,       sct::GlauberObservable::PP3Ecc,
      sct::GlauberObservable::PP4Ecc};

  for (auto& observable : observable_list) {
  }

  gflags::ShutDownCommandLineFlags();
  return 0;
}