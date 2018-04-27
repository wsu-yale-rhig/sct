// sct/analysis/run_systematics.cc

/* Given a set of multiplicity model parameters, will estimate
 * the errors on the glauber weights by comparing the results
 * using differing initial parameters to the glauber model -
 * increasing or decreasing the NN cross section, varying the
 * nucleus radius, etc. Reports systematic errors as well as
 * centrality definitions.
 *
 * The multiplicity model parameters should be taken from the best
 * fit reported from fit_glauber_to_data. The refmult distribution
 * should be the same distribution used in fit_glauber_to_data,
 * and the glauber directory should contain the nominal as well as
 * variations of the glauber model, from run_glauber.
 *
 */

#include "sct/core/base.hh"
#include "sct/core/flags.hh"
#include "sct/core/logging.hh"
#include "sct/core/enumerations.hh"
#include "sct/utils/print_helper.hh"
#include "sct/systematics/glauber_systematics.hh"

#include <string>
#include <vector>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>

// program settings
SCT_DEFINE_string(glauberDir, "glauber", "directory containing results from run_glauber");
SCT_DEFINE_string(outDir, "tmp",
                  "Path to directory to store output and QA");
SCT_DEFINE_string(outputPrefix, "corrected_refmult.root", "output file prefix");
SCT_DEFINE_string(refmultFile, "refmult.root", "path to root file containing refmult distribution");
SCT_DEFINE_string(dataRefmultName, "dataRefmult", "data refmult histogram name");
SCT_DEFINE_string(simuRefmultName, "simuRefmult", "simulated refmult histogram name");

SCT_DEFINE_string(species, "au197", "nucleus species");
SCT_DEFINE_bool(deformation, false, "turn on nucleus deformation");
SCT_DEFINE_int(energy, 200, "collision energy");

SCT_DEFINE_double(npp, 2.0, "multiplicity model parameter Npp (taken from fit)");
SCT_DEFINE_double(k, 2.0, "multiplicity model parameter K (taken from fit)");
SCT_DEFINE_double(x, 0.35, "multiplicity model parameter X (taken from fit)");
SCT_DEFINE_double(ppEff, 0.9, "multiplicity model parameter pp Efficiency (taken from fit)");
SCT_DEFINE_double(aaEff, 0.7, "multiplicity model parameter AA Efficiency (taken from fit)");
SCT_DEFINE_double(aaMult, 540, "multiplicity model parameter AA Average 0-5% multiplicity (taken from fit)");

// function to check if a file should be used by doing some simple substring matching
bool AcceptFile(std::string file, sct::GlauberSpecies species, sct::CollisionEnergy energy) {
  // check to make sure it is a glauber root file
  if (file.find("glauber") == 0 && file.find(".root") == file.size() - 5) {
    
    // check if it has the proper deformation
    if ((FLAGS_deformation == true && file.find("deformed_true") == std::string::npos) ||
        (FLAGS_deformation == false && file.find("deformed_false") == std::string::npos))
      return false;
    
    // now check to make sure it is for the proper nucleus
    if (file.find(sct::speciesString[species]) != std::string::npos &&
        file.find(std::to_string(static_cast<int>(energy))) != std::string::npos) {
      return true;
    }
  }
  return false;
}

int main(int argc, char* argv[]) {
  
  // set help message
  std::string usage = "Runs systematic evaluation of centrality definition and ";
  usage += " refmult reweighting";
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
    case 2760 : energy = sct::CollisionEnergy::E2760; break;
    case 200 :  energy = sct::CollisionEnergy::E200; break;
    case 62 :   energy = sct::CollisionEnergy::E62; break;
    case 39 :   energy = sct::CollisionEnergy::E39; break;
    case 27 :   energy = sct::CollisionEnergy::E27; break;
    case 19 :   energy = sct::CollisionEnergy::E19; break;
    case 14 :   energy = sct::CollisionEnergy::E14; break;
    case 11 :   energy = sct::CollisionEnergy::E11; break;
    case 7 :    energy = sct::CollisionEnergy::E7; break;
    default :
      LOG(ERROR) << "requested unknown energy: " << FLAGS_energy << ": exiting";
      return 1;
  }
  
  // list of the glauber modifications to compare to default
  sct::GlauberMod default_mod = sct::GlauberMod::Nominal;
  std::vector<sct::GlauberMod> modifiers{sct::GlauberMod::Large, sct::GlauberMod::Small,
    sct::GlauberMod::LargeXSec, sct::GlauberMod::SmallXSec, sct::GlauberMod::Gauss};
  std::vector<std::string> mod_string{"large", "small", "largexsec", "smallxsec", "gauss"};
  
  // build output directory if it doesn't exist, using boost::filesystem
  boost::filesystem::path dir(FLAGS_outDir);
  boost::filesystem::create_directories(dir);
  
  // now find all the glauber output files in the specified directory
  boost::filesystem::path glauber_path(FLAGS_glauberDir);
  std::vector<std::string> glauber_input;
  if(boost::filesystem::is_directory(glauber_path)) {
    for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(glauber_path), {})) {
      std::string file = entry.path().filename().string();
      std::transform(file.begin(), file.end(), file.begin(), ::tolower);

      if (AcceptFile(file, species, energy))
        glauber_input.push_back(entry.path().string());
    }
  }
  
  // find the nominal input
  std::string glauber_input_default;
  for (auto& entry : glauber_input) {
    if (entry.find(sct::modToString[sct::GlauberMod::Nominal]) != std::string::npos) {
      glauber_input_default = entry;
      glauber_input.erase(std::find(glauber_input.begin(), glauber_input.end(), entry));
      break;
    }
  }
  
  if (glauber_input_default.empty()) {
    LOG(ERROR) << "could not find glauber file with default (nominal) settings, exiting";
    return 1;
  }
  
  LOG(INFO) << "Default (nominal) glauber results will be taken from: " << glauber_input_default;
  
  // check the size
  if (glauber_input.size() != modifiers.size()) {
    LOG(ERROR) << "Incorrect number of input glauber files found in input directiory: " << FLAGS_glauberDir;
    LOG(ERROR) << "Found(" << glauber_input.size() << ") modified glauber runs: " << glauber_input;
    LOG(ERROR) << "expected files for(" << modifiers.size() << ") modifiers: " << mod_string;
    LOG(ERROR) << "exiting";
    return 1;
  }

  LOG(INFO) << "Found (" << glauber_input.size() << ") systematic modifications: " << glauber_input;

  LOG(INFO) << "Calculating centrality & running systematic evaluation of glauber model";
  sct::GlauberSystematics systematics(glauber_input_default, glauber_input);
  
  gflags::ShutDownCommandLineFlags();
  return 0;
}
