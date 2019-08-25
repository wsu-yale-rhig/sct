#include "sct/systematics/glauber_systematics.h"
#include "sct/lib/enumerations.h"
#include "sct/lib/logging.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/histogram_info.h"

#include <algorithm>
#include <exception>

namespace sct {

GlauberSystematics::GlauberSystematics(const string& out_filename)
    : out_filename_(out_filename),
      do_npp_variations_(false),
      use_unit_weight_(false) {
  if (out_filename_.empty() || !openOutput()) {
    LOG(ERROR) << "could not open output file: " << out_filename_;
  }
}

GlauberSystematics::~GlauberSystematics() {}

bool GlauberSystematics::setOutputFile(const string& out_filename) {
  out_filename_ = out_filename;
  if (out_filename_.empty() || !openOutput()) {
    LOG(ERROR) << "could not open output file: " << out_filename_;
    return false;
  }
  return true;
}

bool GlauberSystematics::load(const string& glauber_filename) {
  GlauberMod tag;
  try {
    tag = findFileTag(glauber_filename);
  } catch (std::invalid_argument e) {
    LOG(ERROR) << "Could not identify glauber systematic variation from "
                  "filename: add tag manually with "
               << "Load(glauber_file, tag). Loading failed for "
               << glauber_filename;
    return false;
  }

  return loadFile(glauber_filename, tag);
}

bool GlauberSystematics::load(const std::set<string>& glauber_filenames) {
  bool status = true;
  GlauberMod tag;
  for (auto& filename : glauber_filenames) {
    try {
      tag = findFileTag(filename);
    } catch (std::invalid_argument e) {
      LOG(ERROR) << "Could not identify glauber systematic variation from "
                    "filename: add tag manually with "
                 << "Load(glauber_file, tag). Loading failed for " << filename;
      status = false;
      continue;
    }

    if (loadFile(filename, tag) == false) {
      status = false;
    }
  }
  return status;
}

bool GlauberSystematics::load(const string& glauber_filename, GlauberMod tag) {
  return loadFile(glauber_filename, tag);
}

void GlauberSystematics::setMultiplicityModel(double npp, double k, double x,
                                              double pp_eff, double aa_eff,
                                              double aa_cent, double trig_eff,
                                              bool const_eff) {
  mult_model_[GlauberMod::Nominal] = make_unique<MultiplicityModel>(
      npp, k, x, pp_eff, aa_eff, aa_cent, trig_eff, const_eff);
  mult_model_[GlauberMod::LargeNpp] = make_unique<MultiplicityModel>(
      npp * 1.05, k, x - 0.02, pp_eff, aa_eff, aa_cent, trig_eff, const_eff);
  mult_model_[GlauberMod::SmallNpp] = make_unique<MultiplicityModel>(
      npp * 0.95, k, x + 0.02, pp_eff, aa_eff, aa_cent, trig_eff, const_eff);
}

void GlauberSystematics::setCentralityBins(std::vector<double> cent_def) {
  if (cent_def.size() + 1 !=
      HistogramInfo::instance().bins(GlauberObservable::Centrality)) {
    LOG(ERROR) << "centrality bin does not match the one defined in "
                  "HistogramInfo: can not set";
    return;
  }
  centrality_lower_bounds_ = cent_def;
  std::sort(centrality_lower_bounds_.begin(), centrality_lower_bounds_.end(),
            std::greater<double>());
}

bool GlauberSystematics::run() {
  // initialize all histogram containers, output directories, etc
  initialize();

  // loop over all files/variations
  for (auto& var : variations_) {
    // load the tree
    GlauberMod key = var.first;
    auto& file_ptr = var.second;
    runFile(key, file_ptr);
  }
  return true;
}

bool GlauberSystematics::write() {
  // check that the file is open
  if (out_file_.get() == nullptr || !out_file_->IsOpen()) {
    LOG(ERROR) << "output file is not open - can not write to disk";
    return false;
  }

  impact_parameter_.write(out_file_.get());
  n_part_.write(out_file_.get());
  n_coll_.write(out_file_.get());
  multiplicity_.write(out_file_.get());
  rp_area_.write(out_file_.get());
  pp_area_.write(out_file_.get());
  rp_ecc_.write(out_file_.get());
  rp_ecc_mult_.write(out_file_.get());

  for (auto& hist : pp_ecc_) hist.write(out_file_.get());
  for (auto& hist : pp_ecc_mult_) hist.write(out_file_.get());

  return true;
}

void GlauberSystematics::initialize() {
  // clear all histograms before (re-)initializing
  clearHistograms();

  // set specific observables to calculate the 2nd, 4th and 6th order cumulant
  for (auto& container : pp_ecc_) container.calculateCumulants();
  for (auto& container : pp_ecc_mult_) container.calculateCumulants();

  // initialize the y axis for our observables
  impact_parameter_.init(GlauberObservable::B);
  n_part_.init(GlauberObservable::Npart);
  n_coll_.init(GlauberObservable::Ncoll);
  multiplicity_.init(GlauberObservable::Multiplicity);
  rp_area_.init(GlauberObservable::RPArea);
  pp_area_.init(GlauberObservable::PPArea);
  rp_ecc_.init(GlauberObservable::RP2Ecc);
  rp_ecc_mult_.init(GlauberObservable::RP2Ecc);
  pp_ecc_[0].init(GlauberObservable::PP2Ecc);
  pp_ecc_[1].init(GlauberObservable::PP3Ecc);
  pp_ecc_[2].init(GlauberObservable::PP4Ecc);
  pp_ecc_mult_[0].init(GlauberObservable::PP2Ecc);
  pp_ecc_mult_[1].init(GlauberObservable::PP3Ecc);
  pp_ecc_mult_[2].init(GlauberObservable::PP4Ecc);

  for (auto& var : variations_) {
    auto& key = var.first;
    auto* file = var.second.get();

    // add the variation to the histogram containers
    impact_parameter_.add(key);
    n_part_.add(key);
    n_coll_.add(key);
    multiplicity_.add(key);
    rp_area_.add(key);
    pp_area_.add(key);
    rp_ecc_.add(key);
    rp_ecc_mult_.add(key);
    for (auto& entry : pp_ecc_) entry.add(key);
    for (auto& entry : pp_ecc_mult_) entry.add(key);
  }
}

void GlauberSystematics::clearHistograms() {
  impact_parameter_.clear();
  n_part_.clear();
  n_coll_.clear();
  multiplicity_.clear();
  rp_area_.clear();
  pp_area_.clear();
  rp_ecc_.clear();
  rp_ecc_mult_.clear();
  for (auto& container : pp_ecc_) container.clear();
  for (auto& container : pp_ecc_mult_) container.clear();
}

bool GlauberSystematics::loadFile(const string& filename, GlauberMod tag) {
  shared_ptr<TFile> file = make_shared<TFile>(filename.c_str(), "READ");
  if (file.get() == nullptr || !file->IsOpen() ||
      file->Get("event") == nullptr || file->Get("header") == nullptr) {
    LOG(ERROR) << "Could not open glauber file: " << filename
               << ": file loading failed";
    return false;
  }

  if (variations_.find(tag) != variations_.end()) {
    LOG(WARNING) << "variation " << glauberModToString[tag]
                 << " already exists as key to " << variations_[tag]->GetFile()
                 << ", overwriting this with filename";
  }

  variations_.insert({tag, std::move(file)});

  return true;
}

bool findCaseInsensitive(std::string data, std::string toSearch,
                         size_t pos = 0) {
  // Convert complete given String to lower case
  std::transform(data.begin(), data.end(), data.begin(), ::tolower);
  // Convert complete given Sub String to lower case
  std::transform(toSearch.begin(), toSearch.end(), toSearch.begin(), ::tolower);
  // Find sub string in given string
  return data.find(toSearch, pos) != string::npos;
}

GlauberMod GlauberSystematics::findFileTag(const string& filename) {
  string tag;
  for (auto& mod : glauberModString) {
    if (findCaseInsensitive(filename, mod)) {
      if (!tag.empty()) {
        LOG(WARNING) << "warning - filename: " << filename
                     << " is ambiguous. String matching returns "
                     << " multiple possible variations. Overwriting " << tag
                     << " with " << mod;
      }
      tag = mod;
    }
  }

  if (tag.empty())
    throw std::invalid_argument(MakeString(
        "Could not identify a Glauber Modification from filename ", filename));
  return stringToGlauberMod[tag];
}

int GlauberSystematics::getCentrality(double multiplicity) {
  for (int i = 0; i < centrality_lower_bounds_.size(); ++i) {
    if (multiplicity > centrality_lower_bounds_[i]) return i;
  }
  return centrality_lower_bounds_.size();
}

bool GlauberSystematics::openOutput() {
  out_file_ = make_unique<TFile>(out_filename_.c_str(), "READ");
  if (out_file_.get() == nullptr || !out_file_->IsOpen()) return false;
  return true;
}

void GlauberSystematics::runFile(GlauberMod key, shared_ptr<TFile> file) {
  GlauberTree reader(GlauberTree::TreeMode::Read);
  reader.read(file);

  // do event loop
  for (int i = 0; i < reader.getEntries(); ++i) {
    reader.getEntry(i);

    // create an event dictionary
    EventDict<double> event_dict;
    event_dict[GlauberObservable::B] = reader.B();
    event_dict[GlauberObservable::Ncoll] = reader.nColl();
    event_dict[GlauberObservable::Npart] = reader.nPart();
    event_dict[GlauberObservable::RP2Ecc] = reader.RP2Ecc(GlauberWeight::NPart);
    event_dict[GlauberObservable::PP2Ecc] = reader.PP2Ecc(GlauberWeight::NPart);
    event_dict[GlauberObservable::PP3Ecc] = reader.PP3Ecc(GlauberWeight::NPart);
    event_dict[GlauberObservable::PP4Ecc] = reader.PP4Ecc(GlauberWeight::NPart);
    event_dict[GlauberObservable::RPArea] = reader.RPArea(GlauberWeight::NPart);
    event_dict[GlauberObservable::PPArea] = reader.PPArea(GlauberWeight::NPart);

    // if reader multiplicity is set, we can use that - unless we are varying
    // npp
    double mult;
    double cent;
    if (key != GlauberMod::LargeNpp && key != GlauberMod::SmallNpp &&
        reader.multiplicity() != 0) {
      mult = reader.multiplicity();
      cent = getCentrality(mult);
    } else if (mult_model_[GlauberMod::Nominal].get() != nullptr &&
               (key == GlauberMod::LargeNpp || key == GlauberMod::SmallNpp)) {
      mult = mult_model_[key]->multiplicity(reader.nColl(), reader.nPart());
      cent = getCentrality(mult);
    } else {
      mult = -1;
      cent = -1;
    }

    event_dict[GlauberObservable::Multiplicity] = mult;
    event_dict[GlauberObservable::Centrality] = cent;

    // weights
    double unit_weight = 1.0;
    double weight = use_unit_weight_ ? 1.0 : reader.multiplicity();

    // fill in histograms with unit weight
    impact_parameter_.fillEvent(key, event_dict, unit_weight);
    n_part_.fillEvent(key, event_dict, unit_weight);
    n_coll_.fillEvent(key, event_dict, unit_weight);
    multiplicity_.fillEvent(key, event_dict, unit_weight);

    // fill histograms with multiplicity weight
    // if use_unit_weight_ is false
    rp_area_.fillEvent(key, event_dict, weight);
    pp_area_.fillEvent(key, event_dict, weight);
    rp_ecc_.fillEvent(key, event_dict, weight);
    for (auto& hist : pp_ecc_) hist.fillEvent(key, event_dict, weight);

    // fill event dict with multiplicity weighted values
    event_dict[GlauberObservable::RP2Ecc] =
        reader.RP2Ecc(GlauberWeight::Multiplicity);
    event_dict[GlauberObservable::PP2Ecc] =
        reader.PP2Ecc(GlauberWeight::Multiplicity);
    event_dict[GlauberObservable::PP3Ecc] =
        reader.PP3Ecc(GlauberWeight::Multiplicity);
    event_dict[GlauberObservable::PP4Ecc] =
        reader.PP4Ecc(GlauberWeight::Multiplicity);

    // fill histograms for multiplicity weighted eccentricities
    rp_ecc_mult_.fillEvent(key, event_dict, weight);
    for (auto& hist : pp_ecc_mult_) hist.fillEvent(key, event_dict, weight);
  }
}

}  // namespace sct
