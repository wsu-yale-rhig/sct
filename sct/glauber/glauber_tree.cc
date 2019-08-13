#include "sct/glauber/glauber_tree.h"

#include "sct/lib/logging.h"
#include "sct/lib/math.h"

namespace sct {

GlauberTree::GlauberTree(TreeMode mode, const string &filename)
    : file_(), header_(nullptr), event_(nullptr), nameA_(new string()),
      nameB_(new string()) {
  switch (mode) {
  case TreeMode::Write:
    createBranches();
    return;
  case TreeMode::Read:
    if (!filename.empty()) {
      VLOG(1) << "Tree set to read mode - opening " << filename;
      open(filename);
    }
    return;
  default:
    LOG(ERROR) << "library error: TreeMode not recognized in GlauberTree";
  }
}

GlauberTree::~GlauberTree() {
  if (file_.get() == nullptr) {
    delete header_;
    delete event_;
  }
}

void GlauberTree::clearEvent() {
  b_ = 0;
  nPart_ = 0;
  nColl_ = 0;

  for (unsigned i = 0; i < n_nuclei; ++i) {
    theta_[i] = 0;
    phi_[i] = 0;
  }

  for (unsigned i = 0; i < nGlauberWeights; ++i) {
    sumX_[i] = 0;
    sumY_[i] = 0;
    sumX2_[i] = 0;
    sumY2_[i] = 0;
    sumXY_[i] = 0;
    rp2ecc_[i] = 0;
    pp2ecc_[i] = 0;
    pp3ecc_[i] = 0;
    pp4ecc_[i] = 0;
    pp2_[i] = 0;
    pp3_[i] = 0;
    pp4_[i] = 0;
  }
}

bool GlauberTree::open(const string &filename) {
  VLOG(1) << "sct::GlauberTree opening file to read";
  file_ = make_shared<TFile>(filename.c_str(), "READ");
  if (!file_ || !file_->IsOpen()) {
    LOG(ERROR) << "failed to open root file at " << filename
               << " does path exist?";
    return false;
  }
  event_ = (TTree *)file_->Get("event");
  header_ = (TTree *)file_->Get("header");
  if (event_ == nullptr) {
    LOG(ERROR) << "event tree not found in root file: was the"
               << " file written by a sct GlauberTree?";
    return false;
  }
  if (header_ == nullptr) {
    LOG(ERROR) << "header tree not found in root file: was the"
               << " file written by a sct GlauberTree?";
    return false;
  }
  loadBranches();
  return true;
}

bool GlauberTree::read(shared_ptr<TFile> file) {
  // sanity checks on input
  if (file.get() == nullptr) {
    LOG(ERROR) << "can not read from uninitialized file";
    return false;
  }
  VLOG(1) << "sct::GlauberTree reading from file: " << file->GetName();
  if (!file->IsOpen()) {
    LOG(ERROR) << "can not read from unopened file";
    return 0;
  }

  // load trees
  file_ = file;
  event_ = (TTree *)file_->Get("event");
  header_ = (TTree *)file_->Get("header");
  if (event_ == nullptr) {
    LOG(ERROR) << "event tree not found in root file: was the"
               << " file written by a sct GlauberTree?";
    return false;
  }
  if (header_ == nullptr) {
    LOG(ERROR) << "header tree not found in root file: was the"
               << " file written by a sct GlauberTree?";
    return false;
  }
  loadBranches();
  return true;
}

void GlauberTree::write() {
  if (header_ != nullptr)
    header_->Write();
  if (event_ != nullptr)
    event_->Write();
}

double GlauberTree::RPArea(GlauberWeight id) {
  // sigma x^2 = <x^2> - <x>^2
  int index = static_cast<int>(id);
  double sigmax2 = sumX2_[index] - pow(sumX_[index], 2.0);
  double sigmay2 = sumY2_[index] - pow(sumY_[index], 2.0);
  return pi * sqrt(sigmax2 * sigmay2);
}

double GlauberTree::PPArea(GlauberWeight id) {
  // sigma x^2 = <x^2> - <x>^2
  int index = static_cast<int>(id);
  double sigmax2 = sumX2_[index] - pow(sumX_[index], 2.0);
  double sigmay2 = sumY2_[index] - pow(sumY_[index], 2.0);
  double sigmaxy = sumXY_[index] - sumX_[index] * sumY_[index];
  return pi * sqrt(sigmax2 * sigmay2 - pow(sigmaxy, 2.0));
}

void GlauberTree::loadBranches() {
  // load in the event TTree first
  VLOG(1) << "Loading event tree from file";
  event_->SetBranchAddress("b", &b_);
  event_->SetBranchAddress("npart", &nPart_);
  event_->SetBranchAddress("ncoll", &nColl_);
  event_->SetBranchAddress("nspec", &nSpec_);
  event_->SetBranchAddress("multiplicity", &multiplicity_);
  event_->SetBranchAddress("theta", theta_);
  event_->SetBranchAddress("phi", phi_);
  event_->SetBranchAddress("sumx", sumX_);
  event_->SetBranchAddress("sumy", sumY_);
  event_->SetBranchAddress("sumx2", sumX2_);
  event_->SetBranchAddress("sumy2", sumY2_);
  event_->SetBranchAddress("sumxy", sumXY_);
  event_->SetBranchAddress("rp2ecc", rp2ecc_);
  event_->SetBranchAddress("pp2ecc", pp2ecc_);
  event_->SetBranchAddress("pp3ecc", pp3ecc_);
  event_->SetBranchAddress("pp4ecc", pp4ecc_);
  event_->SetBranchAddress("pp2", pp2_);
  event_->SetBranchAddress("pp3", pp3_);
  event_->SetBranchAddress("pp4", pp4_);

  VLOG(1) << "Loading header tree from file";
  header_->SetBranchAddress("nameA", &nameA_);
  header_->SetBranchAddress("nameB", &nameB_);
  header_->SetBranchAddress("massnumberA", &massNumberA_);
  header_->SetBranchAddress("massnumberB", &massNumberB_);
  header_->SetBranchAddress("radiusA", &radiusA_);
  header_->SetBranchAddress("radiusB", &radiusB_);
  header_->SetBranchAddress("skindepthA", &skinDepthA_);
  header_->SetBranchAddress("skindepthB", &skinDepthB_);
  header_->SetBranchAddress("beta2A", &beta2A_);
  header_->SetBranchAddress("beta2B", &beta2B_);
  header_->SetBranchAddress("beta4A", &beta4A_);
  header_->SetBranchAddress("beta4B", &beta4B_);
  header_->SetBranchAddress("sigmann", &sigmaNN_);
  header_->SetBranchAddress("sqrtsnn", &sqrtSNN_);
  header_->SetBranchAddress("repulsiond", &repulsionD_);
  header_->SetBranchAddress("totalxsec", &totalXsec_);
  header_->SetBranchAddress("totalxsecerror", &totalXsecError_);
  header_->SetBranchAddress("smearhardcore", &smearHardCore_);
  header_->SetBranchAddress("smeargaussian", &smearGaussian_);
  header_->SetBranchAddress("collisionhardcore", &collisionHardCore_);
  header_->SetBranchAddress("collisiongaussian", &collisionGaussian_);
  header_->SetBranchAddress("bmax", &bMax_);
  header_->SetBranchAddress("bmin", &bMin_);
  header_->SetBranchAddress("eventsaccepted", &eventsAccepted_);
  header_->SetBranchAddress("eventsthrown", &eventsThrown_);
}

void GlauberTree::createBranches() {
  // clear any old trees
  delete event_;
  delete header_;

  // create the event tree
  event_ = new TTree("event", "event-wise sct record");
  event_->Branch("b", &b_);
  event_->Branch("npart", &nPart_);
  event_->Branch("ncoll", &nColl_);
  event_->Branch("nspec", &nSpec_);
  event_->Branch("multiplicity", &multiplicity_);
  event_->Branch("theta", theta_, "theta[2]/D");
  event_->Branch("phi", phi_, "phi[2]/D");
  event_->Branch("sumx", sumX_, "sumx[3]/D");
  event_->Branch("sumy", sumY_, "sumy[3]/D");
  event_->Branch("sumx2", sumX2_, "sumx2[3]/D");
  event_->Branch("sumy2", sumY2_, "sumy2[3]/D");
  event_->Branch("sumxy", sumXY_, "sumxy[3]/D");
  event_->Branch("rp2ecc", rp2ecc_, "rp2ecc[3]/D");
  event_->Branch("pp2ecc", pp2ecc_, "pp2ecc[3]/D");
  event_->Branch("pp3ecc", pp3ecc_, "pp3ecc[3]/D");
  event_->Branch("pp4ecc", pp4ecc_, "pp4ecc[3]/D");
  event_->Branch("pp2", pp2_, "pp2[3]/D");
  event_->Branch("pp3", pp3_, "pp3[3]/D");
  event_->Branch("pp4", pp4_, "pp4[3]/D");

  // create new header tree
  header_ = new TTree("header", "generator settings");
  header_->Branch("nameA", &nameA_);
  header_->Branch("nameB", &nameB_);
  header_->Branch("massnumberA", &massNumberA_);
  header_->Branch("massnumberB", &massNumberB_);
  header_->Branch("radiusA", &radiusA_);
  header_->Branch("radiusB", &radiusB_);
  header_->Branch("skindepthA", &skinDepthA_);
  header_->Branch("skindepthB", &skinDepthB_);
  header_->Branch("beta2A", &beta2A_);
  header_->Branch("beta2B", &beta2B_);
  header_->Branch("beta4A", &beta4A_);
  header_->Branch("beta4B", &beta4B_);
  header_->Branch("sigmann", &sigmaNN_);
  header_->Branch("sqrtsnn", &sqrtSNN_);
  header_->Branch("repulsiond", &repulsionD_);
  header_->Branch("totalxsec", &totalXsec_);
  header_->Branch("totalxsecerror", &totalXsecError_);
  header_->Branch("smearhardcore", &smearHardCore_);
  header_->Branch("smeargaussian", &smearGaussian_);
  header_->Branch("collisionhardcore", &collisionHardCore_);
  header_->Branch("collisiongaussian", &collisionGaussian_);
  header_->Branch("bmax", &bMax_);
  header_->Branch("bmin", &bMin_);
  header_->Branch("eventsaccepted", &eventsAccepted_);
  header_->Branch("eventsthrown", &eventsThrown_);
}

unsigned GlauberTree::getEntries() const {
  if (event_ != nullptr)
    return event_->GetEntries();
  return 0;
}

bool GlauberTree::getEntry(unsigned int idx) {
  if (event_ != nullptr)
    return event_->GetEntry(idx);
  return false;
}

bool GlauberTree::getHeaderEntry(unsigned int idx) {
  if (header_ != nullptr)
    return header_->GetEntry(idx);
  return false;
}

void GlauberTree::fill() {
  if (event_ != nullptr) {
    event_->Fill();
  } else {
    LOG(ERROR) << "attempting to fill uninitialized event tree: "
               << "have you called Open(), and did it succeed?";
  }
}

void GlauberTree::fillHeader() {
  if (header_ != nullptr) {
    header_->Fill();
  } else {
    LOG(ERROR) << "attempting to fill uninitialized header tree: "
               << "have you called Open(), and did it succeed?";
  }
}

} // namespace sct
