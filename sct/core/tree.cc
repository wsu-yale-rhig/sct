// sct/core/tree.cc

#include "sct/core/tree.hh"

#include "sct/core/logging.hh"

namespace sct {
  
  Tree::Tree(const string& filename)
  : file_(), header_(nullptr), event_(nullptr) {
    if (filename == "") {
      createBranches();
    }
    else {
      open(filename);
    }
  }
  
  Tree::~Tree() {
    if (file_.get() == nullptr) {
      delete header_;
      delete event_;
    }
  }
  
  void Tree::clearEvent() {
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
      eccRP2_[i] = 0;
      eccPP2_[i] = 0;
      eccPP3_[i] = 0;
      eccPP4_[i] = 0;
      pp2_[i] = 0;
      pp3_[i] = 0;
      pp4_[i] = 0;
    }
  }
  
  bool Tree::open(const string& filename) {
    
    VLOG(1) << "sct::tree opening file to read";
    file_ = make_unique<TFile>(filename.c_str(), "READ");
    if (!file_ || !file_->IsOpen()) {
      LOG(ERROR) << "failed to open root file at " << filename
                   << " does path exist?";
      return false;
    }
    event_ = (TTree*) file_->Get("event");
    header_ = (TTree*) file_->Get("header");
    if (event_ == nullptr) {
      LOG(ERROR) << "event tree not found in root file: was the"
                 << " file written by a sct tree?";
      return false;
    }
    if (header_ == nullptr) {
      LOG(ERROR) << "header tree not found in root file: was the"
                 << " file written by a sct tree?";
      return false;
    }
    loadBranches();
    return true;
  }
  
  void Tree::write() {
    if (header_ != nullptr)
      header_->Write();
    if (event_ != nullptr)
      event_->Write();
  }
  
  void Tree::loadBranches() {
    // load in the event TTree first
    VLOG(1) << "Loading event tree from file";
    event_->SetBranchAddress("b", &b_);
    event_->SetBranchAddress("npart", &nPart_);
    event_->SetBranchAddress("ncoll", &nColl_);
    event_->SetBranchAddress("nspec", &nSpec_);
    event_->SetBranchAddress("theta", theta_);
    event_->SetBranchAddress("phi", phi_);
    event_->SetBranchAddress("sumx", sumX_);
    event_->SetBranchAddress("sumy", sumY_);
    event_->SetBranchAddress("sumx2", sumX2_);
    event_->SetBranchAddress("sumy2", sumY2_);
    event_->SetBranchAddress("sumxy", sumXY_);
    event_->SetBranchAddress("eccrp2", eccRP2_);
    event_->SetBranchAddress("eccpp2", eccPP2_);
    event_->SetBranchAddress("eccpp3", eccPP3_);
    event_->SetBranchAddress("eccpp4", eccPP4_);
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
    header_->SetBranchAddress("skinDepthB", &skinDepthB_);
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
  
  void Tree::createBranches() {
    // clear any old trees
    delete event_;
    delete header_;
    
    // create the event tree
    event_ = new TTree("event", "event-wise sct record");
    event_->Branch("b", &b_);
    event_->Branch("npart", &nPart_);
    event_->Branch("ncoll", &nColl_);
    event_->Branch("nspec", &nSpec_);
    event_->Branch("theta", theta_, "theta[2]/D");
    event_->Branch("phi", phi_, "phi[2]/D");
    event_->Branch("sumx", sumX_, "sumx[3]/D");
    event_->Branch("sumy", sumY_, "sumy[3]/D");
    event_->Branch("sumx2", sumX2_, "sumx2[3]/D");
    event_->Branch("sumy2", sumY2_, "sumy2[3]/D");
    event_->Branch("sumxy", sumXY_, "sumxy[3]/D");
    event_->Branch("eccrp2", eccRP2_, "eccrp2[3]/D");
    event_->Branch("eccpp2", eccPP2_, "eccpp2[3]/D");
    event_->Branch("eccpp3", eccPP3_, "eccpp3[3]/D");
    event_->Branch("eccpp4", eccPP4_, "eccpp4[3]/D");
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
  
  unsigned Tree::getEntries() const {
    return event_->GetEntries();
  }
  
  bool Tree::getEntry(unsigned int idx) {
    return event_->GetEntry(idx);
  }
  
  bool Tree::getHeaderEntry(unsigned int idx) {
    return header_->GetEntry(idx);
  }
  
  void Tree::fill() {
    if (event_ != nullptr) {
      event_->Fill();
    }
    else {
      LOG(ERROR) << "attempting to fill uninitialized event tree: "
                 << "have you called Open(), and did it succeed?";
    }
  }
  
  void Tree::fillHeader() {
    if (header_ != nullptr) {
      header_->Fill();
    }
    else {
      LOG(ERROR) << "attempting to fill uninitialized header tree: "
                 << "have you called Open(), and did it succeed?";
    }
  }
  
} // namespace sct
