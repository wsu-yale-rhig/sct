// sct/core/mc_glauber.cc

#include "sct/core/mc_glauber.hh"

#include "sct/core/logging.hh"
#include "sct/utils/nucleus_info.hh"
#include "sct/utils/random.hh"

namespace sct {
  
  MCGlauber::MCGlauber(GlauberSpecies species1, GlauberSpecies species2,
                       CollisionEnergy energy, GlauberMod mod,
                       bool deformationA, bool deformationB)
  : eventsGenerated_(0), eventsAccepted_(0), bMin_(0), bMax_(20),
  energy_(static_cast<double>(energy)), modification_(mod) {
    
    // get nominal settings for nuclear parameters
    unsigned massNumberA = NucleusInfo::instance().massNumber(species1);
    unsigned massNumberB = NucleusInfo::instance().massNumber(species2);
    double radiusA = NucleusInfo::instance().radius(species1);
    double radiusB = NucleusInfo::instance().radius(species2);
    double skinDepthA = NucleusInfo::instance().skinDepth(species1);
    double skinDepthB = NucleusInfo::instance().skinDepth(species2);
    double beta2A = 0.0;
    double beta2B = 0.0;
    double beta4A = 0.0;
    double beta4B = 0.0;
    if (deformationA) {
      beta2A = NucleusInfo::instance().beta2(species1);
      beta4A = NucleusInfo::instance().beta4(species1);
    }
    if (deformationB) {
      beta2B = NucleusInfo::instance().beta2(species2);
      beta4B = NucleusInfo::instance().beta4(species2);
    }
    double xsec = lookupXSec(energy);
    
    // apply any requested systematic variation of the parameters
    switch (mod) {
      case GlauberMod::Large :
        radiusA += NucleusInfo::instance().radiusError(species1);
        radiusB += NucleusInfo::instance().radiusError(species2);
        skinDepthA -= NucleusInfo::instance().skinDepthError(species1);
        skinDepthB -= NucleusInfo::instance().skinDepthError(species2);
        break;
      case GlauberMod::Small :
        radiusA -= NucleusInfo::instance().radiusError(species1);
        radiusB -= NucleusInfo::instance().radiusError(species2);
        skinDepthA += NucleusInfo::instance().skinDepthError(species1);
        skinDepthB += NucleusInfo::instance().skinDepthError(species2);
        break;
      case GlauberMod::LargeXSec :
        xsec += 0.1;
        break;
      case GlauberMod::SmallXSec :
        xsec -= 0.1;
        break;
      default :
        break;
    }
    
    // initialize nuclei and output tree
    init(massNumberA, radiusA, skinDepthA, beta2A, beta4A, massNumberB,
         radiusB, skinDepthB, beta2B, beta4B);
    
    initOutput();
    initQA();
    
    // initialize collision
    collision_.setNNCrossSection(xsec);
    nucleusA_->setName(NucleusInfo::instance().name(species1));
    nucleusB_->setName(NucleusInfo::instance().name(species2));
  }
  
  
  MCGlauber::MCGlauber(unsigned massNumber, double radius, double skinDepth,
                       double beta2, double beta4, double inelasticXSec, double energy)
  : eventsGenerated_(0), eventsAccepted_(0), bMin_(0), bMax_(20), energy_(energy),
  modification_(GlauberMod::Nominal) {
    
    init(massNumber, radius, skinDepth, beta2, beta4, massNumber, radius, skinDepth,
         beta2, beta4);
    initOutput();
    initQA();
    collision_.setNNCrossSection(inelasticXSec);
    nucleusA_->setName(MakeString(massNumber));
    nucleusB_->setName(MakeString(massNumber));
  }
  
  MCGlauber::MCGlauber(unsigned massNumberA, double radiusA, double skinDepthA,
                       double beta2A, double beta4A, unsigned massNumberB,
                       double radiusB, double skinDepthB, double beta2B, double beta4B,
                       double inelasticXSec, double energy)
  : eventsGenerated_(0), eventsAccepted_(0), bMin_(0), bMax_(20), energy_(energy),
  modification_(GlauberMod::Nominal) {
    
    init(massNumberA, radiusA, skinDepthA, beta2A, beta4A, massNumberB,
         radiusB, skinDepthB, beta2B, beta4B);
    initOutput();
    initQA();
    collision_.setNNCrossSection(inelasticXSec);
    nucleusA_->setName(MakeString(massNumberA));
    nucleusB_->setName(MakeString(massNumberB));
  }
  
  MCGlauber::~MCGlauber() {
    
  }
  
  void MCGlauber::setImpactParameterRange(double bMin, double bMax) {
    if (bMin < 0) {
      bMin = 0.0;
    }
    if (bMax < bMin) {
      LOG(ERROR) << "Requested impact parameter range doesn't make sense: "
                 <<" setting default range [0, 20] fm.";
      bMin = 0.0;
      bMax = 20.0;
    }
    bMin_ = bMin;
    bMax_ = bMax;
    initQA();
  }
  
  void MCGlauber::setRepulsionDistance(double fm) {
    if (fm < 0) {
      LOG(ERROR) << "Nucleon repulsion distance set negative: no effect,"
                 << "setting to zero";
      nucleusA_->setRepulsionDistance(0);
      nucleusB_->setRepulsionDistance(0);
      return;
    }
    nucleusA_->setRepulsionDistance(fm);
    nucleusB_->setRepulsionDistance(fm);
  }
  
  void MCGlauber::setSmearing(NucleonSmearing smear) {
    nucleusA_->setNucleonSmearing(smear, collision_.NNCrossSection());
    nucleusB_->setNucleonSmearing(smear, collision_.NNCrossSection());
  }
  
  void MCGlauber::setCollisionProfile(CollisionProfile profile) {
    collision_.setCollisionProfile(profile);
  }
  
  void MCGlauber::initOutput() {
    clear();
    tree_ = make_unique<Tree>();
  }
  
  void MCGlauber::initQA() {
    generated_ip_.reset();
    accepted_ip_.reset();
    generated_ip_ = make_unique<TH1D>(MakeString("generated_ip_",
                                                 Random::instance().counter()).c_str(),
                                      ";impact parameter[fm]",
                                      100, bMin_, bMax_);
    generated_ip_->SetDirectory(0);
    accepted_ip_  = make_unique<TH1D>(MakeString("accepted_ip_",
                                                 Random::instance().counter()).c_str(),
                                      ";impact parameter[fm]",
                                      100, bMin_, bMax_);
    accepted_ip_->SetDirectory(0);
  }
  
  void MCGlauber::clear() {
    if (nucleusA_.get() != nullptr)
      nucleusA_->clear();
    if (nucleusB_.get() != nullptr)
      nucleusB_->clear();
    collision_.clear();
  }
  
  // run the glauber MC for N events
  void MCGlauber::run(unsigned N) {
    initOutput();
    
    while (eventsAccepted_ < N) {
      eventsGenerated_++;
      if (generate()) {
        eventsAccepted_++;
        accepted_ip_->Fill(tree_->B());
        if (eventsAccepted_%100 == 0) {
          VLOG(1) << "MCGlauber::run::event = " << eventsAccepted_;
          VLOG(1) << "(accepted/generated, Npart, Ncoll, b) = "
                  << MakeString("(", (double)eventsAccepted_/eventsGenerated_, ", ",
                                tree_->nPart(), ", ", tree_->nColl(),
                                ", ", tree_->B(), ")");
        }
      }
    }
    writeHeader();
  }
  
  // write the tree to file
  void MCGlauber::writeHeader() {
    
    // fill header
    tree_->setNameNucleusA(nucleusA_->name());
    tree_->setNameNucleusB(nucleusB_->name());
    tree_->setMassNumberA(nucleusA_->massNumber());
    tree_->setMassNumberB(nucleusB_->massNumber());
    tree_->setRadiusA(nucleusA_->radius());
    tree_->setRadiusB(nucleusB_->radius());
    tree_->setSkinDepthA(nucleusA_->skinDepth());
    tree_->setSkinDepthB(nucleusB_->skinDepth());
    tree_->setBeta2A(nucleusA_->beta2());
    tree_->setBeta2B(nucleusB_->beta2());
    tree_->setBeta4A(nucleusA_->beta4());
    tree_->setBeta4B(nucleusB_->beta4());
    tree_->setSigmaNN(collision_.NNCrossSection());
    tree_->setSqrtSNN(energy_);
    tree_->setRepulsionD(nucleusA_->repulsionDistance());
    tree_->setSmearHardCore(nucleusA_->nucleonSmearing() == NucleonSmearing::HardCore);
    tree_->setSmearGaussian(nucleusA_->nucleonSmearing() == NucleonSmearing::Gaussian);
    tree_->setCollisionHardCore(collision_.collisionProfile() == CollisionProfile::HardCore);
    tree_->setCollisionGaussian(collision_.collisionProfile() == CollisionProfile::Gaussian);
    tree_->setBMax(bMax_);
    tree_->setBMin(bMin_);
    tree_->setNEventsAccepted(eventsAccepted_);
    tree_->setNEventsThrown(eventsGenerated_);
    
    // calculate total cross section
    double xMin = eventsAccepted_ * pow(bMax_, 2) * pi * 10 / eventsGenerated_;
    double xMax = eventsAccepted_ * pow(bMin_, 2) * pi * 10 / eventsGenerated_;
    double xTotal = xMax - xMin;
    double xErr = (xMax - xMin) * sqrt(1.0 / eventsAccepted_ + 1.0 / eventsGenerated_);
    tree_->setTotalXsec(xTotal);
    tree_->setTotalXsecError(xErr);
    
    // write header, and write tree to file
    tree_->fillHeader();
  }
  
  bool MCGlauber::generate() {
    tree_->clearEvent();
    clear();
    // generate impact parameter
    double b = bMin_ + Random::instance().linear() * (bMax_ - bMin_);
    generated_ip_->Fill(b);
    
    // generate two new nuclei
    nucleusA_->generate(b/2.0);
    nucleusB_->generate(-b/2.0);
    
    // collide the nuclei: if there are no collisions, we are done
    if (!collision_.collide(*nucleusA_.get(), *nucleusB_.get()))
      return false;
  
    // fill the event record in the tree
    tree_->setB(b);
    tree_->setNpart(collision_.nPart());
    tree_->setNcoll(collision_.nColl());
    tree_->setNspectators(collision_.spectators());
    tree_->setTheta(0, nucleusA_->nucleusTheta());
    tree_->setTheta(0, nucleusA_->nucleusPhi());
    tree_->setTheta(1, nucleusB_->nucleusTheta());
    tree_->setTheta(1, nucleusB_->nucleusPhi());
    for (int i = 0; i < nGlauberWeights; ++i) {
      tree_->setSumX(i, collision_.averageX()[i]);
      tree_->setSumY(i, collision_.averageY()[i]);
      tree_->setSumX2(i, collision_.averageX2()[i]);
      tree_->setSumY2(i, collision_.averageY2()[i]);
      tree_->setSumXY(i, collision_.averageXY()[i]);
      tree_->setEccRP2(i, collision_.reactionPlane2Ecc()[i]);
      tree_->setEccPP2(i, collision_.partPlane2Ecc()[i]);
      tree_->setEccPP3(i, collision_.partPlane3Ecc()[i]);
      tree_->setEccPP4(i, collision_.partPlane4Ecc()[i]);
      tree_->setPP2(i, collision_.partPlane2()[i]);
      tree_->setPP3(i, collision_.partPlane3()[i]);
      tree_->setPP4(i, collision_.partPlane4()[i]);
    }
    
    tree_->fill();
    return true;
  }

  double MCGlauber::lookupXSec(CollisionEnergy energy) {
    switch (energy) {
      case CollisionEnergy::E2760 :
        return 6.4;
        break;
      case CollisionEnergy::E200 :
        return 4.2;
        break;
      case CollisionEnergy::E62 :
        return 3.6;
        break;
      case CollisionEnergy::E39 :
        return 3.4;
        break;
      case CollisionEnergy::E27 :
        return 3.3;
        break;
      case CollisionEnergy::E19 :
        return 3.2;
        break;
      case CollisionEnergy::E14 :
        return 3.15;
        break;
      case CollisionEnergy::E11 :
        return 3.12;
        break;
      case CollisionEnergy::E7  :
        return 3.08;
        break;
      default :
        LOG(ERROR) << "No inelastic nucleon-nucleon cross section for collision energy";
        LOG(ERROR) << "Either implement in lookupXSec(CollisionEnergy) or set xsec by hand";
        return 0.0;
    }
  }
  
  void MCGlauber::init(unsigned massNumberA, double radiusA, double skinDepthA,
                       double beta2A, double beta4A, unsigned massNumberB,
                       double radiusB, double skinDepthB, double beta2B,
                       double beta4B) {
    nucleusA_.release();
    nucleusB_.release();
    nucleusA_ = make_unique<Nucleus>(massNumberA, radiusA, skinDepthA, beta2A, beta4A);
    nucleusB_ = make_unique<Nucleus>(massNumberB, radiusB, skinDepthB, beta2B, beta4B);
  }
  
} // namespace sct
