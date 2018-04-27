// sct/core/nucleus.cc

#include "sct/core/nucleus.hh"

#include "sct/core/logging.hh"
#include "sct/utils/functions.hh"
#include "sct/utils/random.hh"

#include "TF2.h"

namespace sct {
  
  Nucleus::Nucleus()
  : name_(""), massNumber_(0), radius_(0), skinDepth_(0),
  beta2_(0), beta4_(0), smear_(NucleonSmearing::None),
  sigmaNN_(0.0), repulsionDistance_(0.0), randomOrientation_(true),
  nucleusTheta_(0.0), nucleusPhi_(0.0), b_(0.0),
  generatedRCosTheta_(nullptr), generatedPosition_(nullptr),
  generatedSmear_(nullptr), smearedPosition_(nullptr) { }
  
  Nucleus::Nucleus(unsigned massNumber, double radius, double skinDepth,
                   double beta2, double beta4)
  : name_(""), massNumber_(massNumber), radius_(radius), skinDepth_(skinDepth),
  beta2_(beta2), beta4_(beta4), smear_(NucleonSmearing::None),
  sigmaNN_(0.0), repulsionDistance_(0.0), randomOrientation_(true),
  nucleusTheta_(0.0), nucleusPhi_(0.0), b_(0.0) {
    
    // reserve space for the correct number of nucleons
    nucleons_.reserve(massNumber);
    
    initWoodsSaxon();
    initHistograms();
  }
  
  Nucleus::Nucleus(const Nucleus& rhs)
  : name_(rhs.name()), massNumber_(rhs.massNumber()), radius_(rhs.radius()),
  beta2_(rhs.beta2()), beta4_(rhs.beta4()), smear_(rhs.nucleonSmearing()),
  sigmaNN_(rhs.nnCrossSection()), repulsionDistance_(rhs.repulsionDistance()),
  randomOrientation_(rhs.randomOrientation()), nucleusTheta_(rhs.nucleusTheta()),
  nucleusPhi_(rhs.nucleusPhi()), b_(rhs.impactParameter()) {
    
    // reserve space for the correct number of nucleons
    nucleons_.reserve(massNumber_);
    for (int i = 0; i < rhs.size(); ++i)
      nucleons_.push_back(rhs[i]);
    
    initWoodsSaxon();
    initHistograms();
  }
  
  Nucleus::~Nucleus() {
    
  }
  
  void Nucleus::clear() {
    nucleons_.clear();
    nucleusTheta_ = 0.0;
    nucleusPhi_ = 0.0;
    b_ = 0.0;
  }
  
  void Nucleus::setParameters(unsigned massNumber, double radius,
                              double skinDepth, double beta2, double beta4) {
    if (massNumber == 0) {
      LOG(ERROR) << "mass number can not be zero, nucleus in invalid state";
      return;
    }
    if (radius == 0) {
      LOG(ERROR) << "mass number can not be zero, nucleus in invalid state";
      return;
    }
    
    clear();
    nucleons_.reserve(massNumber);
    
    massNumber_ = massNumber;
    radius_ = radius;
    skinDepth_ = skinDepth;
    beta2_ = beta2;
    beta4_ = beta4;
    
    initWoodsSaxon();
    initHistograms();
  }
  
  bool Nucleus::generate(double b) {
    // clear the last results
    clear();
    
    // set the impact parameter
    b_ = b;
    
    // get an orientation for the nucleus (theta & phi), if it is deformed,
    // else, keep at zero
    if ((beta2_ != 0.0 || beta4_ != 0.0) && randomOrientation_) {
      nucleusTheta_ = acos(Random::instance().centeredUniform());
      nucleusPhi_ = Random::instance().zeroToPi() * 2.0 - pi;
    }
    else {
      nucleusTheta_ = 0.0;
      nucleusPhi_ = 0.0;
    }
    
    // since repulsion can lead to effectively an infinite loop
    // (if 4/3 * pi * radius^3 <<< 4/3 * pi * repulsionDistance^3 * massNumber_)
    // we will fail out after 10 * massNumber_ of entries
    int tries = 0;
    while (nucleons_.size() < massNumber_) {
      if ( tries > massNumber_ * 10) {
        LOG(ERROR) << "Repeated failure to create nucleus of radius: " << radius_;
        LOG(ERROR) << "with " << massNumber_ << " nucleons, with a nucleon repulsion "
                   << "distance of " << repulsionDistance_;
        return false;
      }
      
      addNucleon(b);
      ++tries;
    }
    return true;
  }
  
  void Nucleus::setNucleonSmearing(NucleonSmearing smear, double sigmaNN) {
    if (sigmaNN < 0) {
      LOG(ERROR) << "Nucleon smearing requested with negative inelastic"
                 << " cross section, setting to zero";
      smearingProfile_.release();
      smear_ = NucleonSmearing::None;
      sigmaNN_ = 0;
      return;
    }
    
    smear_ = smear;
    sigmaNN_ = sigmaNN;
    smearingProfile_.release();
    
    switch (smear) {
      case NucleonSmearing::HardCore : {
        double rMax = sqrt(sigmaNN_ / pi);
        smearingProfile_ = make_unique<TF3>(MakeString("hardCoreSmear_", Random::instance().counter()).c_str(),
                                            StepFunction, -rMax, rMax, -rMax, rMax, -rMax, rMax, 1);
        smearingProfile_->SetParameter(0, sigmaNN_);
        break;
      }
      case NucleonSmearing::Gaussian : {
        double rMax = sigmaNN_ * 5.0;
        double sigma = 0.79 / sqrt(3.0);
        smearingProfile_ = make_unique<TF3>(MakeString("gaussianSmear_", Random::instance().counter()).c_str(),
                                            Gaussian, -rMax, rMax, -rMax, rMax, -rMax, rMax, 1);
        smearingProfile_->SetParameter(0, sigma);
        break;
      }
      default:
        smearingProfile_ = nullptr;
        break;
    }
  }
  
  void Nucleus::setRepulsionDistance(double fm) {
    if (fm >= 0)
      repulsionDistance_ = fm;
    else {
      LOG(ERROR) << "Negative nucleon repulsion distance does nothing, setting to zero";
      repulsionDistance_ = 0.0;
    }
  }
  
  void Nucleus::addNucleon(double b) {
    
    // attempt to add one nucleon to the nucleus
    Nucleon nucleon;
    
    // now try to add the nucleon to the nucleus - "try" because if there
    // is a repulsion it can fail, so try 5 times, and fail out if it doesn't fit
    int tries = 0;
    
    while (tries < 5) {
      tries++;

      // generate a random position
      TVector3 position = generateNucleonPosition();
      // smear, if nucleon position smearing is turned on
      TVector3 smearing = smear();
    
      // get x, y, z
      TVector3 smearedPosition = position + smearing;
    
      // create the nucleon
      nucleon.set(smearedPosition, b, nucleusTheta_, nucleusPhi_, true);
      
      // if no repulsion, we don't need to check against other nucleon
      // positions, just add and break out
      if (repulsionDistance_ <= 0.0 || nucleons_.size() == 0) {
        nucleons_.push_back(nucleon);
        
        // record QA data
        generatedPosition_->Fill(position.X(), position.Y(), position.Z());
        generatedSmear_->Fill(smearing.X(), smearing.Y(), smearing.Z());
        smearedPosition_->Fill(smearedPosition.X(), smearedPosition.Y(), smearedPosition.Z());
        break;
      }
      else {
        // otherwise, we have to check every generated nucleon to see
        // if it overlaps
        bool collision = false;
        for (auto nucleon2 : nucleons_) {
          if (nucleon.deltaR(nucleon2) <= repulsionDistance_)
            collision = true;
        }
        if (collision == false) {
          // if there was no overlap, add the nucleon
          nucleons_.push_back(nucleon);
          
          // record QA data
          generatedPosition_->Fill(position.X(), position.Y(), position.Z());
          generatedSmear_->Fill(smearing.X(), smearing.Y(), smearing.Z());
          smearedPosition_->Fill(smearedPosition.X(), smearedPosition.Y(), smearedPosition.Z());
          break;
        }
      }
    }
  }
  
  void Nucleus::initWoodsSaxon() {
    woodsSaxon_.release();
    
    // either 1D or 2D depending on the deformation parameters
    if (beta2_ == 0.0 && beta4_ == 0.0) {
      woodsSaxon_ = make_unique<TF1>(MakeString("woodsSaxon_", Random::instance().counter()).c_str(),
                                     WoodsSaxonSpherical, 0, 20, 2);
      woodsSaxon_->SetParName(0, "Radius");
      woodsSaxon_->SetParameter(0, radius_);
      woodsSaxon_->SetParName(1, "Skin depth");
      woodsSaxon_->SetParameter(1, skinDepth_);
      woodsSaxon_->SetNpx(400);
    }
    else {
      woodsSaxon_ = make_unique<TF2>(MakeString("woodsSaxon_", Random::instance().counter()).c_str(),
                                     WoodsSaxonDeformed, 0, 20, -1.0, 1.0, 4);
      woodsSaxon_->SetParName(0, "Radius");
      woodsSaxon_->SetParameter(0, radius_);
      woodsSaxon_->SetParName(1, "Skin depth");
      woodsSaxon_->SetParameter(1, skinDepth_);
      woodsSaxon_->SetParName(2, "#beta_{2}");
      woodsSaxon_->SetParameter(2, beta2_);
      woodsSaxon_->SetParName(3, "#beta_{4}");
      woodsSaxon_->SetParameter(3, beta4_);
      woodsSaxon_->SetNpx(400);
      ((TF2*)woodsSaxon_.get())->SetNpy(400);
    }
  }
  
  void Nucleus::initHistograms() {
    generatedRCosTheta_.release();
    generatedPosition_.release();
    generatedSmear_.release();
    smearedPosition_.release();
    
    // initialize QA histograms
    generatedRCosTheta_   = make_unique<TH2D>(MakeString("r_costheta_",
                                              Random::instance().counter()).c_str(),
                                              ";R;cos(#theta)",
                                              400, 0, 3 * radius_,
                                              100, -1.0, 1.0);
    generatedRCosTheta_->SetDirectory(0);
    generatedPosition_ = make_unique<TH3D>(MakeString("nucleonpos_",
                                           Random::instance().counter()).c_str(),
                                           ";dx;dy;dz",
                                           100, -3.0 * radius_, 3.0 * radius_,
                                           100, -3.0 * radius_, 3.0 * radius_,
                                           100, -3.0 * radius_, 3.0 * radius_);
    generatedPosition_->SetDirectory(0);
    generatedSmear_    = make_unique<TH3D>(MakeString("nucleonsmear_",
                                           Random::instance().counter()).c_str(),
                                           ";dx;dy;dz",
                                           100, -3.0 * radius_, 3.0 * radius_,
                                           100, -3.0 * radius_, 3.0 * radius_,
                                           100, -3.0 * radius_, 3.0 * radius_);
    generatedSmear_->SetDirectory(0);
    smearedPosition_   = make_unique<TH3D>(MakeString("nucleonsmearedpos_",
                                           Random::instance().counter()).c_str(),
                                           ";dx;dy;dz",
                                           100, -3.0 * radius_, 3.0 * radius_,
                                           100, -3.0 * radius_, 3.0 * radius_,
                                           100, -3.0 * radius_, 3.0 * radius_);
    smearedPosition_->SetDirectory(0);
  }
  
  TVector3 Nucleus::generateNucleonPosition() {
    double r, theta, phi;
    // if the nucleus is not deformed, we are using a 1D woods-saxon
    if (beta2_ == 0.0 && beta4_ == 0.0) {
      r = woodsSaxon_->GetRandom();
      theta = acos(Random::instance().centeredUniform());
    }
    // otherwise we use a 2D woods-saxon
    else {
      double cosTheta;
      ((TF2*) woodsSaxon_.get())->GetRandom2(r, cosTheta);
      theta = acos(cosTheta);
    }
    
    // record the generated R & cos theta for QA
    generatedRCosTheta_->Fill(r, cos(theta), 1.0 / pow(r, 2.0));
    
    // finally, generate phi flat in -pi < phi < pi
    phi = Random::instance().zeroToPi() * 2.0 - pi;
    
    TVector3 ret;
    ret.SetMagThetaPhi(r, theta, phi);
    
    return ret;
  }
  
  TVector3 Nucleus::smear() {
    // if smearing is turned on, then we will smear the nucleon position
    // in r, theta, phi
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    if (smearingProfile_ != nullptr) {
      smearingProfile_->GetRandom3(dx, dy, dz);
    }
    return TVector3(dx, dy, dz);
  }
  
  const Nucleon& Nucleus::operator[](unsigned idx) const throw (const char *) {
    if (idx > nucleons_.size())
      throw "Out of bounds access";
    return nucleons_[idx];
  }
  
  Nucleon& Nucleus::operator[](unsigned idx) throw (const char *) {
    if (idx > nucleons_.size())
      throw "Out of bounds access";
    return nucleons_[idx];
  }
} // namespace sct
