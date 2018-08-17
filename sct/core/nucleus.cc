// sct/core/nucleus.cc

#include "sct/core/nucleus.hh"

#include "sct/core/logging.hh"
#include "sct/utils/functions.hh"
#include "sct/utils/random.hh"

#include "TF2.h"

namespace sct {
  
  Nucleus::Nucleus()
  : name_(""), mass_number_(0), radius_(0), skin_depth_(0),
  beta2_(0), beta4_(0), smear_(NucleonSmearing::None),
  sigmaNN_(0.0), repulsion_distance_(0.0), random_orientation_(true),
  nucleus_theta_(0.0), nucleus_phi_(0.0), b_(0.0),
  generated_rcostheta_(nullptr), generated_position_(nullptr),
  generated_smear_(nullptr), smeared_position_(nullptr) { }
  
  Nucleus::Nucleus(unsigned massNumber, double radius, double skinDepth,
                   double beta2, double beta4)
  : name_(""), mass_number_(massNumber), radius_(radius), skin_depth_(skinDepth),
  beta2_(beta2), beta4_(beta4), smear_(NucleonSmearing::None),
  sigmaNN_(0.0), repulsion_distance_(0.0), random_orientation_(true),
  nucleus_theta_(0.0), nucleus_phi_(0.0), b_(0.0) {
    
    // reserve space for the correct number of nucleons
    nucleons_.reserve(massNumber);
    
    initWoodsSaxon();
    initHistograms();
  }
  
  Nucleus::Nucleus(const Nucleus& rhs)
  : name_(rhs.name()), mass_number_(rhs.massNumber()), radius_(rhs.radius()),
  beta2_(rhs.beta2()), beta4_(rhs.beta4()), smear_(rhs.nucleonSmearing()),
  sigmaNN_(rhs.nnCrossSection()), repulsion_distance_(rhs.repulsionDistance()),
  random_orientation_(rhs.randomOrientation()), nucleus_theta_(rhs.nucleusTheta()),
  nucleus_phi_(rhs.nucleusPhi()), b_(rhs.impactParameter()) {
    
    // reserve space for the correct number of nucleons
    nucleons_.reserve(mass_number_);
    for (int i = 0; i < rhs.size(); ++i)
      nucleons_.push_back(rhs[i]);
    
    initWoodsSaxon();
    initHistograms();
  }
  
  Nucleus::~Nucleus() {
    
  }
  
  void Nucleus::clear() {
    nucleons_.clear();
    nucleus_theta_ = 0.0;
    nucleus_phi_ = 0.0;
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
    
    mass_number_ = massNumber;
    radius_ = radius;
    skin_depth_ = skinDepth;
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
    if ((beta2_ != 0.0 || beta4_ != 0.0) && random_orientation_) {
      nucleus_theta_ = acos(Random::instance().centeredUniform());
      nucleus_phi_ = Random::instance().zeroToPi() * 2.0 - pi;
    }
    else {
      nucleus_theta_ = 0.0;
      nucleus_phi_ = 0.0;
    }
    
    // since repulsion can lead to effectively an infinite loop
    // (if 4/3 * pi * radius^3 <<< 4/3 * pi * repulsionDistance^3 * mass_number_)
    // we will fail out after 10 * mass_number_ of entries
    int tries = 0;
    while (nucleons_.size() < mass_number_) {
      if ( tries > mass_number_ * 10) {
        LOG(ERROR) << "Repeated failure to create nucleus of radius: " << radius_;
        LOG(ERROR) << "with " << mass_number_ << " nucleons, with a nucleon repulsion "
                   << "distance of " << repulsion_distance_;
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
      smearing_profile_.release();
      smear_ = NucleonSmearing::None;
      sigmaNN_ = 0;
      return;
    }
    
    smear_ = smear;
    sigmaNN_ = sigmaNN;
    smearing_profile_.release();
    
    switch (smear) {
      case NucleonSmearing::HardCore : {
        double rMax = sqrt(sigmaNN_ / pi);
        smearing_profile_ = make_unique<TF3>(MakeString("hardCoreSmear_", Random::instance().counter()).c_str(),
                                            StepFunction, -rMax, rMax, -rMax, rMax, -rMax, rMax, 1);
        smearing_profile_->SetParameter(0, sigmaNN_);
        break;
      }
      case NucleonSmearing::Gaussian : {
        double rMax = sigmaNN_ * 5.0;
        double sigma = 0.79 / sqrt(3.0);
        smearing_profile_ = make_unique<TF3>(MakeString("gaussianSmear_", Random::instance().counter()).c_str(),
                                            Gaussian, -rMax, rMax, -rMax, rMax, -rMax, rMax, 1);
        smearing_profile_->SetParameter(0, sigma);
        break;
      }
      default:
        smearing_profile_ = nullptr;
        break;
    }
  }
  
  void Nucleus::setRepulsionDistance(double fm) {
    if (fm >= 0)
      repulsion_distance_ = fm;
    else {
      LOG(ERROR) << "Negative nucleon repulsion distance does nothing, setting to zero";
      repulsion_distance_ = 0.0;
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
      nucleon.set(smearedPosition, b, nucleus_theta_, nucleus_phi_, true);
      
      // if no repulsion, we don't need to check against other nucleon
      // positions, just add and break out
      if (repulsion_distance_ <= 0.0 || nucleons_.size() == 0) {
        nucleons_.push_back(nucleon);
        
        // record QA data
        generated_position_->Fill(position.X(), position.Y(), position.Z());
        generated_smear_->Fill(smearing.X(), smearing.Y(), smearing.Z());
        smeared_position_->Fill(smearedPosition.X(), smearedPosition.Y(), smearedPosition.Z());
        break;
      }
      else {
        // otherwise, we have to check every generated nucleon to see
        // if it overlaps
        bool collision = false;
        for (auto nucleon2 : nucleons_) {
          if (nucleon.deltaR(nucleon2) <= repulsion_distance_)
            collision = true;
        }
        if (collision == false) {
          // if there was no overlap, add the nucleon
          nucleons_.push_back(nucleon);
          
          // record QA data
          generated_position_->Fill(position.X(), position.Y(), position.Z());
          generated_smear_->Fill(smearing.X(), smearing.Y(), smearing.Z());
          smeared_position_->Fill(smearedPosition.X(), smearedPosition.Y(), smearedPosition.Z());
          break;
        }
      }
    }
  }
  
  void Nucleus::initWoodsSaxon() {
    woods_saxon_.release();
    
    // either 1D or 2D depending on the deformation parameters
    if (beta2_ == 0.0 && beta4_ == 0.0) {
      woods_saxon_ = make_unique<TF1>(MakeString("woods_saxon_", Random::instance().counter()).c_str(),
                                     WoodsSaxonSpherical, 0, 20, 2);
      woods_saxon_->SetParName(0, "Radius");
      woods_saxon_->SetParameter(0, radius_);
      woods_saxon_->SetParName(1, "Skin depth");
      woods_saxon_->SetParameter(1, skin_depth_);
      woods_saxon_->SetNpx(400);
    }
    else {
      woods_saxon_ = make_unique<TF2>(MakeString("woods_saxon_", Random::instance().counter()).c_str(),
                                     WoodsSaxonDeformed, 0, 20, -1.0, 1.0, 4);
      woods_saxon_->SetParName(0, "Radius");
      woods_saxon_->SetParameter(0, radius_);
      woods_saxon_->SetParName(1, "Skin depth");
      woods_saxon_->SetParameter(1, skin_depth_);
      woods_saxon_->SetParName(2, "#beta_{2}");
      woods_saxon_->SetParameter(2, beta2_);
      woods_saxon_->SetParName(3, "#beta_{4}");
      woods_saxon_->SetParameter(3, beta4_);
      woods_saxon_->SetNpx(400);
      ((TF2*)woods_saxon_.get())->SetNpy(400);
    }
  }
  
  void Nucleus::initHistograms() {
    generated_rcostheta_.release();
    generated_position_.release();
    generated_smear_.release();
    smeared_position_.release();
    
    // initialize QA histograms
    generated_rcostheta_   = make_unique<TH2D>(MakeString("r_costheta_",
                                              Random::instance().counter()).c_str(),
                                              ";R;cos(#theta)",
                                              400, 0, 3 * radius_,
                                              100, -1.0, 1.0);
    generated_rcostheta_->SetDirectory(0);
    generated_position_ = make_unique<TH3D>(MakeString("nucleonpos_",
                                           Random::instance().counter()).c_str(),
                                           ";dx;dy;dz",
                                           100, -3.0 * radius_, 3.0 * radius_,
                                           100, -3.0 * radius_, 3.0 * radius_,
                                           100, -3.0 * radius_, 3.0 * radius_);
    generated_position_->SetDirectory(0);
    generated_smear_    = make_unique<TH3D>(MakeString("nucleonsmear_",
                                           Random::instance().counter()).c_str(),
                                           ";dx;dy;dz",
                                           100, -3.0 * radius_, 3.0 * radius_,
                                           100, -3.0 * radius_, 3.0 * radius_,
                                           100, -3.0 * radius_, 3.0 * radius_);
    generated_smear_->SetDirectory(0);
    smeared_position_   = make_unique<TH3D>(MakeString("nucleonsmearedpos_",
                                           Random::instance().counter()).c_str(),
                                           ";dx;dy;dz",
                                           100, -3.0 * radius_, 3.0 * radius_,
                                           100, -3.0 * radius_, 3.0 * radius_,
                                           100, -3.0 * radius_, 3.0 * radius_);
    smeared_position_->SetDirectory(0);
  }
  
  TVector3 Nucleus::generateNucleonPosition() {
    double r, theta, phi;
    // if the nucleus is not deformed, we are using a 1D woods-saxon
    if (beta2_ == 0.0 && beta4_ == 0.0) {
      r = woods_saxon_->GetRandom();
      theta = acos(Random::instance().centeredUniform());
    }
    // otherwise we use a 2D woods-saxon
    else {
      double cosTheta;
      ((TF2*) woods_saxon_.get())->GetRandom2(r, cosTheta);
      theta = acos(cosTheta);
    }
    
    // record the generated R & cos theta for QA
    generated_rcostheta_->Fill(r, cos(theta), 1.0 / pow(r, 2.0));
    
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
    if (smearing_profile_ != nullptr) {
      smearing_profile_->GetRandom3(dx, dy, dz);
    }
    return TVector3(dx, dy, dz);
  }
  
  const Nucleon& Nucleus::operator[](unsigned idx) const {
    if (idx > nucleons_.size())
      throw "Out of bounds access";
    return nucleons_[idx];
  }
  
  Nucleon& Nucleus::operator[](unsigned idx) {
    if (idx > nucleons_.size())
      throw "Out of bounds access";
    return nucleons_[idx];
  }
} // namespace sct
