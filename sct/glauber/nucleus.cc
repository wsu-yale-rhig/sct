#include "sct/glauber/nucleus.h"
#include "sct/lib/assert.h"
#include "sct/lib/logging.h"
#include "sct/lib/math.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/functions.h"
#include "sct/utils/nucleus_info.h"
#include "sct/utils/random.h"

#include "TF2.h"

namespace sct {

Nucleus::Nucleus()
    : name_(""), mass_number_(0), smear_(NucleonSmearing::None),
      repulsion_distance_(0.0), random_orientation_(true), nucleus_theta_(0.0),
      nucleus_phi_(0.0), b_(0.0), generated_rcos_theta_(nullptr),
      generated_position_(nullptr), generated_smear_(nullptr),
      smeared_position_(nullptr) {}

Nucleus::Nucleus(GlauberSpecies species, GlauberMod mod, bool deformed)
    : name_(""), mass_number_(0), smear_(NucleonSmearing::None),
      repulsion_distance_(0.0), random_orientation_(true), nucleus_theta_(0.0),
      nucleus_phi_(0.0), b_(0.0), generated_rcos_theta_(nullptr),
      generated_position_(nullptr), generated_smear_(nullptr),
      smeared_position_(nullptr) {
  name_ = NucleusInfo::instance().name(species);
  setParameters(species, mod, deformed);
}

Nucleus::Nucleus(unsigned mass_number, parameter_list params,
                 NucleonPDF::PDF pdf)
    : name_(""), mass_number_(mass_number), smear_(NucleonSmearing::None),
      repulsion_distance_(0.0), random_orientation_(true), nucleus_theta_(0.0),
      nucleus_phi_(0.0), b_(0.0) {
  setParameters(mass_number, params, pdf);
}

Nucleus::Nucleus(const Nucleus &rhs)
    : name_(rhs.name()), mass_number_(rhs.massNumber()),
      smear_(rhs.nucleonSmearing()),
      repulsion_distance_(rhs.repulsionDistance()),
      random_orientation_(rhs.randomOrientation()),
      nucleus_theta_(rhs.nucleusTheta()), nucleus_phi_(rhs.nucleusPhi()),
      b_(rhs.impactParameter()) {
  // reserve space for the correct number of nucleons
  nucleons_.reserve(mass_number_);
  for (int i = 0; i < rhs.size(); ++i)
    nucleons_.push_back(rhs[i]);

  nucleon_pdf_.init(rhs.nucleon_pdf_.PDFForm(), rhs.nucleon_pdf_.parameters());
  initHistograms();
}

Nucleus::~Nucleus() {}

void Nucleus::clear() {
  nucleons_.clear();
  nucleus_theta_ = 0.0;
  nucleus_phi_ = 0.0;
  b_ = 0.0;
}

bool Nucleus::setParameters(GlauberSpecies species, GlauberMod mod,
                            bool deformed) {
  // reserve space for correct number of nucleons
  nucleons_.reserve(NucleusInfo::instance().massNumber(species));
  mass_number_ = NucleusInfo::instance().massNumber(species);

  initHistograms();

  // initialize the pdf
  return nucleon_pdf_.init(species, mod, deformed);
}

bool Nucleus::setParameters(unsigned mass_number, parameter_list params,
                            NucleonPDF::PDF pdf) {
  if (mass_number == 0) {
    LOG(ERROR) << "mass number can not be zero, nucleus in invalid state";
    return false;
  }

  clear();
  nucleons_.reserve(mass_number);
  mass_number_ = mass_number;

  initHistograms();

  // reset the nuclear PDF
  return nucleon_pdf_.init(pdf, params);
}

bool Nucleus::generate(double b) {
  // clear the last results
  clear();

  // set the impact parameter
  b_ = b;

  // get an orientation for the nucleus (theta & phi), if it is deformed,
  // else, keep at zero
  if (nucleon_pdf_.deformed() && random_orientation_) {
    nucleus_theta_ = acos(Random::instance().centeredUniform());
    nucleus_phi_ = Random::instance().zeroToPi() * 2.0 - pi;
  } else {
    nucleus_theta_ = 0.0;
    nucleus_phi_ = 0.0;
  }

  // special case for deuteron - we need to place the nuclei opposite to each
  // other, so we do not randomly sample each nucleon. Instead, we sample a
  // single nucleon and place the other opposite
  if (mass_number_ == 0 && nucleon_pdf_.PDFForm() == NucleonPDF::PDF::Hulthen)
    return generateDeuteron();

  // since repulsion can lead to effectively an infinite loop
  // (if 4/3 * pi * radius^3 ~ 4/3 * pi * repulsionDistance^3 * mass_number_)
  // we will fail out after 10 * mass_number_ of entries

  return generateNucleus();
}

void Nucleus::setNucleonSmearing(NucleonSmearing smear, double smear_area) {
  if (smear_area < 0) {
    LOG(ERROR) << "Nucleon smearing requested with negative inelastic"
               << " cross section, setting to zero";
    smearing_profile_.reset(nullptr);
    smear_ = NucleonSmearing::None;
    return;
  }

  smear_ = smear;

  switch (smear) {
  case NucleonSmearing::HardCore: {
    double rMax = sqrt(smear_area / pi);
    smearing_profile_ = make_unique<TF3>(
        MakeString("hardCoreSmear_", Random::instance().counter()).c_str(),
        StepFunction, -rMax, rMax, -rMax, rMax, -rMax, rMax, 1);
    smearing_profile_->SetParameter(0, smear_area);
    break;
  }
  case NucleonSmearing::Gaussian: {
    double sigma = 0.79 / sqrt(3.0);
    double rMax = sigma * 5.0;
    smearing_profile_ = make_unique<TF3>(
        MakeString("gaussianSmear_", Random::instance().counter()).c_str(),
        Gaussian, -rMax, rMax, -rMax, rMax, -rMax, rMax, 1);
    smearing_profile_->SetParameter(0, sigma);
    break;
  }
  default:
    smearing_profile_.reset(nullptr);
    break;
  }
}

void Nucleus::setRepulsionDistance(double fm) {
  if (fm >= 0)
    repulsion_distance_ = fm;
  else {
    LOG(ERROR)
        << "Negative nucleon repulsion distance does nothing, setting to zero";
    repulsion_distance_ = 0.0;
  }
}

void Nucleus::rotateAndOffset(Nucleon &n) {
  TVector3 pos = n.position();

  // rotate in the frame of the nucleus
  pos.RotateY(nucleusTheta());
  pos.RotateZ(nucleusPhi());

  // translate by impact parameter along the x axis - moves the nucleon into the
  // collision CM frame
  pos.SetX(pos.X() + b_);

  n.set(pos);
}

void Nucleus::addNucleon(double b) {
  // attempt to add one nucleon to the nucleus
  Nucleon nucleon;

  // now try to add the nucleon to the nucleus - "try" because if there
  // is a repulsion it can fail, so try 5 times, and fail out if it doesn't
  // "fit"
  int tries = 0;

  while (tries < 5) {
    tries++;

    // generate a random position
    TVector3 position = generateNucleonPosition();

    // smear, if nucleon position smearing is turned on
    TVector3 smearing = smear();

    // the nominal position is the addition of raw position + smearing
    TVector3 smeared_position = position + smearing;

    // create the nucleon
    nucleon.set(smeared_position);

    // rotate by the nucleus orientation and then offset by the impact parameter
    rotateAndOffset(nucleon);

    // if a minimum repulsion distance is set, we have to check each existing
    // nucleon for overlap - else, we can skip that relatively inefficient step
    bool collision = false;
    if (repulsion_distance_ > 0.0) {
      for (auto nucleon2 : nucleons_) {
        if (nucleon.deltaR(nucleon2) <= repulsion_distance_)
          collision = true;
      }
    }

    if (collision == false) {
      // if there was no overlap, add the nucleon
      nucleons_.push_back(nucleon);

      // record QA data
      generated_position_->Fill(position.X(), position.Y(), position.Z());
      generated_smear_->Fill(smearing.X(), smearing.Y(), smearing.Z());
      smeared_position_->Fill(smeared_position.X(), smeared_position.Y(),
                              smeared_position.Z());
      break;
    }
  }
}

void Nucleus::initHistograms() {

  // initialize QA histograms
  generated_rcos_theta_ = make_unique<TH2D>(
      MakeString("r_cos_theta_", Random::instance().counter()).c_str(),
      ";R;cos(#theta)", 400, 0, 20, 100, -1.0, 1.0);
  generated_rcos_theta_->SetDirectory(0);
  generated_position_ = make_unique<TH3D>(
      MakeString("nucleonpos_", Random::instance().counter()).c_str(),
      ";dx;dy;dz", 100, -20.0, 20.0, 100, -20.0, 20.0, 100, -20.0, 20.0);
  generated_position_->SetDirectory(0);
  generated_smear_ = make_unique<TH3D>(
      MakeString("nucleonsmear_", Random::instance().counter()).c_str(),
      ";dx;dy;dz", 100, -20.0, 20.0, 100, -20.0, 20.0, 100, -20.0, 20.0);
  generated_smear_->SetDirectory(0);
  smeared_position_ = make_unique<TH3D>(
      MakeString("nucleonsmearedpos_", Random::instance().counter()).c_str(),
      ";dx;dy;dz", 100, -20.0, 20.0, 100, -20.0, 20.0, 100, -20.0, 20.0);
  smeared_position_->SetDirectory(0);
}

TVector3 Nucleus::generateNucleonPosition() {
  double r, theta, phi;

  nucleon_pdf_.sample(r, theta, phi);

  // record the generated R & cos theta for QA
  generated_rcos_theta_->Fill(r, cos(theta), 1.0 / pow(r, 2.0));

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

const Nucleon &Nucleus::operator[](unsigned idx) const {
  SCT_ASSERT(idx < nucleons_.size(), "Out of bounds access");
  return nucleons_[idx];
}

Nucleon &Nucleus::operator[](unsigned idx) {
  SCT_ASSERT(idx < nucleons_.size(), "Out of bounds access");
  return nucleons_[idx];
}

bool Nucleus::generateDeuteron() {
  int tries = 0;
  while (tries < 5) {
    Nucleon nucleon_a;
    Nucleon nucleon_b;
    // generate a random position
    TVector3 position_a = generateNucleonPosition();
    // smear, if nucleon position smearing is turned on
    TVector3 smearing_a = smear();

    // the nominal position is the addition of raw position + smearing
    TVector3 smeared_position_a = position_a + smearing_a;
    TVector3 smeared_position_b = -smeared_position_a;

    // create the nucleon
    nucleon_a.set(smeared_position_a);
    nucleon_b.set(smeared_position_b);

    if (repulsion_distance_ > 0.0 &&
        nucleon_a.deltaR(nucleon_b) >= repulsion_distance_) {
      // rotate by the nucleus orientation and then offset by the impact
      // parameter
      rotateAndOffset(nucleon_a);
      rotateAndOffset(nucleon_b);
      nucleons_.push_back(nucleon_a);
      nucleons_.push_back(nucleon_b);
      return true;
    }
  }
  LOG(ERROR) << "Repeated failure to generate a deuteron nucleus with a "
                "Hulthen PDF: nucleon-nucleon repulsion distance ("
             << repulsion_distance_
             << ") is too large for the specified PDF parameters";
  return false;
}

bool Nucleus::generateNucleus() {
  int tries = 0;

  while (nucleons_.size() < mass_number_) {
    if (tries > mass_number_ * 10) {
      LOG(ERROR) << "Repeated failure to create nucleus "
                 << " with " << mass_number_
                 << " nucleons, with a nucleon repulsion "
                 << "distance of " << repulsion_distance_;
      clear();
      return false;
    }

    addNucleon(b_);
    ++tries;
  }
  return true;
}

} // namespace sct
