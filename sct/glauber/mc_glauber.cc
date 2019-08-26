#include "sct/glauber/mc_glauber.h"
#include "sct/lib/assert.h"
#include "sct/lib/logging.h"
#include "sct/lib/math.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/nucleus_info.h"
#include "sct/utils/random.h"

namespace sct {

MCGlauber::MCGlauber(GlauberSpecies species_A, GlauberSpecies species_B,
                     CollisionEnergy energy, GlauberMod mod, bool deformation_A,
                     bool deformation_B)
    : events_generated_(0), events_accepted_(0), b_min_(0), b_max_(20),
      energy_(static_cast<double>(energy)), modification_(mod) {

  double xsec = lookupXSec(energy);

  // apply any requested systematic variation of the parameters
  switch (mod) {
  case GlauberMod::LargeXSec:
    xsec += 0.1;
    break;
  case GlauberMod::SmallXSec:
    xsec -= 0.1;
    break;
  default:
    break;
  }

  init(species_A, species_B, mod, deformation_A, deformation_B);

  initOutput();
  initQA();

  collision_.setNNCrossSection(xsec);
  nucleusA_->setName(NucleusInfo::instance().name(species_A));
  nucleusB_->setName(NucleusInfo::instance().name(species_B));
}

MCGlauber::MCGlauber(unsigned mass_number, NucleonPDF::PDF pdf,
                     parameter_list params, double inelastic_xsec,
                     double energy)
    : events_generated_(0), events_accepted_(0), b_min_(0), b_max_(20),
      energy_(energy), modification_(GlauberMod::Nominal) {

  init(mass_number, pdf, params, mass_number, pdf, params);

  initOutput();
  initQA();

  collision_.setNNCrossSection(inelastic_xsec);
  nucleusA_->setName(MakeString(mass_number));
  nucleusB_->setName(MakeString(mass_number));
}

MCGlauber::MCGlauber(unsigned mass_number_A, NucleonPDF::PDF pdf_A,
                     parameter_list params_A, unsigned mass_number_B,
                     NucleonPDF::PDF pdf_B, parameter_list params_B,
                     double inelastic_xsec, double energy)
    : events_generated_(0), events_accepted_(0), b_min_(0), b_max_(20),
      energy_(energy), modification_(GlauberMod::Nominal) {

  initOutput();
  initQA();

  init(mass_number_A, pdf_A, params_A, mass_number_B, pdf_B, params_B);

  collision_.setNNCrossSection(inelastic_xsec);
  nucleusA_->setName(MakeString(mass_number_A));
  nucleusB_->setName(MakeString(mass_number_B));
}

MCGlauber::~MCGlauber() {}

void MCGlauber::setImpactParameterRange(double b_min, double b_max) {
  if (b_min < 0) {
    b_min = 0.0;
  }
  if (b_max < b_min) {
    LOG(ERROR) << "Requested impact parameter range doesn't make sense: "
               << " setting default range [0, 20] fm.";
    b_min = 0.0;
    b_max = 20.0;
  }
  b_min_ = b_min;
  b_max_ = b_max;
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

void MCGlauber::setMultiplicityModel(double npp, double k, double x,
                                     double pp_eff, double aa_eff,
                                     double aa_cent, double trig_eff,
                                     bool const_eff) {
  collision_.setMultiplicityModel(npp, k, x, pp_eff, aa_eff, aa_cent, trig_eff,
                                  const_eff);
}

void MCGlauber::initOutput() {
  clear();
  tree_ = make_unique<GlauberTree>(GlauberTree::TreeMode::Write);
}

void MCGlauber::initQA() {
  generated_ip_.reset();
  accepted_ip_.reset();
  generated_ip_ = make_unique<TH1D>(
      MakeString("generated_ip_", Counter::instance().counter()).c_str(),
      ";impact parameter[fm]", 100, b_min_, b_max_);
  generated_ip_->SetDirectory(0);
  accepted_ip_ = make_unique<TH1D>(
      MakeString("accepted_ip_", Counter::instance().counter()).c_str(),
      ";impact parameter[fm]", 100, b_min_, b_max_);
  accepted_ip_->SetDirectory(0);
}

void MCGlauber::clear() {
  if (nucleusA_.get() != nullptr)
    nucleusA_->clear();
  if (nucleusB_.get() != nullptr)
    nucleusB_->clear();
}

// run the glauber MC for N events
void MCGlauber::run(unsigned N) {
  initOutput();
  
  int errors = 0;
  int max_errors = 10;
  while (events_accepted_ < N) {
    EventStatus status = generate();
  
    if (status == EventStatus::Error) {
      SCT_ASSERT(++errors < max_errors,
                 "Repeated errors generating nuclei - aborting");
      continue;
    }

    events_generated_++;
    if (status == EventStatus::Hit) {
      events_accepted_++;
      accepted_ip_->Fill(tree_->B());
      if (events_accepted_ % 100 == 0) {
        VLOG(1) << "MCGlauber::run::event = " << events_accepted_;
        VLOG(1) << "(accepted/generated, Npart, Ncoll, b) = "
                << MakeString("(", (double)events_accepted_ / events_generated_,
                              ", ", tree_->nPart(), ", ", tree_->nColl(), ", ",
                              tree_->B(), ")");
      }
    }
  }
  writeHeader();
}

void MCGlauber::writeHeader() {
  // fill header
  tree_->setNameNucleusA(nucleusA_->name());
  tree_->setNameNucleusB(nucleusB_->name());
  tree_->setMassNumberA(nucleusA_->massNumber());
  tree_->setMassNumberB(nucleusB_->massNumber());
  tree_->setSigmaNN(collision_.NNCrossSection());
  tree_->setSqrtSNN(energy_);
  tree_->setRepulsionD(nucleusA_->repulsionDistance());
  tree_->setSmearHardCore(nucleusA_->nucleonSmearing() ==
                          NucleonSmearing::HardCore);
  tree_->setSmearGaussian(nucleusA_->nucleonSmearing() ==
                          NucleonSmearing::Gaussian);
  tree_->setCollisionHardCore(collision_.collisionProfile() ==
                              CollisionProfile::HardCore);
  tree_->setCollisionGaussian(collision_.collisionProfile() ==
                              CollisionProfile::Gaussian);
  tree_->setBMax(b_max_);
  tree_->setBMin(b_min_);
  tree_->setNEventsAccepted(events_accepted_);
  tree_->setNEventsThrown(events_generated_);

  // these parameters are all possible parameters depending on the form of the
  // NucleonPDF
  parameter_list param_a = nucleusA_->nuclearPDF().parameters();
  parameter_list param_b = nucleusB_->nuclearPDF().parameters();

  for (auto &par : param_a) {
    string par_name = par.first;
    double par_val = par.second;
    if (par_name == "radius") {
      tree_->setRadiusA(par_val);
    } else if (par_name == "skin_depth") {
      tree_->setSkinDepthA(par_val);
    } else if (par_name == "beta2") {
      tree_->setBeta2A(par_val);
    } else if (par_name == "beta4") {
      tree_->setBeta4A(par_val);
    } else if (par_name == "a") {
      tree_->setHulthenAA(par_val);
    } else if (par_name == "b") {
      tree_->setHulthenBA(par_val);
    } else {
      LOG(ERROR) << "MCGlauber/GlauberTree not aware of PDF parameter: "
                 << par_name << ", it will not be written to file";
    }
  }

  for (auto &par : param_b) {
    string par_name = par.first;
    double par_val = par.second;
    if (par_name == "radius") {
      tree_->setRadiusB(par_val);
    } else if (par_name == "skin_depth") {
      tree_->setSkinDepthB(par_val);
    } else if (par_name == "beta2") {
      tree_->setBeta2B(par_val);
    } else if (par_name == "beta4") {
      tree_->setBeta4B(par_val);
    } else if (par_name == "a") {
      tree_->setHulthenAB(par_val);
    } else if (par_name == "b") {
      tree_->setHulthenBB(par_val);
    } else {
      LOG(ERROR) << "MCGlauber/GlauberTree not aware of PDF parameter: "
                 << par_name << ", it will not be written to file";
    }
  }

  // calculate total cross section
  double x_min =
      events_accepted_ * pow(b_min_, 2) * pi * 10 / events_generated_;
  double x_max =
      events_accepted_ * pow(b_max_, 2) * pi * 10 / events_generated_;
  double x_total = x_max - x_min;
  double xErr =
      (x_max - x_min) * sqrt(1.0 / events_accepted_ + 1.0 / events_generated_);
  tree_->setTotalXsec(x_total);
  tree_->setTotalXsecError(xErr);

  // write header, and write tree to file
  tree_->fillHeader();
}

EventStatus MCGlauber::generate() {
  tree_->clearEvent();
  clear();

  // generate impact parameter
  double b = b_min_ + Random::instance().linear() * (b_max_ - b_min_);
  generated_ip_->Fill(b);
  
  // generate two new nuclei
  bool status_1 = nucleusA_->generate(b / 2.0);
  bool status_2 = nucleusB_->generate(-b / 2.0);

  if (!status_1 || !status_2)
    return EventStatus::Error;

  // collide the nuclei: if there are no collisions, we are done
  if (!collision_.collide(*nucleusA_.get(), *nucleusB_.get()))
    return EventStatus::Miss;

  // fill the event record in the tree
  tree_->setB(b);
  tree_->setNpart(collision_.nPart());
  tree_->setNcoll(collision_.nColl());
  tree_->setNspectators(collision_.spectators());
  tree_->setTheta(0, nucleusA_->nucleusTheta());
  tree_->setTheta(0, nucleusA_->nucleusPhi());
  tree_->setTheta(1, nucleusB_->nucleusTheta());
  tree_->setTheta(1, nucleusB_->nucleusPhi());

  for (auto &weight : glauberWeightSet) {
    int index = static_cast<int>(weight);
    tree_->setSumX(weight, collision_.averageX()[index]);
    tree_->setSumY(weight, collision_.averageY()[index]);
    tree_->setSumX2(weight, collision_.averageX2()[index]);
    tree_->setSumY2(weight, collision_.averageY2()[index]);
    tree_->setSumXY(weight, collision_.averageXY()[index]);
    tree_->setRP2Ecc(weight, collision_.reactionPlane2Ecc()[index]);
    tree_->setPP2Ecc(weight, collision_.partPlane2Ecc()[index]);
    tree_->setPP3Ecc(weight, collision_.partPlane3Ecc()[index]);
    tree_->setPP4Ecc(weight, collision_.partPlane4Ecc()[index]);
    tree_->setPP2(weight, collision_.partPlane2()[index]);
    tree_->setPP3(weight, collision_.partPlane3()[index]);
    tree_->setPP4(weight, collision_.partPlane4()[index]);
  }

  tree_->fill();
  return EventStatus::Hit;
}

double MCGlauber::lookupXSec(CollisionEnergy energy) {
  switch (energy) {
  case CollisionEnergy::E2760:
    return 6.4;
    break;
  case CollisionEnergy::E200:
    return 4.2;
    break;
  case CollisionEnergy::E62:
    return 3.6;
    break;
  case CollisionEnergy::E39:
    return 3.4;
    break;
  case CollisionEnergy::E27:
    return 3.3;
    break;
  case CollisionEnergy::E19:
    return 3.2;
    break;
  case CollisionEnergy::E14:
    return 3.15;
    break;
  case CollisionEnergy::E11:
    return 3.12;
    break;
  case CollisionEnergy::E7:
    return 3.08;
    break;
  default:
    LOG(ERROR)
        << "No inelastic nucleon-nucleon cross section for collision energy";
    LOG(ERROR) << "Either implement in lookupXSec(CollisionEnergy) or set "
                  "xsec by hand";
    return 0.0;
  }
}

void MCGlauber::init(GlauberSpecies species_A, GlauberSpecies species_B,
                     GlauberMod mod, bool deformed_a, bool deformed_b) {
  nucleusA_ = make_unique<Nucleus>(species_A, mod, deformed_a);
  nucleusB_ = make_unique<Nucleus>(species_B, mod, deformed_b);
}

void MCGlauber::init(unsigned mass_number_A, NucleonPDF::PDF pdf_A,
                     parameter_list params_A, unsigned mass_number_B,
                     NucleonPDF::PDF pdf_B, parameter_list params_B) {
  nucleusA_ = make_unique<Nucleus>(mass_number_A, params_A, pdf_A);
  nucleusB_ = make_unique<Nucleus>(mass_number_B, params_B, pdf_B);
}

} // namespace sct
