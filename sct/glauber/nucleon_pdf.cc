#include "sct/glauber/nucleon_pdf.h"
#include "sct/lib/assert.h"
#include "sct/lib/logging.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/functions.h"
#include "sct/utils/nucleus_info.h"
#include "sct/utils/random.h"
#include "sct/lib/math.h"

namespace sct {

NucleonPDF::NucleonPDF() : pdf_1d_(nullptr), pdf_2d_(nullptr) {}

NucleonPDF::~NucleonPDF() {}

void NucleonPDF::clear() {
  pdf_form_ = NucleonPDF::PDF::Undetermined;
  pdf_1d_.reset();
  pdf_2d_.reset();
  parameters_.clear();
}

bool NucleonPDF::init(GlauberSpecies species, GlauberMod mod, bool deformed) {
  clear();

  double radius = NucleusInfo::instance().radius(species);
  double skin_depth = NucleusInfo::instance().skinDepth(species);
  switch (mod) {
  case GlauberMod::Large:
    radius += NucleusInfo::instance().radiusError(species);
    skin_depth -= NucleusInfo::instance().skinDepthError(species);
    break;
  case GlauberMod::Small:
    radius -= NucleusInfo::instance().radiusError(species);
    skin_depth += NucleusInfo::instance().skinDepthError(species);
    break;
  default:
    break;
  }
  double beta2 = NucleusInfo::instance().beta2(species);
  double beta4 = NucleusInfo::instance().beta4(species);
  double hulthen_a = NucleusInfo::instance().hulthenA(species);
  double hulthen_b = NucleusInfo::instance().hulthenB(species);

  parameter_list params;
  switch (species) {
  case GlauberSpecies::p1:
    pdf_form_ = PDF::StepFunction1D;
    params = {{"radius", radius}};
    initStepFunction(params);
    break;
  case GlauberSpecies::d2:
    pdf_form_ = PDF::Hulthen;
    params = {{"a", hulthen_a}, {"b", hulthen_b}};
    initHulthen(params);
    break;
  default:
    if (deformed) {
      pdf_form_ = PDF::WoodsSaxon2D;
      params = {{"radius", radius},
                {"skin_depth", skin_depth},
                {"beta2", beta2},
                {"beta4", beta4}};
      initWS2D(params);
    } else {
      pdf_form_ = PDF::WoodsSaxon1D;
      params = {{"radius", radius}, {"skin_depth", skin_depth}};
      initWS1D(params);
    }
    break;
  }
  return true;
}

bool NucleonPDF::init(PDF pdf, parameter_list parameters) {
  clear();

  if (pdf == PDF::Undetermined) {
    SCT_THROW("Can not select Undetermined for initializing the PDF");
  }
  pdf_form_ = pdf;

  switch (pdf) {
  case PDF::StepFunction1D:
    initStepFunction(parameters);
    break;
  case PDF::Hulthen:
    initHulthen(parameters);
    break;
  case PDF::WoodsSaxon1D:
    initWS1D(parameters);
    break;
  case PDF::WoodsSaxon2D:
    initWS2D(parameters);
    break;
  default:
    SCT_THROW("Should not be here: this is an error in the SCT implementation");
  }

  return true;
}

bool NucleonPDF::sample(double &r, double &theta, double &phi) {
  
  if (pdf_2d_.get() != nullptr) {
    return sample_2d(r, theta, phi);
  } else if (pdf_1d_.get() != nullptr) {
    return sample_1d(r, theta, phi);
  } else {
    SCT_THROW("NucleonPDF is not initialized, can not sample");
  }
}

void NucleonPDF::setParameter(string par_name, double val) {
  if (parameters_.count(par_name) > 0) {
    if (deformed()) {
      pdf_2d_->SetParameter(par_name.c_str(), val);
    } else {
      pdf_1d_->SetParameter(par_name.c_str(), val);
    }
    parameters_[par_name] = val;
  }
}

bool NucleonPDF::initWS1D(parameter_list &params) {
  pdf_1d_ = make_unique<TF1>(
      sct::MakeString("woodsaxon_", Random::instance().counter()).c_str(),
      WoodsSaxonSpherical, 0.0, 20.0, WoodsSaxonSpherical_npar);
  pdf_1d_->SetNpx(400);
  initParameters(params, WoodsSaxonSpherical_params, pdf_1d_.get());
  return true;
}
bool NucleonPDF::initWS2D(parameter_list &params) {
  pdf_2d_ = make_unique<TF2>(
      sct::MakeString("woodsaxon_", Random::instance().counter()).c_str(),
      WoodsSaxonDeformed, 0.0, 20.0, -1.0, 1.0, WoodsSaxonDeformed_npar);
  
  pdf_2d_->SetNpx(400);
  pdf_2d_->SetNpy(400);
  
  initParameters(params, WoodsSaxonDeformed_params, (TF1 *)pdf_2d_.get());
  
  return true;
}
bool NucleonPDF::initStepFunction(parameter_list &params) {
  pdf_1d_ = make_unique<TF1>(
      sct::MakeString("stepfunction_", Random::instance().counter()).c_str(),
      StepFunction1D, 0.0, 20.0, StepFunction1D_npar);
  pdf_1d_->SetNpx(10000);
  initParameters(params, StepFunction1D_params, pdf_1d_.get());
  return true;
}
bool NucleonPDF::initHulthen(parameter_list &params) {
  pdf_1d_ = make_unique<TF1>(
      sct::MakeString("hulthen_", Random::instance().counter()).c_str(),
      HulthenPDF, 0.0, 20.0, HulthenPDF_npar);
  pdf_1d_->SetNpx(500);
  initParameters(params, HulthenPDF_params, pdf_1d_.get());
  return true;
}

bool NucleonPDF::sample_1d(double &r, double &theta, double &phi) {
  r = pdf_1d_->GetRandom();
  theta = acos(Random::instance().centeredUniform());
  phi = Random::instance().zeroToPi() * 2.0 - pi;
  return true;
}

bool NucleonPDF::sample_2d(double &r, double &theta, double &phi) {
  double cos_theta = 0.0;
  pdf_2d_.get()->GetRandom2(r, cos_theta);
  theta = acos(cos_theta);
  phi = Random::instance().zeroToPi() * 2.0 - pi;
  return true;
}

bool NucleonPDF::initParameters(parameter_list &params,
                                const parameter_template &par_template,
                                TF1 *f) {
  for (auto &par_def : par_template) {
    string par_name = par_def.first;
    unsigned par_idx = par_def.second.first;
    double par_val = par_def.second.second;
    try {
      par_val = params.at(par_name);
    } catch (std::exception &e) {
      LOG(ERROR) << "Could not find a parameter value for " << par_name
                 << " using default value: " << par_val;
    }
    f->SetParName(par_idx, par_name.c_str());
    f->SetParameter(par_name.c_str(), par_val);
    parameters_[par_name] = par_val;
  }
  return true;
}

} // namespace sct