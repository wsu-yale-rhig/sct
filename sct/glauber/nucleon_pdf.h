#ifndef SCT_GLAUBER_NUCLEON_PDF_H
#define SCT_GLAUBER_NUCLEON_PDF_H

/* Implements the nucleon PDF for generation of the nucleus. Provides
 * implementations for Woods-Saxon models, (symmetric or deformed), and a
 * functional form for deuteron (using the Hulthen form, see
 * sct/utils/functions.h). Proton distribution is a delta function.
 */

#include "sct/lib/enumerations.h"
#include "sct/lib/map.h"
#include "sct/lib/memory.h"
#include "sct/lib/string/string.h"
#include "sct/utils/functions.h"

#include <vector>

#include "TF1.h"
#include "TF2.h"

namespace sct {

typedef sct_map<string, double> parameter_list;

class NucleonPDF {
public:
  enum class PDF {
    Undetermined,
    WoodsSaxon1D,
    WoodsSaxon2D,
    StepFunction1D,
    Hulthen
  };

  // By default the nucleon pdf is uninitialized - user must call Init() before
  // sample() can be called
  NucleonPDF();
  ~NucleonPDF();

  // clears any stored distributions
  void clear();

  bool init(GlauberSpecies species, GlauberMod mod = GlauberMod::Nominal,
            bool deformed = false);
  bool init(PDF pdf, parameter_list &parameters);

  bool sample(double &r, double &theta, double &phi);

  // returns true if the theta distribution is non-uniform
  bool deformed() { return pdf_form_ == PDF::WoodsSaxon2D; }

  PDF PDFForm() { return pdf_form_; }

  void setParameter(string par_name, double val);

  sct_map<string, double> parameters() { return parameters_; }

private:
  bool initWS1D(parameter_list &params);
  bool initWS2D(parameter_list &params);
  bool initStepFunction(parameter_list &params);
  bool initHulthen(parameter_list &params);

  bool initParameters(parameter_list &params,
                      const parameter_template &par_template, TF1 *f);

  bool sample_1d(double &r, double &theta, double &phi);
  bool sample_2d(double &r, double &theta, double &phi);

  unique_ptr<TF1> pdf_1d_;
  unique_ptr<TF2> pdf_2d_;

  PDF pdf_form_;
  parameter_list parameters_;
};
} // namespace sct

#endif // SCT_GLAUBER_NUCLEON_PDF_H
