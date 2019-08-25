#ifndef SCT_GLAUBER_NUCLEON_PDF_H
#define SCT_GLAUBER_NUCLEON_PDF_H

/* Implements the nucleon PDF for generation of the nucleus. Provides
 * implementations for Woods-Saxon models, (symmetric or deformed), and a
 * functional form for deuteron (using the Hulthen form, see
 * sct/utils/functions.h). Proton distribution is a delta function.
 *
 * Available PDFs:
 *
 * Woods Saxon 1D: a mean field potential for nucleons when the mass number is
 * "large enough". In the 1D case, theta and phi are sampled from uniform
 * distributions. See sct/utils/functions.h for details
 * parameters:
 * radius - the width of the distribution
 * skin_depth - controls how fast the potential "dies off" at the edge of the
 * nucleus
 *
 * Woods Saxon 2D: same as 1D, but this time the potential is a function of
 * radius and theta - this is used for non-spherical nuclei. Phi is still
 * sampled uniformly. See sct/utils/functions.h for details
 * parameters:
 * radius - the width of the distribution
 * skin_depth - controls how fast the potential "dies off" at the edge of the
 * nucleus
 * beta2 - controls the magnitude of the second order spherical harmonic
 * contribution  (modulates theta)
 * beta4 - controls the magnitude of the fourth order spherical harmonic
 * contribution (modulates theta)
 *
 * Step function 1D: a step function in r, used to sample position for a single
 * proton. The width of the function is taken to be small, generally - almost
 * delta function-like. Theta and phi are sampled uniformly
 * parameters:
 * radius - controls the width of the step function.
 *
 * Hulthen PDF: functional form for deuteron. Generally sampled once, with theta
 * and phi being selected at random, and the second nucleon being placed
 * opposite.
 * parameters:
 * a - controls the first exponential width
 * b - controls the second exponential width
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
    Undetermined,   // parameters: none
    WoodsSaxon1D,   // parameters: "radius", "skin_depth"
    WoodsSaxon2D,   // parameters: "radius", "skin_depth", "beta2", "beta4"
    StepFunction1D, // parameters: "d"
    Hulthen         // parameters: "a", "b"
  };

  // By default the nucleon pdf is uninitialized - user must call Init() before
  // sample() can be called
  NucleonPDF();
  ~NucleonPDF();

  // clears any stored distributions
  void clear();

  bool init(GlauberSpecies species, GlauberMod mod = GlauberMod::Nominal,
            bool deformed = false);
  bool init(PDF pdf, parameter_list parameters);

  bool sample(double &r, double &theta, double &phi);

  // returns true if the theta distribution is non-uniform
  bool deformed() { return pdf_form_ == PDF::WoodsSaxon2D; }

  PDF PDFForm() const { return pdf_form_; }

  void setParameter(string par_name, double val);

  sct_map<string, double> parameters() const { return parameters_; }

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
