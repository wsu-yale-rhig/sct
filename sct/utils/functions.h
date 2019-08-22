#ifndef SCT_UTIL_FUNCTIONS_H
#define SCT_UTIL_FUNCTIONS_H

// Implementation of a few functions that are used in the MC glauber generator.
// This includes the functional forms of the Nucleon PDFS

#include "sct/lib/map.h"
#include "sct/lib/string/string.h"

namespace sct {

// we specify default values for all function parameters - these defaults are
// not physically meaningful, they are selected so that the function is in a
// valid state. These defaults, along with the relevant indices, are stored in a
// dictionary
typedef std::pair<unsigned, double> param_def;
typedef sct_map<string, param_def> parameter_template;

// given a position x, is that point less than a distance d away from the origin
// - if so, returns 1/d
//  x[0] = x
//  par[0] = d
const unsigned StepFunction1D_npar = 1;
const parameter_template StepFunction1D_params{{"d", {0, 1.0}}};
double StepFunction1D(double *x, double *par);

// given a position (dx, dy, dz), is that point in a sphere
// of area sqrt(xsec/pi)^3*4/3*pi, and if so, return 1/V
// x[0] = dx
// x[1] = dy
// x[2] = dz
// par[0] = sigma
const unsigned StepFunction_npar = 1;
const parameter_template StepFunction_params{{"sigma", {0, 1.0}}};
double StepFunction(double *x, double *par);

// evaluate a 1D gaussian in r = sqrt(x^2 + y^2 + z^2)
// x[0] = x
// x[1] = y
// x[2] = z
// par[0] = sigma
const unsigned Gaussian_npar = 1;
const parameter_template Gaussian_params{{"sigma", {0, 1.0}}};
double Gaussian(double *x, double *par);

// woods-saxon distribution for a spherical nuclei
// returns r^2 / (1 + exp((r-R)/a))
// x[0] = radial distance
// par[0] = nuclear radius
// par[1] = surface thickness
const unsigned WoodsSaxonSpherical_npar = 2;
const parameter_template WoodsSaxonSpherical_params{{"radius", {0, 5.0}},
                                                    {"skin_depth", {1, 0.5}}};
double WoodsSaxonSpherical(double *x, double *par);

// 2D (r, cos(theta)) woods-saxon distribution for deformed nuclei
// spherical harmonics Y(l=2, m=0) = Y20, Y(l=4, m=0) = Y40
// Y20 = sqrt(5.0 / pi) / 4.0 * (3.0 * cosTheta^2 - 1.0)
// Y40 = sqrt(1.0 / pi) * 3.0 / 16.0 * (35.0 * cosTheta^4 - 30.0 * cosTheta^2
// + 3.0) R = R0 * (1.0 + beta2 * Y20 + beta4 * Y40); returns r^2 / (1 +
// exp((r-R)/a)) x[0] = radial distance x[1] = cos(theta) theta = polar angle
// par[0] = nuclear radius (R0)
// par[1] = surface thickness (a)
// par[2] = 2nd order deformation parameter (beta2)
// par[3] = 4th order deformation parameter (beta4)
const unsigned WoodsSaxonDeformed_npar = 4;
const parameter_template WoodsSaxonDeformed_params{{"radius", {0, 5.0}},
                                                   {"skin_depth", {1, 0.5}},
                                                   {"beta2", {2, 0.1}},
                                                   {"beta4", {3, 0.1}}};
double WoodsSaxonDeformed(double *x, double *par);

// Hulthen form for deuteron PDF (L. Hulth ́en andM.Sagawara, Handbuch der Physik
// 39 1 (1957)), as implemented in the PHOBOS Glauber arXiv:1408.2549 [nucl-ex]
// x[0] = r
// par[0] = a
// par[1] = b
// functional form ((e^(-a*r') - e^(-b*r')) / r')^2 where r' is the distance
// between the neutron and proton: r' = 2*r
const unsigned HulthenPDF_npar = 2;
const parameter_template HulthenPDF_params{{"a", {0, 0.228}},
                                           {"b", {1, 1.177}}};
double HulthenPDF(double *x, double *par);

} // namespace sct

#endif // SCT_UTIL_FUNCTIONS_H
