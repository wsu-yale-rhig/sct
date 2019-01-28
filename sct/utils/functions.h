#ifndef SCT_UTIL_FUNCTIONS_H
#define SCT_UTIL_FUNCTIONS_H

// Implementation of a few functions that are used in the MC glauber
// generator.

namespace sct {

// given a position (dx, dy, dz), is that point in a sphere
// of area sqrt(xsec/pi)^3*4/3*pi, and if so, return 1/V
// x[0] = dx
// x[1] = dy
// x[2] = dz
// par[0] = sigma
double StepFunction(double *x, double *par);

// evaluate a 1D gaussian in r = sqrt(x^2 + y^2 + z^2)
// x[0] = x
// x[1] = y
// x[2] = z
// par[0] = sigma
double Gaussian(double *x, double *par);

// woods-saxon distribution for a spherical nuclei
// returns r^2 / (1 + exp((r-R)/a))
// x[0] = radial distance
// par[0] = nuclear radius
// par[1] = surface thickness
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
double WoodsSaxonDeformed(double *x, double *par);

}  // namespace sct

#endif  // SCT_UTIL_FUNCTIONS_H
