#include "sct/utils/functions.h"

#include <cmath>

#include "sct/lib/math.h"

namespace sct {
    
  double StepFunction(double *x, double *par) {
      
    double x_ = x[0];
    double y_ = x[1];
    double z_ = x[2];
    double xsec = par[0];
      
    double distance = std::sqrt(std::pow(x_, 2) +
                                std::pow(y_, 2) +
                                std::pow(z_, 2));
      
    double r2 = xsec / pi;
    double r  = std::sqrt(r2);
    double V  = std::pow(r, 3.0) * pi * 4.0 / 3.0;
      
    return (r - distance < 0.0) ? 0.0 : 1 / V;
  }
    
  double Gaussian(double *x, double *par) {
      
    double x_ = x[0];
    double y_ = x[1];
    double z_ = x[2];
    double sigma = par[0];
      
    double r = std::sqrt(std::pow(x_, 2) +
                         std::pow(y_, 2) +
                         std::pow(z_, 2));
    double sigma2 = std::pow(sigma, 2.0);
    double norm = 2.0 * pi * sigma2;
      
    return 1.0 / std::pow(norm, 3.0 / 2.0) *
           std::exp( - 0.5 * std::pow(r, 2.0) / sigma2);
  }
    
  double WoodsSaxonSpherical(double *x, double *par) {
    double r = x[0];
    double R = par[0];
    double a = par[1];
    return std::pow(r, 2) / (1.0 + std::exp((r - R) / a));
  }
    
  double WoodsSaxonDeformed(double *x, double *par) {
    double r = x[0];
    double cosTheta = x[1];
    double R0 = par[0];
    double a = par[1];
    double beta2 = par[2];
    double beta4 = par[3];
      
    // powers of cosTheta
    double cos_theta2 = std::pow(cosTheta, 2.0);
    double cos_theta4 = std::pow(cosTheta, 4.0);
    
    // spherical harmonics Y(l=2, m=0), Y(l=4, m=0)
    double Y20 = std::sqrt(5.0 / pi) / 4.0 * (3.0 * cos_theta2 - 1.0);
    double Y40 = std::sqrt(1.0 / pi) * 3.0 / 16.0 * (35.0 * cos_theta4 - 30.0 * cos_theta2 + 3.0);
      
    double R = R0 * (1.0 + beta2 * Y20 + beta4 * Y40);
      
    return std::pow(r, 2) / (1.0 + std::exp((r - R) / a));
  }
} // namespace sct
