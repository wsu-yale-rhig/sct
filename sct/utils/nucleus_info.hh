// sct/utils/nucleus_info.hh

#ifndef SCT_UTILS_NUCLEUS_INFO_HH
#define SCT_UTILS_NUCLEUS_INFO_HH

/* Table of nucleus parameters -
 * mass number
 * charge number
 * charge radius
 * skin depth
 * charge radius absolute error (2*sigma)
 * skin depth absolute error (2*sigma)
 * beta2 (second order deformation parameter)
 * beta4 (fourth order deformation parameter)
 *
 * SOURCES:
 * PHOBOS MC Glauber:
 * https://arxiv.org/abs/0805.4411
 * https://arxiv.org/abs/1408.2549
 * charge radius & skin depth:
 * H. De Vries, C.W. De Jager, and C. De Vries, Atomic Data and Nuclear Data Tables 36 495 (1987).
 * http://faculty.virginia.edu/ncd/dldata/1987_source.pdf
 * deformation parameters (beta2 & beta4)
 * P. Moller, J.R. Nix, W.D. Myers and W.J. Swiatecki, Atomic Data and Nuclear Data Tables 59 185 (1995).
 * https://www.sciencedirect.com/science/article/pii/S0092640X85710029
 * *NEW TABLE* https://arxiv.org/pdf/1508.06294.pdf need to check it out
 */

#include "sct/core/base.hh"
#include "sct/core/enumerations.hh"

namespace sct {
  class NucleusInfo {
  public:
    
    static NucleusInfo& instance();
    virtual ~NucleusInfo();
    
    unsigned massNumber(GlauberSpecies sp) {return mass_number_[sp];}
    unsigned chargeNumber(GlauberSpecies sp) {return charge_number_[sp];}
    double radius(GlauberSpecies sp) {return radius_[sp];}
    double skinDepth(GlauberSpecies sp) {return skin_depth_[sp];}
    double radiusError(GlauberSpecies sp) {return radius_error_[sp];}
    double skinDepthError(GlauberSpecies sp) {return skin_depth_error_[sp];}
    double beta2(GlauberSpecies sp) {return beta2_[sp];}
    double beta4(GlauberSpecies sp) {return beta4_[sp];}
    string name(GlauberSpecies sp) {return name_[sp];}
    
  private:
    
    std::unordered_map<GlauberSpecies, unsigned, EnumClassHash> mass_number_;
    std::unordered_map<GlauberSpecies, unsigned, EnumClassHash> charge_number_;
    std::unordered_map<GlauberSpecies, double, EnumClassHash> radius_;
    std::unordered_map<GlauberSpecies, double, EnumClassHash> skin_depth_;
    std::unordered_map<GlauberSpecies, double, EnumClassHash> radius_error_;
    std::unordered_map<GlauberSpecies, double, EnumClassHash> skin_depth_error_;
    std::unordered_map<GlauberSpecies, double, EnumClassHash> beta2_;
    std::unordered_map<GlauberSpecies, double, EnumClassHash> beta4_;
    std::unordered_map<GlauberSpecies, string, EnumClassHash> name_;
    
    NucleusInfo();
    NucleusInfo(const NucleusInfo&);
    
  };
} // namespace sct

#endif // SCT_UTILS_NUCLEUS_INFO_HH
