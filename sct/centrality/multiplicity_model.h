#ifndef SCT_CENTRALITY_MULTIPLICITY_MODEL_H
#define SCT_CENTRALITY_MULTIPLICITY_MODEL_H

/*  Implements the two-component multilpicity model used in the
 *  glauber MC. multiplicity per collision is sampled using a
 *  simple ((1 - x_) * npart / 2.0) + (x_ * ncoll) relationship,
 *  where x is described as the fraction coming from "hard processes",
 *  (1-x) the fraction coming from other processes. X, however, stays
 *  approximately constant over a wide range of collision energies, so
 *  it probably has nothing to do with "hard" and "soft" processes, but
 *  it works well in the model.
 */

#include "sct/core/base.h"
#include "sct/utils/negative_binomial.h"

namespace sct {
  class MultiplicityModel : public NegativeBinomial {
  public:
    
    // creates a negative binomial distribution with parameters (npp, k), and
    // sets efficiencies & the trigger bias.
    MultiplicityModel(double npp = 2.38, double k = 2.00, double x = 0.13,
                      double ppEff = 0.98, double centEff = 0.84,
                      double centMult = 540, double triggerBias = 1.0,
                      bool const_eff = false);
    MultiplicityModel(const MultiplicityModel& rhs);
    virtual ~MultiplicityModel();
    
    // returns two component multiplicity (1-x) * npart/2 + (x) * ncoll
    double twoComponentMultiplicity(double npart, double ncoll) const;
    // returns convolution with NBD
    double multiplicity(double npart, double ncoll) const;
    // return multiplicity distribution with scaled NBD with mult*npp, k*mult
    TH1D* multiplicity(double npart, double ncoll, double weight) const;
    
    // reset NBD probability distribution
    void setNBD(double npp, double k);
    
    // depending on if ConstEfficiency() is true or false, this will return
    // different values - either the efficiency at 0-5% centrality if true, or
    // a parameterized efficiency that runs from pp_efficiency_ at 0 mult to
    // central_efficiency_ at 0-5% central mult
    double evalEfficiency(unsigned mult) const;
    
    // sets the fraction of multiplicity from "hard processes" in the two
    // component multiplicity model
    inline void setX(double val) {x_ = val;}
    inline double x() {return x_;}
    
    // multiplicity calculations take into account average detector efficiency
    // at the given "true" multiplicity, this sets the zero multiplicity limit
    // and should be set to pp or peripheral AuAu
    inline void setppEfficiency(double val) {pp_efficiency_ = val;}
    inline double ppEfficiency() {return pp_efficiency_;}
    
    // multiplicity calculations take into account average detector efficiency
    // at the given "true" multiplicity, this sets the high multiplicity limit
    // and should be set to the average efficiency in 0-5%
    inline  void setCentralEfficiency(double val) {central_efficiency_ = val;}
    inline double centralEfficiency() {return central_efficiency_;}
    
    // the average central (0-5%) multiplicity should be specified to define
    // the slope of the linear efficiency curve
    inline void setCentralMultiplicity(unsigned mult) {cent_mult_ = mult;}
    inline unsigned centralMultiplicity() const {return cent_mult_;}
    
    // if set to true, always use central_efficiency_, instead of a multiplicity
    // dependent parameterization
    inline void setConstEfficiency(bool flag) {const_efficiency_ = flag;}
    inline double constEfficiency() {return const_efficiency_;}
    
    // triggers can bias you towards higher multiplicities, this attempts
    // to correct that - if set to a value between [0, 1), multiplicity calculation
    // will add extra tracks with this value as the probability, so you end up with
    // between uncorrected multiplicity & 2 * uncorrected multiplicity, with an
    // average of (1 + trigger_bias_) * multiplicity
    inline void setTriggerBias(double val) {trigger_bias_ = val;}
    inline double triggerBias() {return trigger_bias_;}
    
  private:
    
    double pp_efficiency_;         // pp_efficiency
    double central_efficiency_;    // efficiency at 0-5% centrality for the given
                                   // collision system
    double cent_mult_;             // average multiplicity in 0-5% centrality
    double trigger_bias_;          // trigger bias
    double x_;                     // fraction of production from hard processes
    bool const_efficiency_;        // flag for constant or non-constant efficiency
  };
} // namespace sct

#endif // SCT_CENTRALITY_MULTIPLICITY_MODEL_H
