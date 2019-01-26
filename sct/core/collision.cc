// sct/core/collision.cc

#include "sct/core/collision.hh"

#include "sct/core/logging.hh"
#include "sct/utils/random.hh"

namespace sct {
  Collision::Collision()
  : inelasticXSec_(0), profile_(CollisionProfile::HardCore){
    
  }
  
  Collision::~Collision() { }

  void Collision::setMultiplicityModel(double npp, double k, double x, double ppEff,
                                       double aaEff, double aaCent, double trigEff,
                                       bool constEff) {
    mult_model_ = make_unique<MultiplicityModel>(npp, k, x, ppEff, aaEff, aaCent, trigEff, constEff);
  }
  
  template<typename Container>
  bool Collision::collide(Container& nucleusA, Container& nucleusB) {
    clear();
    
    // count the number of binary collisions
    for (int i = 0; i < nucleusA.size(); ++i) {
      for (int j = 0; j < nucleusB.size(); ++j) {
        if (nucleonCollision(nucleusA[i], nucleusB[j])) {
          nColl_++;
          nucleusA[i].incrementNColl();
          nucleusB[j].incrementNColl();
        }
      }
    }
    
    // if there are no binary collisions, return
    if (nColl_ == 0)
      return false;
    
    // get indices for the arrays
    int nPart_index = static_cast<int>(GlauberWeight::NPart);
    int nColl_index = static_cast<int>(GlauberWeight::NColl);
    int spectator_index = static_cast<int>(GlauberWeight::Spectators);
    int multiplicity_index = static_cast<int>(GlauberWeight::Multiplicity);
    
    // calculate kinematic event averages
    for (auto nucleus : {&nucleusA, &nucleusB}) {
      for (int i = 0; i < nucleus->size(); ++i) {
        Nucleon nucleon = (*nucleus)[i];
        
        if (nucleon.nColl() > 0) {
          // particle weighted averages
          count_[nPart_index]++;
          avgX_[nPart_index] += nucleon.x();
          avgY_[nPart_index] += nucleon.y();
          avgX2_[nPart_index] += nucleon.x2();
          avgY2_[nPart_index] += nucleon.y2();
          avgXY_[nPart_index] += nucleon.xy();
          
          // collision weighted averages
          count_[nColl_index] += nucleon.nColl();
          avgX_[nColl_index] += nucleon.nColl() * nucleon.x();
          avgY_[nColl_index] += nucleon.nColl() * nucleon.y();
          avgX2_[nColl_index] += nucleon.nColl() * nucleon.x2();
          avgY2_[nColl_index] += nucleon.nColl() * nucleon.y2();
          avgXY_[nColl_index] += nucleon.nColl() * nucleon.xy();

          if (mult_model_.get() != nullptr) {
            int mult = mult_model_->twoComponentMultiplicity(1.0, nucleon.nColl());
            nucleon.setMultiplicity(mult);
            count_[multiplicity_index] += mult;
            avgX_[multiplicity_index] += mult * nucleon.x();
            avgY_[multiplicity_index] += mult * nucleon.y();
            avgX2_[multiplicity_index] += mult * nucleon.x2();
            avgY2_[multiplicity_index] += mult * nucleon.y2();
            avgXY_[multiplicity_index] += mult * nucleon.xy();
          }
        }
        else {
          // spectator weighted averages
          count_[spectator_index]++;
          avgX_[spectator_index] += nucleon.x();
          avgY_[spectator_index] += nucleon.y();
          avgX2_[spectator_index] += nucleon.x2();
          avgY2_[spectator_index] += nucleon.y2();
          avgXY_[spectator_index] += nucleon.xy();
        }
      }
    }
    
    // now calculate eccentricity for each class
    for (int idx : {nPart_index, nColl_index, spectator_index, multiplicity_index}) {
      if (count_[idx] == 0)
        continue;
      
      if (idx == multiplicity_index && mult_model_ == nullptr)
        continue;
      
      avgX_[idx] /= count_[idx];
      avgY_[idx] /= count_[idx];
      avgX2_[idx] /= count_[idx];
      avgY2_[idx] /= count_[idx];
      avgXY_[idx] /= count_[idx];
      
      double sigmaX2 = avgX2_[idx] - pow(avgX_[idx], 2.0);
      double sigmaY2 = avgY2_[idx] - pow(avgY_[idx], 2.0);
      double sigmaSum2 = sigmaX2 + sigmaY2;
      
      // for the vanishingly small chance that the denominator is zero
      if (sigmaSum2 == 0) {
        eccRP2_[idx] = -999;
      }
      else {
        // reaction plane eccentricity
        eccRP2_[idx] = (sigmaY2 - sigmaX2) / sigmaSum2;
      }
      
      // calculate nth order participant plane eccentricity
      for (int order : {2, 3, 4}) {
        std::pair<double, double>
        pp_ecc = participantPlaneEcc(nucleusA, nucleusB, order, idx);
        if (order == 2) {
          pp2_[idx] = pp_ecc.first;
          eccPP2_[idx] = pp_ecc.second;
        }
        if (order == 3) {
          pp3_[idx] = pp_ecc.first;
          eccPP3_[idx] = pp_ecc.second;
        }
        if (order == 4) {
          pp4_[idx] = pp_ecc.first;
          eccPP4_[idx] = pp_ecc.second;
        }
      } // participant plane order
    } // weight index
    return true;
  }
  
  // instantiate for common containers
  template bool Collision::collide<vector<Nucleon>>(vector<Nucleon>& nucleusA,
                                                    vector<Nucleon>& nucleusB);
  template bool Collision::collide<Nucleus>(Nucleus& nucleusA, Nucleus& nucleusB);
  
  bool Collision::nucleonCollision(Nucleon nucleonA, Nucleon nucleonB) {
    double dR = nucleonA.deltaXY(nucleonB);
    double dRMax = sqrt(inelasticXSec_/pi);
    
    switch (profile_) {
      case CollisionProfile::HardCore :
        if (dR <= dRMax)
          return true;
        return false;
        break;
      case CollisionProfile::Gaussian :
        return Random::instance().uniform() <= exp(-pow(dR/dRMax, 2.0) / 2.0);
        break;
    }
  }

  void Collision::clear() {
    nColl_ = 0;
    count_.fill(0);
    avgX_.fill(0);
    avgY_.fill(0);
    avgX2_.fill(0);
    avgY2_.fill(0);
    avgXY_.fill(0);
    eccRP2_.fill(0);
    eccPP2_.fill(0);
    eccPP3_.fill(0);
    eccPP4_.fill(0);
    pp2_.fill(0);
    pp3_.fill(0);
    pp4_.fill(0);
  }
  
  template<typename Container>
  std::pair<double, double>
  Collision::participantPlaneEcc(Container& nucleusA, Container& nucleusB, int order,
                                 int weightIdx) {
    double qx = 0;
    double qy = 0;
    double qw = 0;
    
    if (nPart() <= 2)
      return {-999.0, -999.0};

    // get indices for the arrays
    constexpr unsigned nPart_index = static_cast<int>(GlauberWeight::NPart);
    constexpr unsigned nColl_index = static_cast<int>(GlauberWeight::NColl);
    constexpr unsigned spectator_index = static_cast<int>(GlauberWeight::Spectators);
    constexpr unsigned multiplicity_index = static_cast<int>(GlauberWeight::Multiplicity);
    
    for (auto nucleus : {&nucleusA, &nucleusB}) {
      for (int i = 0; i < nucleus->size(); ++i) {
        Nucleon nucleon = (*nucleus)[i];
        double nColl = nucleon.nColl();
        
        // spectators have zero collisions
        if (weightIdx == spectator_index && nColl > 0)
          continue;
        
        // participants have one or more collisions
        if (weightIdx != spectator_index && nColl == 0)
          continue;
        
        TVector3 tmp(nucleon.x() - avgX_[weightIdx],
                     nucleon.y() - avgY_[weightIdx],
                     nucleon.z());
        double weight = 1.0;
        if (weightIdx == nColl_index)
          weight = nColl;
        else if (weightIdx == multiplicity_index)
          weight = nucleon.multiplicity();
        double r = tmp.Perp();
        double phi = tmp.Phi();
        
        qx += weight * r * r * cos(order * phi);
        qy += weight * r * r * sin(order * phi);
        qw += weight * r * r;
        
      }
    }
    
    qx /= count_[weightIdx];
    qy /= count_[weightIdx];
    qw /= count_[weightIdx];
    
    double participant_plane = atan2(qy, -qx);
    double eccentricity = sqrt(qx * qx + qy * qy) / qw;
    
    return {participant_plane, eccentricity};
  }
} // namespace sct
