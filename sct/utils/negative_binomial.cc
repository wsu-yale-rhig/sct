#include "sct/utils/negative_binomial.h"

#include <limits>

#include "sct/utils/random.h"
#include "sct/lib/logging.h"
#include "sct/lib/string/string_utils.h"

#include "TMath.h"

namespace sct {

  NegativeBinomial::NegativeBinomial(double npp, double k)
  : npp_(npp), k_(k) {
    initNBD();
  }

  NegativeBinomial::~NegativeBinomial() {}
  
  double NegativeBinomial::evaluateNBD(int i, double m) const {
    double term_1 = TMath::Exp(TMath::LnGamma(i + k_ * m) -
                               TMath::LnGamma(i + 1) -
                               TMath::LnGamma(k_ * m));
    double term_2 = i * TMath::Log(npp_ / k_) - ( i + k_ * m ) *
                    TMath::Log(npp_ / k_ + 1.0);
    
    
    return term_1 * TMath::Exp(term_2);
  }
  
  void NegativeBinomial::setParameters(double npp, double k) {
    npp_ = npp;
    k_ = k;
    initNBD();
  }
  
  void NegativeBinomial::initNBD() {
    nbd_.release();
    
    int nBins = 100;
    nbd_ = make_unique<TH1D>(MakeString("nbd",
                             Random::instance().counter()).c_str(),
                             "", nBins, 0, nBins);
    nbd_->SetDirectory(0);
    for (int i = 0; i < nBins; ++i) {
      nbd_->SetBinContent(i+1, evaluateNBD(i));
    }
    nbd_->Scale(1.0 / nbd_->Integral());
  }
} // namespace sct
