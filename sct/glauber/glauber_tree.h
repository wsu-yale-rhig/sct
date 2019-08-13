#ifndef SCT_GLAUBER_GLAUBER_TREE_H
#define SCT_GLAUBER_GLAUBER_TREE_H

#include "sct/lib/enumerations.h"
#include "sct/lib/memory.h"

#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"

namespace sct {
class GlauberTree {
 public:
  static const unsigned n_nuclei = 2;  // 2 nucleus collision (A, B)
  enum class TreeMode { Read, Write };

  // if filename is empty, assumes user wants to create a new tree
  GlauberTree(TreeMode mode, const string& filename = "");
  virtual ~GlauberTree();

  // zeros values for the event tree
  void clearEvent();

  // to read a Tree from file, use open(filename)
  bool open(const string& filename);

  // convenience function for the Tree to read from a file
  // that is opened elsewhere
  bool read(shared_ptr<TFile> file);

  // writes current values to event tree
  void fill();

  // writes current values to header tree
  void fillHeader();

  // write to current TFile
  void write();

  unsigned getEntries() const;
  bool getEntry(unsigned int idx);
  bool getHeaderEntry(unsigned int idx = 0);

  // event-wise variables
  inline void setB(double val) { b_ = val; }
  inline void setNpart(unsigned val) { nPart_ = val; }
  inline void setNcoll(unsigned val) { nColl_ = val; }
  inline void setNspectators(unsigned val) { nSpec_ = val; }
  inline void setMultiplicity(unsigned val) { multiplicity_ = val; }
  inline void setTheta(unsigned id, double val) { theta_[id] = val; }
  inline void setPhi(unsigned id, double val) { phi_[id] = val; }
  inline void setSumX(GlauberWeight id, double val) {
    sumX_[static_cast<unsigned>(id)] = val;
  }
  inline void setSumY(GlauberWeight id, double val) {
    sumY_[static_cast<unsigned>(id)] = val;
  }
  inline void setSumX2(GlauberWeight id, double val) {
    sumX2_[static_cast<unsigned>(id)] = val;
  }
  inline void setSumY2(GlauberWeight id, double val) {
    sumY2_[static_cast<unsigned>(id)] = val;
  }
  inline void setSumXY(GlauberWeight id, double val) {
    sumXY_[static_cast<unsigned>(id)] = val;
  }
  inline void setRP2Ecc(GlauberWeight id, double val) {
    rp2ecc_[static_cast<unsigned>(id)] = val;
  }
  inline void setPP2Ecc(GlauberWeight id, double val) {
    pp2ecc_[static_cast<unsigned>(id)] = val;
  }
  inline void setPP3Ecc(GlauberWeight id, double val) {
    pp3ecc_[static_cast<unsigned>(id)] = val;
  }
  inline void setPP4Ecc(GlauberWeight id, double val) {
    pp4ecc_[static_cast<unsigned>(id)] = val;
  }
  inline void setPP2(GlauberWeight id, double val) {
    pp2_[static_cast<unsigned>(id)] = val;
  }
  inline void setPP3(GlauberWeight id, double val) {
    pp3_[static_cast<unsigned>(id)] = val;
  }
  inline void setPP4(GlauberWeight id, double val) {
    pp4_[static_cast<unsigned>(id)] = val;
  }

  inline double B() const { return b_; }
  inline unsigned nPart() const { return nPart_; }
  inline unsigned nColl() const { return nColl_; }
  inline unsigned nSpectators() const { return nSpec_; }
  inline unsigned multiplicity() const { return multiplicity_; }
  inline double theta(unsigned id) const { return theta_[id]; }
  inline double phi(unsigned id) const { return phi_[id]; }
  inline double sumX(GlauberWeight id) const {
    return sumX_[static_cast<unsigned>(id)];
  }
  inline double sumY(GlauberWeight id) const {
    return sumY_[static_cast<unsigned>(id)];
  }
  inline double sumX2(GlauberWeight id) const {
    return sumX2_[static_cast<unsigned>(id)];
  }
  inline double sumY2(GlauberWeight id) const {
    return sumY2_[static_cast<unsigned>(id)];
  }
  inline double sumXY(GlauberWeight id) const {
    return sumXY_[static_cast<unsigned>(id)];
  }
  inline double RP2Ecc(GlauberWeight id) const {
    return rp2ecc_[static_cast<unsigned>(id)];
  }
  inline double PP2Ecc(GlauberWeight id) const {
    return pp2ecc_[static_cast<unsigned>(id)];
  }
  inline double PP3Ecc(GlauberWeight id) const {
    return pp3ecc_[static_cast<unsigned>(id)];
  }
  inline double PP4Ecc(GlauberWeight id) const {
    return pp4ecc_[static_cast<unsigned>(id)];
  }
  inline double PP2(GlauberWeight id) const {
    return pp2_[static_cast<unsigned>(id)];
  }
  inline double PP3(GlauberWeight id) const {
    return pp3_[static_cast<unsigned>(id)];
  }
  inline double PP4(GlauberWeight id) const {
    return pp4_[static_cast<unsigned>(id)];
  }

  // area calculations
  double RPArea(GlauberWeight id);
  double PPArea(GlauberWeight id);

  // header
  inline void setNameNucleusA(const string& val) { *nameA_ = val; }
  inline void setNameNucleusB(const string& val) { *nameB_ = val; }
  inline void setMassNumberA(unsigned val) { massNumberA_ = val; }
  inline void setMassNumberB(unsigned val) { massNumberB_ = val; }
  inline void setRadiusA(double val) { radiusA_ = val; }
  inline void setRadiusB(double val) { radiusB_ = val; }
  inline void setSkinDepthA(double val) { skinDepthA_ = val; }
  inline void setSkinDepthB(double val) { skinDepthB_ = val; }
  inline void setBeta2A(double val) { beta2A_ = val; }
  inline void setBeta2B(double val) { beta2B_ = val; }
  inline void setBeta4A(double val) { beta4A_ = val; }
  inline void setBeta4B(double val) { beta4B_ = val; }
  inline void setSigmaNN(double val) { sigmaNN_ = val; }
  inline void setSqrtSNN(double val) { sqrtSNN_ = val; }
  inline void setRepulsionD(double val) { repulsionD_ = val; }
  inline void setTotalXsec(double val) { totalXsec_ = val; }
  inline void setTotalXsecError(double val) { totalXsecError_ = val; }
  inline void setSmearHardCore(bool val) { smearHardCore_ = val; }
  inline void setSmearGaussian(bool val) { smearGaussian_ = val; }
  inline void setCollisionHardCore(bool val) { collisionHardCore_ = val; }
  inline void setCollisionGaussian(bool val) { collisionGaussian_ = val; }
  inline void setBMax(double val) { bMax_ = val; }
  inline void setBMin(double val) { bMin_ = val; }
  inline void setNEventsAccepted(unsigned val) { eventsAccepted_ = val; }
  inline void setNEventsThrown(unsigned val) { eventsThrown_ = val; }

  inline string nameNucleusA() const { return *nameA_; }
  inline string nameNucleusB() const { return *nameB_; }
  inline unsigned massNumberA() const { return massNumberA_; }
  inline unsigned massNumberB() const { return massNumberB_; }
  inline double radiusA() const { return radiusA_; }
  inline double radiusB() const { return radiusB_; }
  inline double skinDepthA() const { return skinDepthA_; }
  inline double skinDepthB() const { return skinDepthB_; }
  inline double beta2A() const { return beta2A_; }
  inline double beta2B() const { return beta2B_; }
  inline double beta4A() const { return beta4A_; }
  inline double beta4B() const { return beta4B_; }
  inline double sigmaNN() const { return sigmaNN_; }
  inline double sqrtSNN() const { return sqrtSNN_; }
  inline double repulsionD() const { return repulsionD_; }
  inline double totalXsec() const { return totalXsec_; }
  inline double totalXsecError() const { return totalXsecError_; }
  inline bool smearHardCore() const { return smearHardCore_; }
  inline bool smearGaussian() const { return smearGaussian_; }
  inline bool collisionHardCore() const { return collisionHardCore_; }
  inline bool collisionGaussian() const { return collisionGaussian_; }
  inline double BMax() const { return bMax_; }
  inline double BMin() const { return bMin_; }
  inline unsigned NEventsAccepted() const { return eventsAccepted_; }
  inline unsigned NEventsThrown() const { return eventsThrown_; }

 private:
  void loadBranches();
  void createBranches();

  shared_ptr<TFile> file_;
  TTree* header_;
  TTree* event_;

  double b_;
  unsigned nPart_;
  unsigned nColl_;
  unsigned nSpec_;
  unsigned multiplicity_;
  double theta_[n_nuclei];
  double phi_[n_nuclei];
  double sumX_[nGlauberWeights];
  double sumY_[nGlauberWeights];
  double sumX2_[nGlauberWeights];
  double sumY2_[nGlauberWeights];
  double sumXY_[nGlauberWeights];
  double rp2ecc_[nGlauberWeights];
  double pp2ecc_[nGlauberWeights];
  double pp3ecc_[nGlauberWeights];
  double pp4ecc_[nGlauberWeights];
  double pp2_[nGlauberWeights];
  double pp3_[nGlauberWeights];
  double pp4_[nGlauberWeights];

  string* nameA_;
  string* nameB_;
  unsigned massNumberA_;
  unsigned massNumberB_;
  double radiusA_;
  double radiusB_;
  double skinDepthA_;
  double skinDepthB_;
  double beta2A_;
  double beta2B_;
  double beta4A_;
  double beta4B_;
  double sigmaNN_;
  double sqrtSNN_;
  double repulsionD_;
  double totalXsec_;
  double totalXsecError_;
  bool smearHardCore_;
  bool smearGaussian_;
  bool collisionHardCore_;
  bool collisionGaussian_;
  double bMax_;
  double bMin_;
  unsigned eventsAccepted_;
  unsigned eventsThrown_;
};
}  // namespace sct

#endif  // SCT_GLAUBER_GLAUBER_TREE_H
