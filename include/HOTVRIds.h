#pragma once

#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"

/** \brief HOTVR top tagger
 * 
 * Cuts on the pt fraction of the leading subjet, the minimum pairwise
 * subjet mass, and on the top jet mass.  A cut on n_subjets >= 3 is
 * always applied.
 * 
 * The default values of the constructor correspond to the default values used in HOTVR.
 */
class HOTVRTopTag {
 public:
  explicit HOTVRTopTag(double fpt_upper = 0.8, double mjet_lower = 140., double mjet_upper = 220., double mmin_lower = 50.);
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;

 private:
  double m_fpt_upper;
  double m_mjet_lower;
  double m_mjet_upper;
  double m_mmin_lower;

};

class LepOverlap {
 public: 
  explicit LepOverlap(bool clacRadius = true);
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;

 private:
  bool m_calcRadius;

};

class DeltaPhiCut {
 public: 
  explicit DeltaPhiCut(double deltaphi_lower = M_PI/2);
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;

 private:
  double m_deltaphi_lower;

};

class Tau32Groomed {
 public: 
  explicit Tau32Groomed(double tau32_upper);
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;

 private:
  double m_tau32_upper;
};
  
class bTag {
 public:
  explicit bTag(JetId jetid);
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;

 private:
  JetId m_jetid;
};

class GenMatch {
 public:
  explicit GenMatch(uhh2::Context &ctx);
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;
 private:
  uhh2::Event::Handle<BstarToTWGen> h_BstarToTWGen;
};

class GenDeltaPhiCut {
 public: 
  explicit GenDeltaPhiCut(uhh2::Context &ctx, double deltaphi_lower = M_PI/2);
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;

 private:
  uhh2::Event::Handle<BstarToTWGen> h_BstarToTWGen;
  double m_deltaphi_lower;

};

class GenLepOverlap {
 public: 
  explicit GenLepOverlap(uhh2::Context &ctx, bool clacRadius = true);
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;

 private:
  bool m_calcRadius;
  uhh2::Event::Handle<BstarToTWGen> h_BstarToTWGen;

};
