#pragma once

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/GenTopJet.h"
#include "UHH2/common/include/ObjectIdUtils.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"

/** \brief HOTVR top tagger
 * 
 * Cuts on the pt fraction of the leading subjet, the minimum pairwise
 * subjet mass, and on the top jet mass.  A cut on n_subjets >= 3 is
 * always applied.
 * 
 * The default values of the constructor correspond to the default values used in HOTVR.
 */

typedef std::function<bool (const GenTopJet &, const uhh2::Event &)> GenTopJetId;

class GenHOTVRTopTag {
 public:
  explicit GenHOTVRTopTag(double fpt_upper = 0.8, double mjet_lower = 140., double mjet_upper = 220., double mmin_lower = 50.);
  bool operator()(const GenTopJet &gentopjet, const uhh2::Event &event) const;

 private:
  double m_fpt_upper;
  double m_mjet_lower;
  double m_mjet_upper;
  double m_mmin_lower;

};

class GenTopMuonCleaner {
 public: 
  explicit GenTopMuonCleaner(uhh2::Context &ctx, bool clacRadius = true);
  bool operator()(const GenTopJet &topjet, const uhh2::Event &event) const;

 private:
  bool m_calcRadius;
  uhh2::Event::Handle<BstarToTWGen> h_BstarToTWGen;

};
