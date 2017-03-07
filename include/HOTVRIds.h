#pragma once

#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"


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

class DeltaPhiCut {
 public: 
  explicit DeltaPhiCut(double deltaphi_lower = M_PI/2);
  bool operator()(const TopJet &topjet, const uhh2::Event &event) const;

 private:
  double m_deltaphi_lower;
};
  
