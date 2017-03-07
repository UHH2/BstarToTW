#include "UHH2/core/include/LorentzVector.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/BstarToTW/include/HOTVRIds.h"

using namespace uhh2;
using namespace std;

HOTVRTopTag::HOTVRTopTag(double fpt_upper, double mjet_lower, double mjet_upper, double mmin_lower):
  m_fpt_upper(fpt_upper), m_mjet_lower(mjet_lower), m_mjet_upper(mjet_upper), m_mmin_lower(mmin_lower) {}

bool HOTVRTopTag::operator()(const TopJet &topjet, const Event &event) const {
  vector<Jet> subjets = topjet.subjets();
  if(subjets.size() < 3) return false;

  double fpt = subjets.at(0).v4().pt() / topjet.v4().pt();
  if (fpt > m_fpt_upper) return false;

  double mjet = topjet.v4().M();
  if(mjet < m_mjet_lower) return false;
  if(mjet > m_mjet_upper) return false;

  double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
  double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
  double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
  double mmin = min(min(m12, m13), m23);

  if(mmin < m_mmin_lower) return false;

  return true;
}

DeltaPhiCut::DeltaPhiCut(double deltaphi_lower):
  m_deltaphi_lower(deltaphi_lower) {}

bool DeltaPhiCut::operator()(const TopJet &topjet, const Event &event) const {
  LorentzVector muon = event.muons->at(0).v4();
  return ( deltaPhi(muon, topjet.v4()) > m_deltaphi_lower );
}
