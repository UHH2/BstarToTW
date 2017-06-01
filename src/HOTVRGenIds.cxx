#include "UHH2/core/include/LorentzVector.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/BstarToTW/include/HOTVRGenIds.h"

using namespace uhh2;
using namespace std;

GenHOTVRTopTag::GenHOTVRTopTag(double fpt_upper, double mjet_lower, double mjet_upper, double mmin_lower):
  m_fpt_upper(fpt_upper), m_mjet_lower(mjet_lower), m_mjet_upper(mjet_upper), m_mmin_lower(mmin_lower) {}

bool GenHOTVRTopTag::operator()(const GenTopJet &gentopjet, const Event &event) const {
  vector<Particle> subjets = gentopjet.subjets();
  if(subjets.size() < 3) return false;

  double fpt = subjets.at(0).v4().pt() / gentopjet.v4().pt();
  if (fpt > m_fpt_upper) return false;

  double mjet = gentopjet.v4().M();
  if(mjet < m_mjet_lower) return false;
  if(mjet > m_mjet_upper) return false;

  double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
  double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
  double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
  double mmin = min(min(m12, m13), m23);

  if(mmin < m_mmin_lower) return false;

  return true;
}


GenTopMuonCleaner::GenTopMuonCleaner(Context &ctx, bool calcRadius):
  m_calcRadius(calcRadius){
  h_BstarToTWGen = ctx.get_handle<BstarToTWGen>("BstarToTWgen");
}

bool GenTopMuonCleaner::operator()(const GenTopJet &topjet, const Event &event) const {
  BstarToTWGen gen = event.get(h_BstarToTWGen);
  const auto lep = gen.ChargedLepton();
  double radius = 0.8;
  if (m_calcRadius)
    {
      radius = 600./topjet.v4().pt();
      if (radius < 0.1) radius = 0.1;
      else if (radius > 1.5) radius = 1.5;
    }
  return ( deltaR(lep, topjet.v4()) > radius );
}
