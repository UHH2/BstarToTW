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

LepOverlap::LepOverlap(bool calcRadius):
m_calcRadius(calcRadius){}

bool LepOverlap::operator()(const TopJet &topjet, const Event &event) const {
  for (Muon muon : *event.muons)
    {
      double radius = 0.8;
      if (m_calcRadius)
	{
	  radius = 600./topjet.v4().pt();
	  if (radius < 0.1) radius = 0.1;
	  else if (radius > 1.5) radius = 1.5;
	}
      if (deltaR(muon.v4(), topjet.v4()) <= radius) return false;
    }
  return true;
}

DeltaPhiCut::DeltaPhiCut(double deltaphi_lower):
  m_deltaphi_lower(deltaphi_lower) {}

bool DeltaPhiCut::operator()(const TopJet &topjet, const Event &event) const {
  if(event.muons->size() < 1) return false;
  LorentzVector muon = event.muons->at(0).v4();
  return ( abs(deltaPhi(muon, topjet.v4())) > m_deltaphi_lower );
}

Tau32Groomed::Tau32Groomed(double tau32_upper):
  m_tau32_upper(tau32_upper) {}

bool Tau32Groomed::operator()(const TopJet &topjet, const Event &event) const {
  double tau32 = topjet.tau3_groomed() / topjet.tau2_groomed();
  return ( tau32 < m_tau32_upper );
}

bTag::bTag(JetId jetid):
  m_jetid(jetid) {}

bool bTag::operator()(const TopJet &topjet, const Event &event) const {

  double rho = 600.;
  double radius = rho/topjet.v4().pt();
  if (radius < 0.1) radius = 0.1;
  else if (radius > 1.5) radius = 1.5;

  vector<Jet> jets = *event.jets;
  for (Jet jet : jets)
    {
      if (m_jetid(jet, event))
	{
	  if (deltaR(topjet, jet) < radius) return true;
	}
    }
  return false;
}

GenMatch::GenMatch(Context &ctx) {
  h_BstarToTWGen = ctx.get_handle<BstarToTWGen>("BstarToTWgen");
}

bool GenMatch::operator()(const TopJet &topjet, const Event &event) const {
  BstarToTWGen gen = event.get(h_BstarToTWGen);
  const auto top = gen.tbstar();
  double rho = 600.;
  double radius = rho/topjet.v4().pt();
  if (radius < 0.1) radius = 0.1;
  else if (radius > 1.5) radius = 1.5;

  double pt_jet = topjet.pt();
  double pt_top = top.pt();
  double dPt = abs(pt_jet - pt_top) / pt_top;

  return (deltaR(top, topjet) <= radius && dPt < 0.05);

}

GenDeltaPhiCut::GenDeltaPhiCut(Context &ctx, double deltaphi_lower):
  m_deltaphi_lower(deltaphi_lower) {
  h_BstarToTWGen = ctx.get_handle<BstarToTWGen>("BstarToTWgen");
}

bool GenDeltaPhiCut::operator()(const TopJet &topjet, const Event &event) const {
  BstarToTWGen gen = event.get(h_BstarToTWGen);
  const auto lep = gen.ChargedLepton();
  return ( deltaPhi(lep, topjet.v4()) > m_deltaphi_lower );
}

GenLepOverlap::GenLepOverlap(Context &ctx, bool calcRadius):
  m_calcRadius(calcRadius){
  h_BstarToTWGen = ctx.get_handle<BstarToTWGen>("BstarToTWgen");
}

bool GenLepOverlap::operator()(const TopJet &topjet, const Event &event) const {
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
