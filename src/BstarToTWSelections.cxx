#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/common/include/Utils.h"

#include <stdexcept>
#include <vector>

using namespace uhh2;
using namespace std;

NHotvrSelection::NHotvrSelection(unsigned int n_min_, unsigned int n_max_, double pt_min_, double eta_max_) {
  n_min   = n_min_;
  n_max   = n_max_;
  pt_min  = pt_min_;
  eta_max = eta_max_;
}

bool NHotvrSelection::passes(const Event &event) {
  vector<TopJet> topjets = *event.topjets;
  unsigned int n = 0;
  for (TopJet topjet : topjets)
    { 
      double pt = topjet.v4().pt();
      double eta = abs(topjet.v4().eta());
      if ( (pt_min < pt) && (eta < eta_max) ) ++n;
    }
  return (n_min <= n && n <= n_max);
}

METSelection::METSelection(double met_min_) {
  met_min = met_min_;
}

bool METSelection::passes(const Event &event) {
  return (met_min < event.met->pt());
}

STSelection::STSelection(double st_min):
  m_st_min(st_min) {}

bool STSelection::passes(const Event &event) {
  auto met = event.met->pt();

  double ht = 0.0;
  double ht_jets = 0.0;
  double ht_lep = 0.0;
  for(const auto & jet : *event.jets){
    ht_jets += jet.pt();
  }
  for(const auto & electron : *event.electrons){
    ht_lep += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_lep += muon.pt();
  }

  ht = ht_lep + ht_jets + met;

  
  return ht > m_st_min;
}

Chi2Selection::Chi2Selection(Context &ctx, string label, double chi2_max):
  m_chi2_max(chi2_max),
  h_hyp(ctx.get_handle<vector<BstarToTWHypothesis>>(label)) {}

bool Chi2Selection::passes(const Event &event) {
  const BstarToTWHypothesis *hyp = get_best_hypothesis(event.get(h_hyp), "Chi2");
  if (hyp)
    {
      double chi2 = hyp->get_discriminator("Chi2");
      return chi2 < m_chi2_max;
    }
  return false;
}

MassCutSelection::MassCutSelection(double m_min_, double m_max_):
  m_min(m_min_),
  m_max(m_max_) {}

bool MassCutSelection::passes(const Event &event) {
  double m_inv = -1;
  if (event.muons->size() >= 2)
    {
      vector<Muon> muons = *event.muons;
      m_inv = (muons.at(0).v4() + muons.at(1).v4()).M();
    }
  if (event.electrons->size() >= 2)
    {
      vector<Electron> electrons = *event.electrons;
      m_inv = (electrons.at(0).v4() + electrons.at(1).v4()).M();
    }
  return (m_min < m_inv && m_inv < m_max);
}

NGenJetSelection::NGenJetSelection(unsigned int n_min_, unsigned int n_max_):
  n_min(n_min_),
  n_max(n_max_) {}

bool NGenJetSelection::passes(const Event &event) {
  return (n_min < event.genjets->size() && event.genjets->size() < n_max);
}

