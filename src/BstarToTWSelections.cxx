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

Chi2Selection::Chi2Selection(Context &ctx, string label, double chi2_max, const string disc_name):
  m_chi2_max(chi2_max),
  h_hyp(ctx.get_handle<vector<BstarToTWHypothesis>>(label)),
  m_disc_name(disc_name){}

bool Chi2Selection::passes(const Event &event) {
  const BstarToTWHypothesis *hyp = get_best_hypothesis(event.get(h_hyp), m_disc_name);
  if (hyp)
    {
      double chi2 = hyp->get_discriminator(m_disc_name);
      return chi2 < m_chi2_max;
    }
  return false;
}

RecoMassSelection::RecoMassSelection(Context &ctx, double m_min_, string label):
  m_min(m_min_),
  h_hyp(ctx.get_handle<vector<BstarToTWHypothesis>>(label)) {}

bool RecoMassSelection::passes(const Event &event) {
  std::vector<BstarToTWHypothesis> hyps = event.get(h_hyp);
  const BstarToTWHypothesis* hyp = get_best_hypothesis( hyps, "Chi2" );
  if (!hyp)
    {
      cout << "WARNING: RecoMassSelection: Chi2 No hypothesis was valid!" << endl;
      return false;
    }
  double mbstar_reco = 0;
  if((hyp->get_topjet() + hyp->get_w()).isTimelike())
    {    
      mbstar_reco = (hyp->get_topjet() + hyp->get_w()).M();
    }
  else
    {
      mbstar_reco = sqrt(-(hyp->get_topjet()+hyp->get_w()).mass2());
    }
  return mbstar_reco > m_min;
}

JetDeltaPhiSelection::JetDeltaPhiSelection(Context &ctx, double delta_phi_min, const boost::optional<Event::Handle<std::vector<Jet> > > jet_collection) {
  m_delta_phi_min = delta_phi_min;  
  h_jets = jet_collection;
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
}

bool JetDeltaPhiSelection::passes(const Event &event) {
  const vector<Jet> &jets = h_jets ? event.get(*h_jets) : *event.jets;
  const FlavorParticle &muo = event.get(h_primlep);
  if (jets.size() > 0)
    return (deltaPhi(jets.at(0), muo) > m_delta_phi_min);
  else
    return false;
}
LeadingJetDeltaRSelection::LeadingJetDeltaRSelection(Context &ctx, double delta_R_min, const boost::optional<Event::Handle<std::vector<Jet> > > jet_collection) {
  m_delta_R_min = delta_R_min;  
  h_jets = jet_collection;
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
}

bool LeadingJetDeltaRSelection::passes(const Event &event) {
  const vector<Jet> &jets = h_jets ? event.get(*h_jets) : *event.jets;
  const FlavorParticle &lep = event.get(h_primlep);
  return (jets.size() > 0 && deltaR(lep, jets.at(0)) > m_delta_R_min);
}

JetDeltaRSelection::JetDeltaRSelection(Context &ctx, double delta_R_min, const boost::optional<Event::Handle<std::vector<Jet> > > jet_collection) {
  m_delta_R_min = delta_R_min;  
  h_jets = jet_collection;
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
}

bool JetDeltaRSelection::passes(const Event &event) {
  const vector<Jet> &jets = h_jets ? event.get(*h_jets) : *event.jets;
  const FlavorParticle &lep = event.get(h_primlep);
  for (const Jet &jet : jets)
    {
      if (deltaR(lep, jet) < m_delta_R_min)
	return false;
    }
  return true;
}

TopJetDeltaRSelection::TopJetDeltaRSelection(Context &ctx, double delta_R_max, const boost::optional<Event::Handle<std::vector<TopJet> > > topjet_collection) {
  m_delta_R_max = delta_R_max;
  do_R_calc = delta_R_max < 0;
  h_topjets = topjet_collection;
  h_bjets = ctx.get_handle<vector<Jet> >("btag_loose");
}

bool TopJetDeltaRSelection::passes(const Event &event) {
  const vector<TopJet> &topjets = h_topjets ? event.get(*h_topjets) : *event.topjets;
  const vector<Jet> &bjets = event.get(h_bjets);
  if (bjets.size() > 0)
    {
      for (const TopJet topjet : topjets)
	{
	  if (do_R_calc)
	    m_delta_R_max = sqrt(topjet.jetArea() / M_PI);
	  if (deltaR(bjets.at(0), topjet) > m_delta_R_max)
	    return false;
	}
      return true;
    }
  else return false;
}

METDeltaPhiSelection::METDeltaPhiSelection(Context &ctx, double delta_phi_max) {
  m_delta_phi_max = delta_phi_max;  
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
}

bool METDeltaPhiSelection::passes(const Event &event) {
  const auto met = *event.met;
  const FlavorParticle &muo = event.get(h_primlep);
  return (deltaPhi(met, muo) < m_delta_phi_max);
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

