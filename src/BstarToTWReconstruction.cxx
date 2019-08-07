#include "UHH2/BstarToTW/include/BstarToTWReconstruction.h"

using namespace uhh2;
using namespace std;

BstarToTWReconstruction::BstarToTWReconstruction(Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, const string & label, TopJetId tjetid):
  m_neutrinofunction(neutrinofunction),
  h_recohyps(ctx.get_handle<vector<BstarToTWHypothesis>>(label)),
  h_primlep(ctx.get_handle<FlavorParticle>("PrimaryLepton")),
  m_topjetid(tjetid) {}


bool BstarToTWReconstruction::process(uhh2::Event & event) {
  assert(event.topjets);
  assert(event.met);

  vector<BstarToTWHypothesis> recoHyps;

  const Particle& lepton = event.get(h_primlep);
  vector<LorentzVector> neutrinos = m_neutrinofunction(lepton.v4(), event.met->v4());

  for(const auto& topjet : *event.topjets)
    {

      if(!m_topjetid(topjet, event)) continue;

      for(const auto& neutrino : neutrinos)
	{
	  BstarToTWHypothesis hyp;
	  hyp.set_lepton(lepton.v4());
	  hyp.set_neutrino(neutrino);
          hyp.set_topjet(topjet.v4());

          recoHyps.push_back(hyp);
	}
    }

  event.set(h_recohyps, recoHyps);
  return true;
}

LeptonicTopReconstruction::LeptonicTopReconstruction(Context &ctx, const string &hyps_name, const string &bjet_name, const string &label):
  h_recohyps(ctx.get_handle<vector<BstarToTWHypothesis>>(hyps_name)),
  h_bjets(ctx.get_handle<vector<Jet> >(bjet_name)),
  h_toplephyps(ctx.get_handle<vector<LeptonicTopHypothesis> >(label)){}

bool LeptonicTopReconstruction::process(Event &event) {
  
  vector<BstarToTWHypothesis> reco_hyps = event.get(h_recohyps);
  vector<Jet> bjets = event.get(h_bjets);
  vector<LeptonicTopHypothesis> toplephyps;
  bool fast_reco = (bjets.size() >= 2); // fast reconstruction can be done when there are two b jets

  for (BstarToTWHypothesis &bstar_hyp : reco_hyps)
    {
      LeptonicTopHypothesis hyp;
      LorentzVector tophad = bstar_hyp.get_topjet();
      LorentzVector w = bstar_hyp.get_w();
      hyp.set_lepton(bstar_hyp.get_lepton());
      hyp.set_neutrino(bstar_hyp.get_neutrino());
      hyp.set_tophad(tophad);
      if (fast_reco)
	{
	  LorentzVector bjet = (deltaR(bjets.at(0).v4(), tophad) > deltaR(bjets.at(1).v4(), tophad))  ? bjets.at(0).v4() : bjets.at(1).v4();
	  hyp.set_toplep(bjet + w);
	  toplephyps.push_back(hyp);
	}
      else
	{
	  for (Jet &jet : *event.jets)
	    {
	      if (deltaR(jet, tophad) < 1.9) continue;
	      hyp.set_toplep(jet.v4() + w);
	      toplephyps.push_back(hyp);
	    }
	}
    }
  event.set(h_toplephyps, toplephyps);
  return true;
}
