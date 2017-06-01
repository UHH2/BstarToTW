#include "UHH2/BstarToTW/include/BstarToTWReconstruction.h"

using namespace uhh2;
using namespace std;

BstarToTWReconstruction::BstarToTWReconstruction(Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, const string & label, TopJetId tjetid):
  m_neutrinofunction(neutrinofunction),
  h_recohyps(ctx.declare_event_output<vector<BstarToTWHypothesis>>(label)),
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
