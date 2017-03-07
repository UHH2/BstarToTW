#include "UHH2/BstarToTW/include/BstarReconstruction.h"

using namespace uhh2;
using namespace std;

BstarReconstruction::BstarReconstruction(Context &ctx, const NeutrinoReconstructionMethod &NeutrinoReco, const string & name): m_NeutrinoReco(NeutrinoReco) {
  h_BstarHyps = ctx.declare_event_output<vector<BstarReconstructionHypothesis>>(name);
}

bool BstarReconstruction::process(Event &event) {
  vector<BstarReconstructionHypothesis> recohyps;
  LorentzVector topjet = event.topjets->at(0).v4();
  LorentzVector muon = event.muons->at(0).v4();
  vector<LorentzVector> neutrinos = m_NeutrinoReco(muon, event.met->v4());
  
  for( LorentzVector neutrino : neutrinos )
    {
      BstarReconstructionHypothesis hyp;
      hyp.set_Top(topjet);
      hyp.set_Muon(muon);
      hyp.set_Neutrino(neutrino);
      hyp.set_Bstar(neutrino + muon + topjet);
      recohyps.push_back(hyp);
    }

}
