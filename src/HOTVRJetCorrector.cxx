#include "UHH2/BstarToTW/include/HOTVRJetCorrector.h"

using namespace uhh2;
using namespace std;

HOTVRJetCorrector::HOTVRJetCorrector(Context &ctx) {
  mIsMC = (ctx.get("dataset_type") == "MC");
}

bool HOTVRJetCorrector::process(Event &event) {
  double rho = event.rho;
  for (auto & topjet : *event.topjets)
    {
      // Apply Jet corrections on all subjets and combine subjets back to new fatjet
      LorentzVector temp_jet;
      vector<Jet> subjets = topjet.subjets();
      for (auto & subjet : subjets)
	{
	  // Pile-up corrections
	  double jetArea = subjet.jetArea();
	  double pileup_factor = 1 - (jetArea * rho / subjet.pt());
	  subjet.set_JEC_L1factor_raw(pileup_factor);

	  double corr_factor = pileup_factor;
	  subjet.set_JEC_factor_raw(1. / corr_factor);
	  temp_jet += subjet.v4() * corr_factor;
	}
      topjet.set_v4(temp_jet);
      // double jetArea = topjet.jetArea();
      // double pileup_factor = 1 - (jetArea * rho / topjet.pt());

      // reconstruction corrections

      // LorentzVector v4_corr = topjet.v4() * pileup_factor;
      // topjet.set_v4(v4_corr);
      // topjet.set_JEC_L1factor_raw(pileup_factor);
    }
  return true;
}


HOTVRJetCorrector::~HOTVRJetCorrector() {}
