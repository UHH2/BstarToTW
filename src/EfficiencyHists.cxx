#include "UHH2/BstarToTW/include/EfficiencyHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

EfficiencyHists::EfficiencyHists(Context & ctx, const string & dirname): 
  Hists(ctx, dirname),
  h_AK8Jets(ctx.get_handle<vector<TopJet>>("slimmedJetsAK8_SoftDrop"))
{
  // book all histograms here



  // efficiency
  all           = book<TH1F>("all", "N", 2, 0, 2);
  hotvr_counter = book<TH1F>("hotvr_counter", "N", 2, 0, 2); // N_events that survivied HOTVR Cuts
  cms_counter   = book<TH1F>("cms_counter", "N", 2, 0, 2); // N_events that survived CMS Top Tagger Cuts

}

void EfficiencyHists::fill(const Event & event){  

  double weight = event.weight;

  vector<TopJet> hotvrJets = *event.topjets;
  vector<TopJet> ak8Jets = event.get(h_AK8Jets);

  all->Fill(0., weight);

  int nCMS = 0;
  for (TopJet topjet : ak8Jets)
    {
      if (topjet.pt() > 400 && 
	  abs(topjet.v4().Rapidity()) < 2.4 &&
	  topjet.tau3() / topjet.tau2() < 0.69 &&
	  topjet.softdropmass() > 110 && topjet.softdropmass() < 210)
	++nCMS;
    }
  if (nCMS == 1) cms_counter->Fill(0., weight);

  int nHOTVR = 0;
  for (TopJet topjet : hotvrJets)
    {
      vector<Jet> subjets = topjet.subjets();
      if (topjet.pt() > 200 &&
	  abs(topjet.v4().eta()) < 2.4 &&
	  subjets.size() >= 3 &&
	  topjet.v4().M() > 140 && topjet.v4().M() < 220)
	{
	  double fpt = subjets.at(0).pt()/topjet.pt();
	  double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
	  double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
	  double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
	  double mpair = min(min(m12, m13), m23);
	  if (fpt < 0.8 && mpair > 50) ++nHOTVR;
	}
    }
  if (nHOTVR == 1) hotvr_counter->Fill(0., weight);
}

EfficiencyHists::~EfficiencyHists(){}
