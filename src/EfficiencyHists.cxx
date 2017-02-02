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
  hotvr_counter = book<TH1F>("hotvr_counter", "N", 2, 0, 2); // N_events that survivied HOTVR Cuts
  cms_counter   = book<TH1F>("cms_counter", "N", 2, 0, 2); // N_events that survived CMS Top Tagger Cuts

}

void EfficiencyHists::fill(const Event & event){  

  double weight = event.weight;

  vector<TopJet> hotvrJets = *event.topjets;
  vector<TopJet> ak8Jets = event.get(h_AK8Jets);

  int nCMS = 0;
  for (TopJet topjet : ak8Jets)
    {
      if (topjet.pt() > 400 && 
	  topjet.v4().Rapidity() < 2.4 &&
	  topjet.tau3()/topjet.tau2() < 0.69 &&
	  topjet.softdropmass() > 110 && topjet.softdropmass() < 210)
	++nCMS;
    }
  if (nCMS == 1) cms_counter->Fill(0., weight);

  int nHOTVR = 0;
  for (TopJet topjet : hotvrJets)
    {
      if (topjet.v4().M() > 140 && topjet.v4().M() < 220)
	++nHOTVR;
    }
  if (nHOTVR == 1) hotvr_counter->Fill(0., weight);
}

EfficiencyHists::~EfficiencyHists(){}
