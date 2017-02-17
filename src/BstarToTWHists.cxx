#include "UHH2/BstarToTW/include/BstarToTWHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

BstarToTWHists::BstarToTWHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here

  // Bstar
  // yet to come



}


void BstarToTWHists::fill(const Event & event){
  double weight = event.weight;

  // Topjets
  vector<TopJet> topjets = *event.topjets;
  N_top->Fill(topjets.size(), weight);
  for (TopJet topjet : topjets)
    {
      M_top->Fill(topjet.v4().M(), weight);
      pt_top->Fill(topjet.v4().pt(), weight);
      eta_top->Fill(topjet.v4().eta(), weight);
      y_top->Fill(topjet.v4().Rapidity(), weight);
      vector<Jet> subjets = topjet.subjets();
      N_subjets->Fill(subjets.size(), weight);
      double fpt = 0;
      double mmin = 0;
      if (subjets.size() > 1) 
	{
	  fpt = subjets.at(0).v4().pt()/topjet.v4().pt();
	  if (subjets.size() > 2)
	    {
	      double m01 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
	      double m02 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
	      double m12 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
	      mmin = min(min(m01, m02), m12);
	    }
	}
      fpt_1->Fill(fpt, weight);
      mpairwise->Fill(mmin, weight);
    }

}

BstarToTWHists::~BstarToTWHists(){}
