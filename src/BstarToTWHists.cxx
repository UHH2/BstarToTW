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

  // HOTVR Topjets
  N_top     = book<TH1F>("N_topjets", "N_{topjets}", 5, 0, 5);
  M_top     = book<TH1F>("M_top", "M_{top} [GeV/c^{2}]", 40, 0, 400);
  pt_top    = book<TH1F>("pt_top", "p_{T}^{top} [GeV/c]", 80, 0, 1600);
  // substructure variables
  N_subjets = book<TH1F>("N_subjets", "N_{subjets}", 10, 0, 10);
  mpairwise = book<TH1F>("mpairwise", "m_{min}^{pairwise} [GeV/c^{2}]", 20, 0, 200);
  fpt_1     = book<TH1F>("fpt", "f_{pt,1}", 10, 0, 1);


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
