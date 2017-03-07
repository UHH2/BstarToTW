#include "UHH2/BstarToTW/include/EfficiencyHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

EfficiencyHists::EfficiencyHists(Context & ctx, const string & dirname, const boost::optional<Event::Handle<std::vector<TopJet> > > &topjetcollection):
  Hists(ctx, dirname), h_topjetcollection(topjetcollection)
{
  // book all histograms here
  All = book<TH1F>("All", "All", 1, 0, 1); 
  Matched = book<TH1F>("Matched", "Matched", 1, 0, 1); 

  h_bstargen = ctx.get_handle<BstarToTWGen>("BstarToTWgen");
}

void EfficiencyHists::fill(const Event & event){  

  double weight = event.weight;

  BstarToTWGen bstargen = event.get(h_bstargen);

  const LorentzVector top = bstargen.tbstar();
  const vector<TopJet> & topjets = h_topjetcollection ? event.get(*h_topjetcollection) : *event.topjets;
  const LorentzVector topjet = topjets.at(0).v4();
  All->Fill(0., weight);
  // ToDo: find better way for matching
  if(deltaPhi(top, topjet) < M_PI/2) Matched->Fill(0., weight);
}

EfficiencyHists::~EfficiencyHists(){}
