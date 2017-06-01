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

  All->Fill(0., weight);

  if (topjets.size() > 0){
    const LorentzVector topjet = topjets.at(0).v4();

    double radius = 0.3;

    double dR = deltaR(top, topjet);
    double pt_jet = topjet.pt();
    double pt_top = top.pt();
    double dPt = abs(pt_jet - pt_top) / pt_top;
      

    if(dR <= radius) Matched->Fill(0., weight);
    // if(deltaPhi(top, topjet) <= M_PI/2) Matched->Fill(0., weight);
  }
}

EfficiencyHists::~EfficiencyHists(){}

GenEfficiencyHists::GenEfficiencyHists(Context & ctx, const string & dirname, const boost::optional<Event::Handle<std::vector<GenTopJet> > > &gentopjetcollection):
  Hists(ctx, dirname), h_gentopjetcollection(gentopjetcollection)
{
  // book all histograms here
  All = book<TH1F>("All", "All", 1, 0, 1); 
  Matched = book<TH1F>("Matched", "Matched", 1, 0, 1); 

  h_bstargen = ctx.get_handle<BstarToTWGen>("BstarToTWgen");
}

void GenEfficiencyHists::fill(const Event & event){  

  double weight = event.weight;
  BstarToTWGen bstargen = event.get(h_bstargen);


  const LorentzVector top = bstargen.tbstar();
  const vector<GenTopJet> & topjets = h_gentopjetcollection ? event.get(*h_gentopjetcollection) : *event.gentopjets;

  All->Fill(0., weight);

  if (topjets.size() > 0){
    const LorentzVector topjet = topjets.at(0).v4();

    double radius = 0.3;

    double dR = deltaR(top, topjet);
    double pt_jet = topjet.pt();
    double pt_top = top.pt();
    double dPt = abs(pt_jet - pt_top) / pt_top;
      

    if(dR <= radius  && dPt < 0.1 ) Matched->Fill(0., weight);
    // if(deltaPhi(top, topjet) <= M_PI/2) Matched->Fill(0., weight);
  }
}

GenEfficiencyHists::~GenEfficiencyHists(){}
