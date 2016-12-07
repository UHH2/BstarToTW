#include "UHH2/BstarToTW/include/BstarToTWHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

BstarToTWHists::BstarToTWHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  N_jets   = book<TH1F>("N_jets", "N_{jets}", 20, 0, 20);  
  eta_jet1 = book<TH1F>("eta_jet1", "#eta^{jet 1}", 40, -2.5, 2.5);
  eta_jet2 = book<TH1F>("eta_jet2", "#eta^{jet 2}", 40, -2.5, 2.5);
  eta_jet3 = book<TH1F>("eta_jet3", "#eta^{jet 3}", 40, -2.5, 2.5);
  pt_jet1  = book<TH1F>("pt_jet1", "p_{T}^{jet 1}", 40, 0, 200);
  pt_jet2  = book<TH1F>("pt_jet2", "p_{T}^{jet 2}", 40, 0, 200);
  pt_jet3  = book<TH1F>("pt_jet3", "p_{T}^{jet 3}", 40, 0, 200);

  // leptons
  N_mu      = book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  pt_mu     = book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV/c]", 40, 0, 200);
  eta_mu    = book<TH1F>("eta_mu", "#eta^{#mu}", 40, -2.1, 2.1);
  reliso_mu = book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);

  // primary vertices
  N_pv = book<TH1F>("N_pv", "N^{PV}", 50, 0, 50);
}


void BstarToTWHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
  
  std::vector<Jet>* jets = event.jets;
  int Njets = jets->size();
  N_jets->Fill(Njets, weight);
  
  if(Njets>=1){
    eta_jet1->Fill(jets->at(0).eta(), weight);
    pt_jet1->Fill(jets->at(0).pt(), weight);
  }
  if(Njets>=2){
    eta_jet2->Fill(jets->at(1).eta(), weight);
    pt_jet2->Fill(jets->at(1).pt(), weight);
  }
  if(Njets>=3){
    eta_jet3->Fill(jets->at(2).eta(), weight);
    pt_jet3->Fill(jets->at(2).pt(), weight);
  }


  int Nmuons = event.muons->size();
  N_mu->Fill(Nmuons, weight);
  for (const Muon & thismu : *event.muons){
      pt_mu->Fill(thismu.pt(), weight);
      eta_mu->Fill(thismu.eta(), weight);
      reliso_mu->Fill(thismu.relIso(), weight);
  }
  
  int Npvs = event.pvs->size();
  N_pv->Fill(Npvs, weight);
}

BstarToTWHists::~BstarToTWHists(){}
