#include "UHH2/BstarToTW/include/HOTVRPerformanceHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <vector>

using namespace std;
using namespace uhh2;

/*
 * WARNING: Fill this Hists only after cuts are applied to ensure
 * there is >=1 HOTVR TopJet.
 * 
 * This Hists Class implements Histograms with informations about the
 * physics performance of the HOTVR algorithm.
 * 
 */
HOTVRPerformanceHists::HOTVRPerformanceHists(Context & ctx, const string & dirname, const boost::optional<Event::Handle<std::vector<TopJet> > > &topjetcollection): 
Hists(ctx, dirname), 
  h_BstarToTWGen(ctx.get_handle<BstarToTWGen>("BstarToTWgen")),
  h_topjetcollection(topjetcollection)
{
  // book all histograms here

  DeltaR_Top_HotvrTopjet = book<TH1F>("DeltaR_Top_HOTVR", "DeltaR_{gen, jet}", 20,  0, 4);

  EffPt_Top_HotvrTopjet  = book<TH1F>("EffPt_Top_HOTVR", "(p_{T}^{jet} - p_{T}^{top})/p_{T}^{top}", 40,  -0.2, 0.2);

  EffPt_Top_HotvrTopjet_vs_pt_top  = book<TH2F>("EffPt_Top_HOTVR_vs_pt_top", "(p_{T}^{jet} - p_{T}^{top})/p_{T}^{top}", 32, 0, 1600, 40,  -0.2, 0.2);
  EffPt_Top_HotvrTopjet_vs_npv  = book<TH2F>("EffPt_Top_HOTVR_vs_npv", "(p_{T}^{jet} - p_{T}^{top})/p_{T}^{top}", 32, 0, 1600, 40,  -0.2, 0.2);
  // pt_reco_over_pt_top_vs_pt_reco  = book<TH2F>("pt_reco_over_pt_top_vs_pt_reco",  " ; p_{T}^{reco}; p_{T}^{reco} / p_{T}^{top}", 32, 0, 1600, 20, 0, 2.);
  // pt_reco_over_pt_top_vs_pt_top  = book<TH2F>("pt_reco_over_pt_top_vs_pt_top",  " ; p_{T}^{reco}; p_{T}^{reco} / p_{T}^{top}", 32, 0, 1600, 20, 0, 2.);
  // pt_reco_over_pt_top_vs_npv = book<TH2F>("pt_reco_over_pt_top_vs_npv", " ; N_{pv}; p_{T}^{reco} / p_{T}^{top}",       50, 0, 50,   20, 0, 2.);

}

void HOTVRPerformanceHists::fill(const Event & event){  
  // Do not fill histograms if BstarToTWgen information has not been filled
  if(!event.is_valid(h_BstarToTWGen))
    {
      return;
    }
  double weight = event.weight; // event weight

  LorentzVector genTop     = event.get(h_BstarToTWGen).tbstar(); // generator top
  const vector<TopJet> & topjets = h_topjetcollection ? event.get(*h_topjetcollection) : *event.topjets;

  double pt_top   = genTop.pt();       // pt of generator top
  int npv         = event.pvs->size(); // number of primary vertices

  // fill historams for HOTVR topjets
  for (TopJet topjet : topjets)
    {
      double pt_topjet = topjet.v4().pt();
      double EffPt = (pt_topjet-pt_top)/pt_top;

      DeltaR_Top_HotvrTopjet->Fill(deltaR(topjet.v4(), genTop), weight);
      EffPt_Top_HotvrTopjet->Fill(EffPt, weight);
      EffPt_Top_HotvrTopjet_vs_pt_top->Fill(pt_top, EffPt, weight);
      EffPt_Top_HotvrTopjet_vs_npv->Fill(npv, EffPt, weight);

    }

}

HOTVRPerformanceHists::~HOTVRPerformanceHists(){}


GenHOTVRPerformanceHists::GenHOTVRPerformanceHists(Context & ctx, const string & dirname, const boost::optional<Event::Handle<std::vector<GenTopJet> > > &topjetcollection): 
Hists(ctx, dirname), 
  h_BstarToTWGen(ctx.get_handle<BstarToTWGen>("BstarToTWgen")),
  h_topjetcollection(topjetcollection)
{
  // book all histograms here

  DeltaR_Top_HotvrTopjet = book<TH1F>("DeltaR_Top_HOTVR", "DeltaR_{gen, jet}", 20,  0, 4);

  EffPt_Top_HotvrTopjet  = book<TH1F>("EffPt_Top_HOTVR", "(p_{T}^{genjet} - p_{T}^{top})/p_{T}^{top}", 40,  -0.2, 0.2);

  EffPt_Top_HotvrTopjet_vs_pt_top  = book<TH2F>("EffPt_Top_HOTVR_vs_pt_top", "(p_{T}^{genjet} - p_{T}^{top})/p_{T}^{top}", 32, 0, 1600, 40,  -0.2, 0.2);
  EffPt_Top_HotvrTopjet_vs_npv  = book<TH2F>("EffPt_Top_HOTVR_vs_npv", "(p_{T}^{genjet} - p_{T}^{top})/p_{T}^{top}", 32, 0, 1600, 40,  -0.2, 0.2);
  // pt_reco_over_pt_top_vs_pt_reco  = book<TH2F>("pt_reco_over_pt_top_vs_pt_reco",  " ; p_{T}^{reco}; p_{T}^{reco} / p_{T}^{top}", 32, 0, 1600, 20, 0, 2.);
  // pt_reco_over_pt_top_vs_pt_top  = book<TH2F>("pt_reco_over_pt_top_vs_pt_top",  " ; p_{T}^{reco}; p_{T}^{reco} / p_{T}^{top}", 32, 0, 1600, 20, 0, 2.);
  // pt_reco_over_pt_top_vs_npv = book<TH2F>("pt_reco_over_pt_top_vs_npv", " ; N_{pv}; p_{T}^{reco} / p_{T}^{top}",       50, 0, 50,   20, 0, 2.);

}

void GenHOTVRPerformanceHists::fill(const Event & event){  
  // Do not fill histograms if BstarToTWgen information has not been filled
  if(!event.is_valid(h_BstarToTWGen))
    {
      return;
    }
  double weight = event.weight; // event weight

  LorentzVector genTop     = event.get(h_BstarToTWGen).tbstar(); // generator top
  const vector<GenTopJet> & topjets = h_topjetcollection ? event.get(*h_topjetcollection) : *event.gentopjets;

  double pt_top   = genTop.pt();       // pt of generator top
  int npv         = event.pvs->size(); // number of primary vertices

  // fill historams for HOTVR topjets
  for (GenTopJet topjet : topjets)
    {
      double pt_topjet = topjet.v4().pt();
      double EffPt = (pt_topjet-pt_top)/pt_top;

      DeltaR_Top_HotvrTopjet->Fill(deltaR(topjet.v4(), genTop), weight);
      EffPt_Top_HotvrTopjet->Fill(EffPt, weight);
      EffPt_Top_HotvrTopjet_vs_pt_top->Fill(pt_top, EffPt, weight);
      EffPt_Top_HotvrTopjet_vs_npv->Fill(npv, EffPt, weight);

    }

}

GenHOTVRPerformanceHists::~GenHOTVRPerformanceHists(){}

