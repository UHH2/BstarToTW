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
HOTVRPerformanceHists::HOTVRPerformanceHists(Context & ctx, const string & dirname): 
  Hists(ctx, dirname),
  h_BstarToTWGen(ctx.get_handle<BstarToTWGen>("BstarToTWgen")),
  h_TopTagIndexer(ctx.get_handle<TopTagIndexer>("TopTagIndexer")),
  h_AK8Jets(ctx.get_handle<vector<TopJet>>("slimmedJetsAK8_SoftDrop"))
{
  // book all histograms here

  DeltaR_Top_HotvrTopjets = book<TH1F>("DeltaR_Top_HOTVR", "DeltaR_{topjet, top}", 20,  0, 4);

  jet_area_vs_jet_pt         = book<TH2F>("jet_area_vs_jet_pt",         " ; p_T^{reco,top}; A_{reco_top}",             32, 0, 1600, 10, 0, 10);
  pt_reco_over_pt_top_vs_pt  = book<TH2F>("pt_reco_over_pt_top_vs_pt",  " ; p_{T}^{reco}; p_{T}^{reco} / p_{T}^{top}", 32, 0, 1600, 20, 0, 2.);
  pt_reco_over_pt_top_vs_npv = book<TH2F>("pt_reco_over_pt_top_vs_npv", " ; N_{pv}; p_{T}^{reco} / p_{T}^{top}",       50, 0, 50,   20, 0, 2.);

}

void HOTVRPerformanceHists::fill(const Event & event){  
  // Do not fill histograms if BstarToTWgen information has not been filled
  if(!event.is_valid(h_BstarToTWGen) || !event.is_valid(h_TopTagIndexer))
    {
      return;
    }
  double weight = event.weight; // event weight

  LorentzVector genTop     = event.get(h_BstarToTWGen).tbstar(); // generator top
  vector<TopJet> hotvrJets = *event.topjets;                     // hotvr topjets
  vector<TopJet> ak8Jets   = event.get(h_AK8Jets);               // softdropped ak8 topjets

  vector<int> hotvrTopTaggedInd = event.get(h_TopTagIndexer).GetIndex(); // indices of tagged hotvr topjets

  double pt_top   = genTop.pt();       // pt of generator top
  int npv         = event.pvs->size(); // number of primary vertices

  // fill historams for HOTVR topjets
  for (int ind : hotvrTopTaggedInd)
    {
      TopJet topjet = hotvrJets.at(ind);
      vector<Jet> subjets = topjet.subjets();

      double pt_topjet = topjet.v4().pt();

      DeltaR_Top_HotvrTopjets->Fill(deltaR(topjet.v4(), genTop));

      jet_area_vs_jet_pt->Fill(pt_topjet, topjet.jetArea(), weight);
      pt_reco_over_pt_top_vs_pt->Fill(pt_topjet, pt_topjet/pt_top, weight);
      pt_reco_over_pt_top_vs_npv->Fill(npv, pt_topjet/pt_top, weight);
    }

}

HOTVRPerformanceHists::~HOTVRPerformanceHists(){}
