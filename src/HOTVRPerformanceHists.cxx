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

  // HOTVR hists
  N_HotvrTopjets          = book<TH1F>("N_HOTVR",          "N_{topjets}", 10,  0, 10);
  Pt_HotvrTopjets         = book<TH1F>("Pt_HOTVR",         "p_{t}^{topjet}", 160, 0, 1600);
  Eta_HotvrTopjets        = book<TH1F>("Eta_HOTVR",        "#eta^{topjet}", 60, -6, 6);
  M_HotvrTopjets          = book<TH1F>("M_HOTVR",          "M^{topjet}", 60,  0, 300);
  NSub_HotvrTopjets       = book<TH1F>("NSub_HOTVR",       "N_{subjets}", 10,  0, 10);
  Fpt_HotvrTopjets        = book<TH1F>("Fpt_HOTVR",        "f_{pt, 1}", 20,  0, 1);
  Mpair_HotvrTopjets      = book<TH1F>("Mpair_HOTVR",      "M_pair", 40,  0, 200);
  Tau32_HotvrTopjets      = book<TH1F>("tau32_HOTVR",      "#tau_{32}", 20,  0, 1);
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
  N_HotvrTopjets->Fill(hotvrTopTaggedInd.size(), weight);
  for (int ind : hotvrTopTaggedInd)
    {
      TopJet topjet = hotvrJets.at(ind);
      vector<Jet> subjets = topjet.subjets();

      double pt_topjet = topjet.v4().pt();
      double fpt = 1.;
      // fpt can only be calculated if there are subjets
      if (subjets.size() >= 1)
	{
	  fpt = subjets.at(0).v4().pt() / pt_topjet;
	}
      double mpair = 0;
      // mpair can only be calculated if there are at least 3 subjets
      if (subjets.size() >= 3)
	{
	  double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
	  double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
	  double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
	  mpair = min(min(m12, m13), m23);
	}

      // TH1Fs
      Pt_HotvrTopjets->Fill(pt_topjet, weight);
      Eta_HotvrTopjets->Fill(topjet.v4().eta(), weight);
      M_HotvrTopjets->Fill(topjet.v4().M(), weight);
      NSub_HotvrTopjets->Fill(subjets.size(), weight);
      Fpt_HotvrTopjets->Fill(fpt, weight);
      Mpair_HotvrTopjets->Fill(mpair, weight);
      Tau32_HotvrTopjets->Fill(topjet.tau3_groomed()/topjet.tau2_groomed(), weight);
      DeltaR_Top_HotvrTopjets->Fill(deltaR(topjet.v4(), genTop));

      // TH2Fs
      jet_area_vs_jet_pt->Fill(pt_topjet, topjet.jetArea(), weight);
      pt_reco_over_pt_top_vs_pt->Fill(pt_topjet, pt_topjet/pt_top, weight);
      pt_reco_over_pt_top_vs_npv->Fill(npv, pt_topjet/pt_top, weight);
    }

}

HOTVRPerformanceHists::~HOTVRPerformanceHists(){}
