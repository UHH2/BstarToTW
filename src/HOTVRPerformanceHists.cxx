#include "UHH2/BstarToTW/include/HOTVRPerformanceHists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
/*
 * Fill this class only after cuts are applied to ensure there is ==1 HOTVR TopJet.
 */
HOTVRPerformanceHists::HOTVRPerformanceHists(Context & ctx, const string & dirname): 
  Hists(ctx, dirname),
  h_BstarToTWgen(ctx.get_handle<BstarToTWGen>("BstarToTWgen")) 
{
  // book all histograms here

  // performance hists
  h_pt_reco_over_pt_top = book<TH1F>("pt_reco_over_pt_top", "p_{T}^{reco} / p_{T}^{top, gen}", 20, 0, 2);
  h_pt_reco_over_pt_gen = book<TH1F>("pt_reco_over_pt_gen", "p_{T}^{reco} / p_{T}^{gen}", 20, 0, 2);

  h_pt_reco_over_pt_top_vs_pt  = book<TH2F>("pt_reco_over_pt_top_vs_pt", "x= p_{T}^{reco} y=p_{T}^{reco} / p_{T}^{top}", 32, 0, 1600., 20, 0, 2.);
  h_pt_reco_over_pt_gen_vs_pt  = book<TH2F>("pt_reco_over_pt_gen_vs_pt", "x= p_{T}^{reco} y=p_{T}^{reco} / p_{T}^{gen}", 32, 0, 1600., 20, 0, 2.);
  h_pt_reco_over_pt_top_vs_npv = book<TH2F>("pt_reco_over_pt_top_vs_npv", "x= N_{primary vertices} y=p_{T}^{reco} / p_{T}^{top}", 50, 0, 50., 20, 0, 2.);
  h_pt_reco_over_pt_gen_vs_npv = book<TH2F>("pt_reco_over_pt_gen_vs_npv", "x= N_{primary vertices} y=p_{T}^{reco} / p_{T}^{gen}", 50, 0, 50., 20, 0, 2.);

  h_deltaR_reco_top     = book<TH1F>("deltaR_reco_top", "#Delta R_{reco,top}", 20, 0, 4);
  h_deltaR_reco_gen     = book<TH1F>("deltaR_reco_gen", "#Delta R_{reco,gen}", 20, 0, 4);
  h_deltaR_gen_top      = book<TH1F>("deltaR_gen_top", "#Delta R_{gen,top}", 20, 0, 4);

  h_matched_pt_reco     = book<TH1F>("matched_pt", "p_{T}^{reco,top} [GeV/c]", 80, 0, 1600);
  h_matched_M_reco      = book<TH1F>("matched_M", "M^{reco,top} [GeV/c^2]", 40, 0, 400);
  h_unmatched_pt_reco   = book<TH1F>("unmatched_pt", "p_{T}^{reco,top} [GeV/c]", 80, 0, 1600);
  h_unmatched_M_reco    = book<TH1F>("unmatched_M", "M^{reco,top} [GeV/c^2]", 40, 0, 400);

  h_jet_area_vs_jet_pt  = book<TH2F>("jet_area_vs_jet_pt", "x= p_T^{reco,top} y= A_{reco_top}", 32, 0, 1600, 10, 0, 10);  

}

void HOTVRPerformanceHists::fill(const Event & event){  
  // Do not fill histograms if BstarToTWgen information has not been filled
  if(!event.is_valid(h_BstarToTWgen))
    {
      return;
    }
  double weight = event.weight;
  TopJet hotvrTop = event.topjets->at(0);
  GenTopJet hotvrGenTop = event.gentopjets->at(0);
  LorentzVector genTop = event.get(h_BstarToTWgen).tbstar();


  // To Do: pt and npv binning, probably with TH2F
  double pt_reco = hotvrTop.pt();
  double pt_gen = hotvrGenTop.pt();
  double pt_top = genTop.pt();
  double pt_reco_over_pt_top = pt_reco / pt_top;
  double pt_reco_over_pt_gen = pt_reco / pt_gen;
  int npv = event.pvs->size();

  h_pt_reco_over_pt_top->Fill(pt_reco_over_pt_top, weight);
  h_pt_reco_over_pt_gen->Fill(pt_reco_over_pt_gen, weight);

  h_pt_reco_over_pt_top_vs_pt->Fill(pt_reco, pt_reco_over_pt_top, weight);
  h_pt_reco_over_pt_gen_vs_pt->Fill(pt_reco, pt_reco_over_pt_gen, weight);
 
  h_pt_reco_over_pt_top_vs_npv->Fill(npv, pt_reco_over_pt_top, weight);
  h_pt_reco_over_pt_gen_vs_npv->Fill(npv, pt_reco_over_pt_gen, weight);


  double deltaR_reco_top = deltaR( hotvrTop.v4(), genTop);
  double deltaR_reco_gen = deltaR( hotvrTop.v4(), hotvrGenTop.v4());
  double deltaR_gen_top = deltaR( hotvrGenTop.v4(), genTop);
  h_deltaR_reco_top->Fill(deltaR_reco_top, weight);
  h_deltaR_reco_gen->Fill(deltaR_reco_gen, weight);
  h_deltaR_gen_top->Fill(deltaR_gen_top, weight);
  
  if (deltaR_reco_top < 1.5)
    {
      h_matched_pt_reco->Fill(pt_reco, weight);
      h_matched_M_reco->Fill(hotvrTop.v4().M(), weight);
    }
  else 
    {
      h_unmatched_pt_reco->Fill(pt_reco, weight);
      h_unmatched_M_reco->Fill(hotvrTop.v4().M(),weight);      
    }

  h_jet_area_vs_jet_pt->Fill(pt_reco, hotvrTop.jetArea(), weight);

}

HOTVRPerformanceHists::~HOTVRPerformanceHists(){}
