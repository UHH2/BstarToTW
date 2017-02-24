#include "UHH2/BstarToTW/include/HOTVRHists.h"
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
HOTVRHists::HOTVRHists(Context & ctx, const string & dirname): 
  Hists(ctx, dirname),
  h_TopTagIndexer(ctx.get_handle<TopTagIndexer>("TopTagIndexer"))
{
  // book all histograms here

  // HOTVR hists
  N_HotvrTopjets          = book<TH1F>("N_HOTVR",          "N_{topjets}", 10,  0, 10);
  Pt_HotvrTopjets         = book<TH1F>("Pt_HOTVR",         "p_{t}^{topjet}", 80, 0, 1600);
  Eta_HotvrTopjets        = book<TH1F>("Eta_HOTVR",        "#eta^{topjet}", 30, -6, 6);
  M_HotvrTopjets          = book<TH1F>("M_HOTVR",          "M^{topjet}", 60,  0, 600);
  NSub_HotvrTopjets       = book<TH1F>("NSub_HOTVR",       "N_{subjets}", 10,  0, 10);
  Fpt_HotvrTopjets        = book<TH1F>("Fpt_HOTVR",        "f_{pt, 1}", 20,  0, 1);
  Mpair_HotvrTopjets      = book<TH1F>("Mpair_HOTVR",      "M_pair", 40,  0, 200);
  Tau32_HotvrTopjets      = book<TH1F>("tau32_HOTVR",      "#tau_{32}", 20,  0, 1);

  DeltaR_L_HotvrTopjets   = book<TH1F>("DeltaR_L_HOTVR",   "#Delta R_{l,t}", 20, 0, 4);
  DeltaPhi_L_HotvrTopjets = book<TH1F>("DeltaPhi_L_HOTVR", "#Delta #phi_{l,t}", 20, 0, 4);

  // b-jets
  N_bjets_loose           = book<TH1F>("N_bjets_loose",    "N_{bjets}", 10,  0, 10);
  N_bjets_medium          = book<TH1F>("N_bjets_medium",   "N_{bjets}", 10,  0, 10);
  N_bjets_tight           = book<TH1F>("N_bjets_tight",    "N_{bjets}", 10,  0, 10);

}

void HOTVRHists::fill(const Event & event){  
  // Do not fill histograms if BstarToTWgen information has not been filled
  if(!event.is_valid(h_TopTagIndexer))
    {
      return;
    }
  double weight = event.weight; // event weight

  vector<TopJet> hotvrJets = *event.topjets;                     // hotvr topjets
  vector<int> hotvrTopTaggedInd = event.get(h_TopTagIndexer).GetIndex(); // indices of tagged hotvr topjets
  vector<Muon> muons = *event.muons;
  vector<Jet> jets = *event.jets;
  const CSVBTag btag_loose(CSVBTag::WP_LOOSE);
  const CSVBTag btag_medium(CSVBTag::WP_MEDIUM);
  const CSVBTag btag_tight(CSVBTag::WP_TIGHT);

  // fill HOTVR hists
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
      if (muons.size() > 0)
	{
	  DeltaR_L_HotvrTopjets->Fill(deltaR(topjet.v4(), muons.at(0).v4()));
	  DeltaPhi_L_HotvrTopjets->Fill(deltaPhi(topjet.v4(), muons.at(0).v4()));
	}
    }

  // fill bjet hists
  int n_loose  = 0;
  int n_medium = 0;
  int n_tight  = 0;

  for (Jet jet : jets)
    {
      if (btag_loose(jet, event)) ++n_loose;
      if (btag_medium(jet, event)) ++n_medium;
      if (btag_tight(jet, event)) ++n_tight;
    }

  N_bjets_loose->Fill(n_loose, weight);
  N_bjets_medium->Fill(n_medium, weight);
  N_bjets_tight->Fill(n_tight, weight);

}

HOTVRHists::~HOTVRHists(){}
