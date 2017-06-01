#include "UHH2/BstarToTW/include/HOTVRHists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

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
HOTVRHists::HOTVRHists(Context & ctx, const string & dirname, const boost::optional<TopJetId> &topjetid): 
  Hists(ctx, dirname), m_topjetid(topjetid)
{
  // book all histograms here

  // HOTVR hists
  N_HotvrTopjets          = book<TH1F>("N_HOTVR",           "N_{topjets}", 10,  0, 10);
  Pt_HotvrTopjets         = book<TH1F>("Pt_HOTVR",          "p_{T}^{topjet} [GeV/c]", 40, 0, 1600);
  Eta_HotvrTopjets        = book<TH1F>("Eta_HOTVR",         "#eta^{topjet}", 30, -6, 6);
  M_HotvrTopjets          = book<TH1F>("M_HOTVR",           "M^{topjet} [GeV/c^{2}]", 40,  0, 400);
  A_HotvrTopjets          = book<TH1F>("A_HOTVR",           "A^{topjet} [a.u.]", 100, 0, 10);
  NSub_HotvrTopjets       = book<TH1F>("NSub_HOTVR",        "N_{subjets}", 10,  0, 10);
  Fpt_HotvrTopjets        = book<TH1F>("Fpt_HOTVR",         "f_{p_{T}, 1}", 20,  0, 1);
  Mpair_HotvrTopjets      = book<TH1F>("Mpair_HOTVR",       "M_pair [GeV/c^{2}]", 40,  0, 200);
  Tau32_HotvrTopjets      = book<TH1F>("tau32_HOTVR",       "#tau_{32}", 20,  0, 1);
  Pt_HotvrTopjets_Sub1    = book<TH1F>("Pt_HOTVR_Subjet1",  "p_{T}^{subjet1} [GeV/c]", 100, 0, 1000);
  A_HotvrTopjets_Sub1     = book<TH1F>("A_HOTVR_Subjet1",   "A^{subjetjet1} [a.u.]", 50, 0, 5);
  Pt_HotvrTopjets_Sub2    = book<TH1F>("Pt_HOTVR_Subjet2",  "p_{T}^{subjet2} [GeV/c]", 100, 0, 1000);
  A_HotvrTopjets_Sub2     = book<TH1F>("A_HOTVR_Subjet2",   "A^{subjetjet2} [a.u.]", 50, 0, 5);
  Pt_HotvrTopjets_Sub3    = book<TH1F>("Pt_HOTVR_Subjet3",  "p_{T}^{subjet3} [GeV/c]", 100, 0, 1000);
  A_HotvrTopjets_Sub3     = book<TH1F>("A_HOTVR_Subjet3",   "A^{subjetjet3} [a.u.]", 50, 0, 5);

  DeltaR_L_HotvrTopjets   = book<TH1F>("DeltaR_L_HOTVR",   "#Delta R_{l,t}", 20, 0, 4);
  DeltaPhi_L_HotvrTopjets = book<TH1F>("DeltaPhi_L_HOTVR", "#Delta #phi_{l,t}", 20, 0, 4);

  Pt_HotvrTopjet1         = book<TH1F>("Pt_HOTVR1",         "p_{T}^{topjet} [GeV/c]", 40, 0, 1600);
  Eta_HotvrTopjet1        = book<TH1F>("Eta_HOTVR1",        "#eta^{topjet}", 30, -6, 6);
  M_HotvrTopjet1          = book<TH1F>("M_HOTVR1",          "M^{topjet} [GeV/c^{2}]", 60,  0, 400);
  A_HotvrTopjet1          = book<TH1F>("A_HOTVR1",           "A^{topjet} [a.u.]", 100, 0, 10);
  NSub_HotvrTopjet1       = book<TH1F>("NSub_HOTVR1",       "N_{subjets}", 10,  0, 10);
  Fpt_HotvrTopjet1        = book<TH1F>("Fpt_HOTVR1",        "f_{pt, 1}", 20,  0, 1);
  Mpair_HotvrTopjet1      = book<TH1F>("Mpair_HOTVR1",      "M_pair [GeV/c^{2}]", 40,  0, 200);
  Tau32_HotvrTopjet1      = book<TH1F>("tau32_HOTVR1",      "#tau_{32}", 20,  0, 1);
  Pt_HotvrTopjet1_Sub1    = book<TH1F>("Pt_HOTVR1_Subjet1", "p_{T}^{subjet1} [GeV/c]", 100, 0, 1000);
  A_HotvrTopjet1_Sub1     = book<TH1F>("A_HOTVR1_Subjet1",   "A^{subjetjet1} [a.u.]", 50, 0, 5);
  Pt_HotvrTopjet1_Sub2    = book<TH1F>("Pt_HOTVR1_Subjet2", "p_{T}^{subjet2} [GeV/c]", 100, 0, 1000);
  A_HotvrTopjet1_Sub2     = book<TH1F>("A_HOTVR1_Subjet2",   "A^{subjetjet2} [a.u.]", 50, 0, 5);
  Pt_HotvrTopjet1_Sub3    = book<TH1F>("Pt_HOTVR1_Subjet3", "p_{T}^{subjet3} [GeV/c]", 100, 0, 1000);
  A_HotvrTopjet1_Sub3     = book<TH1F>("A_HOTVR1_Subjet3",   "A^{subjetjet3} [a.u.]", 50, 0, 5);

  DeltaR_L_HotvrTopjet1   = book<TH1F>("DeltaR_L_HOTVR1",   "#Delta R_{l,t}", 20, 0, 4);
  DeltaPhi_L_HotvrTopjet1 = book<TH1F>("DeltaPhi_L_HOTVR1", "#Delta #phi_{l,t}", 20, 0, 4);

  // b-jets
  N_bjets_loose           = book<TH1F>("N_bjets_loose",    "N_{bjets}", 10,  0, 10);
  N_bjets_medium          = book<TH1F>("N_bjets_medium",   "N_{bjets}", 10,  0, 10);
  N_bjets_tight           = book<TH1F>("N_bjets_tight",    "N_{bjets}", 10,  0, 10);

}

void HOTVRHists::fill(const Event & event) {  

  double weight = event.weight; // event weight

  vector<TopJet> hotvrJets = *event.topjets;
  vector<Muon> muons = *event.muons;
  sort_by_pt<Muon>(muons);
  vector<Jet> jets = *event.jets;
  
  int n_jets = 0;

  // fill HOTVR hists
  for (TopJet topjet : hotvrJets)
    {
      if (m_topjetid)
	{
	  if (!(*m_topjetid)(topjet, event)) continue;
	}
      ++n_jets;
      vector<Jet> subjets = topjet.subjets();

      double pt_topjet = topjet.v4().pt();
      double a_topjet = topjet.jetArea();
      double fpt = -1;
      // fpt can only be calculated if there are subjets
      if (subjets.size() >= 1)
	{
	  fpt = subjets.at(0).v4().pt() / pt_topjet;
	}
      double mpair = -1;
      double pt_sub1 = -1;
      double a_sub1 = -1;
      double pt_sub2 = -1;
      double a_sub2 = -1;
      double pt_sub3 = -1;
      double a_sub3 = -1;

      // mpair can only be calculated if there are at least 3 subjets
      if (subjets.size() >= 3)
	{
	  double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
	  double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
	  double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
	  mpair = min(min(m12, m13), m23);
	  pt_sub1 = subjets.at(0).v4().pt();
	  a_sub1  = subjets.at(0).jetArea();
	  pt_sub2 = subjets.at(1).v4().pt();
	  a_sub2  = subjets.at(1).jetArea();
	  pt_sub3 = subjets.at(2).v4().pt();
	  a_sub3  = subjets.at(2).jetArea();
	}

      // TH1Fs
      Pt_HotvrTopjets->Fill(pt_topjet, weight);
      Eta_HotvrTopjets->Fill(topjet.v4().eta(), weight);
      M_HotvrTopjets->Fill(topjet.v4().M(), weight);
      A_HotvrTopjets->Fill(a_topjet, weight);
      NSub_HotvrTopjets->Fill(subjets.size(), weight);
      Fpt_HotvrTopjets->Fill(fpt, weight);
      Mpair_HotvrTopjets->Fill(mpair, weight);
      Tau32_HotvrTopjets->Fill(topjet.tau3_groomed()/topjet.tau2_groomed(), weight);
      Pt_HotvrTopjets_Sub1->Fill(pt_sub1, weight);
      A_HotvrTopjets_Sub1->Fill(a_sub1, weight);
      Pt_HotvrTopjets_Sub2->Fill(pt_sub2, weight);
      A_HotvrTopjets_Sub2->Fill(a_sub2, weight);
      Pt_HotvrTopjets_Sub3->Fill(pt_sub3, weight);
      A_HotvrTopjets_Sub3->Fill(a_sub3, weight);

      if (muons.size() > 0)
	{
	  DeltaR_L_HotvrTopjets->Fill(deltaR(topjet.v4(), muons.at(0).v4()), weight);
	  DeltaPhi_L_HotvrTopjets->Fill(deltaPhi(topjet.v4(), muons.at(0).v4()), weight);
	}

      if (n_jets == 1)
	{
	  Pt_HotvrTopjet1->Fill(pt_topjet, weight);
	  Eta_HotvrTopjet1->Fill(topjet.v4().eta(), weight);
	  M_HotvrTopjet1->Fill(topjet.v4().M(), weight);
	  A_HotvrTopjet1->Fill(a_topjet, weight);
	  NSub_HotvrTopjet1->Fill(subjets.size(), weight);
	  Fpt_HotvrTopjet1->Fill(fpt, weight);
	  Mpair_HotvrTopjet1->Fill(mpair, weight);
	  Tau32_HotvrTopjet1->Fill(topjet.tau3_groomed()/topjet.tau2_groomed(), weight);
	  Pt_HotvrTopjet1_Sub1->Fill(pt_sub1, weight);
	  A_HotvrTopjet1_Sub1->Fill(a_sub1, weight);
	  Pt_HotvrTopjet1_Sub2->Fill(pt_sub2, weight);
	  A_HotvrTopjet1_Sub2->Fill(a_sub2, weight);
	  Pt_HotvrTopjet1_Sub3->Fill(pt_sub3, weight);
	  A_HotvrTopjet1_Sub3->Fill(a_sub3, weight);
	  if (muons.size() > 0)
	    {
	      DeltaR_L_HotvrTopjet1->Fill(deltaR(topjet.v4(), muons.at(0).v4()), weight);
	      DeltaPhi_L_HotvrTopjet1->Fill(deltaPhi(topjet.v4(), muons.at(0).v4()), weight);
	    }
	}
      
    }
  N_HotvrTopjets->Fill(n_jets, weight);

  // fill bjet hists

  const CSVBTag btag_loose(CSVBTag::WP_LOOSE);
  const CSVBTag btag_medium(CSVBTag::WP_MEDIUM);
  const CSVBTag btag_tight(CSVBTag::WP_TIGHT);

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

/*
 * WARNING: Fill this Hists only after cuts are applied to ensure
 * there is >=1 HOTVR TopJet.
 * 
 * This Hists Class implements Histograms with informations about the
 * physics performance of the HOTVR algorithm.
 * 
 */
HOTVRGenHists::HOTVRGenHists(Context & ctx, const string & dirname, const boost::optional<GenTopJetId> &gentopjetid): 
  Hists(ctx, dirname), 
  h_bstargen(ctx.get_handle<BstarToTWGen>("BstarToTWgen")),
  m_gentopjetid(gentopjetid)
{
  // book all histograms here

  // HOTVR hists
  N_HotvrTopjets          = book<TH1F>("N_HOTVR",          "N_{topjets}", 10,  0, 10);
  Pt_HotvrTopjets         = book<TH1F>("Pt_HOTVR",         "p_{t}^{topjet} [GeV/c]", 40, 0, 1600);
  Eta_HotvrTopjets        = book<TH1F>("Eta_HOTVR",        "#eta^{topjet}", 30, -6, 6);
  M_HotvrTopjets          = book<TH1F>("M_HOTVR",          "M^{topjet} [GeV/c^{2}]", 40,  0, 400);
  R_HotvrTopjets          = book<TH1F>("R_HOTVR",          "R_{topjet}", 15,  0.1, 1.6);
  NSub_HotvrTopjets       = book<TH1F>("NSub_HOTVR",       "N_{subjets}", 10,  0, 10);
  Fpt_HotvrTopjets        = book<TH1F>("Fpt_HOTVR",        "f_{pt, 1}", 20,  0, 1);
  Mpair_HotvrTopjets      = book<TH1F>("Mpair_HOTVR",      "M_pair [GeV/c^{2}]", 40,  0, 200);
  Pt_HotvrTopjets_Sub1    = book<TH1F>("Pt_HOTVR_Subjet1",         "p_{t}^{subjet1} [GeV/c]", 100, 0, 1000);
  Pt_HotvrTopjets_Sub2    = book<TH1F>("Pt_HOTVR_Subjet2",         "p_{t}^{subjet2} [GeV/c]", 100, 0, 1000);
  Pt_HotvrTopjets_Sub3    = book<TH1F>("Pt_HOTVR_Subjet3",         "p_{t}^{subjet3} [GeV/c]", 100, 0, 1000);

  DeltaR_L_HotvrTopjets   = book<TH1F>("DeltaR_L_HOTVR",   "#Delta R_{l,t}", 20, 0, 4);
  DeltaPhi_L_HotvrTopjets = book<TH1F>("DeltaPhi_L_HOTVR", "#Delta #phi_{l,t}", 20, 0, 4);

  N_HotvrTopjet1          = book<TH1F>("N_HOTVR1",          "N_{topjets}", 10,  0, 10);
  Pt_HotvrTopjet1         = book<TH1F>("Pt_HOTVR1",         "p_{t}^{topjet} [GeV/c]", 40, 0, 1600);
  Eta_HotvrTopjet1        = book<TH1F>("Eta_HOTVR1",        "#eta^{topjet}", 30, -6, 6);
  M_HotvrTopjet1          = book<TH1F>("M_HOTVR1",          "M^{topjet} [GeV/c^{2}]", 40,  0, 400);
  R_HotvrTopjet1          = book<TH1F>("R_HOTVR1",           "R_{topjet}", 15,  0.1, 1.6);
  NSub_HotvrTopjet1       = book<TH1F>("NSub_HOTVR1",       "N_{subjets}", 10,  0, 10);
  Fpt_HotvrTopjet1        = book<TH1F>("Fpt_HOTVR1",        "f_{pt, 1}", 20,  0, 1);
  Mpair_HotvrTopjet1      = book<TH1F>("Mpair_HOTVR1",      "M_pair [GeV/c^{2}]", 40,  0, 200);
  Pt_HotvrTopjet1_Sub1    = book<TH1F>("Pt_HOTVR1_Subjet1",         "p_{t}^{subjet1} [GeV/c]", 100, 0, 1000);
  Pt_HotvrTopjet1_Sub2    = book<TH1F>("Pt_HOTVR1_Subjet2",         "p_{t}^{subjet2} [GeV/c]", 100, 0, 1000);
  Pt_HotvrTopjet1_Sub3    = book<TH1F>("Pt_HOTVR1_Subjet3",         "p_{t}^{subjet3} [GeV/c]", 100, 0, 1000);

  DeltaR_L_HotvrTopjet1   = book<TH1F>("DeltaR_L_HOTVR1",   "#Delta R_{l,t}", 20, 0, 4);
  DeltaPhi_L_HotvrTopjet1 = book<TH1F>("DeltaPhi_L_HOTVR1", "#Delta #phi_{l,t}", 20, 0, 4);

}

void HOTVRGenHists::fill(const Event & event) {  

  double weight = event.weight; // event weight
  vector<GenTopJet> hotvrJets = *event.gentopjets;
  BstarToTWGen gen = event.get(h_bstargen);
  LorentzVector muon = gen.ChargedLepton();
  
  
  int n_jets = 0;

  // fill HOTVR hists
  for (GenTopJet topjet : hotvrJets)
    {
      if (m_gentopjetid)
	{
	  if (!(*m_gentopjetid)(topjet, event)) continue;
	}
      ++n_jets;
      vector<Particle> subjets = topjet.subjets();

      double pt_topjet = topjet.v4().pt();
      double R_topjet  = 600/pt_topjet;
	if (R_topjet < 0.1) R_topjet = 0.1;
	else if (R_topjet > 1.5) R_topjet = 1.5;
      double fpt = -1;
      // fpt can only be calculated if there are subjets
      if (subjets.size() >= 1)
	{
	  fpt = subjets.at(0).v4().pt() / pt_topjet;
	}
      double mpair = -1;
      double pt_sub1 = -1;
      double pt_sub2 = -1;
      double pt_sub3 = -1;

      // mpair can only be calculated if there are at least 3 subjets
      if (subjets.size() >= 3)
	{
	  double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
	  double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
	  double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
	  mpair = min(min(m12, m13), m23);
	  pt_sub1 = subjets.at(0).v4().pt();
	  pt_sub2 = subjets.at(1).v4().pt();
	  pt_sub3 = subjets.at(2).v4().pt();
	}

      // TH1Fs
      Pt_HotvrTopjets->Fill(pt_topjet, weight);
      Eta_HotvrTopjets->Fill(topjet.v4().eta(), weight);
      M_HotvrTopjets->Fill(topjet.v4().M(), weight);
      R_HotvrTopjets->Fill(R_topjet, weight);
      NSub_HotvrTopjets->Fill(subjets.size(), weight);
      Fpt_HotvrTopjets->Fill(fpt, weight);
      Mpair_HotvrTopjets->Fill(mpair, weight);
      Pt_HotvrTopjets_Sub1->Fill(pt_sub1, weight);
      Pt_HotvrTopjets_Sub2->Fill(pt_sub2, weight);
      Pt_HotvrTopjets_Sub3->Fill(pt_sub3, weight);

      DeltaR_L_HotvrTopjets->Fill(deltaR(topjet.v4(), muon), weight);
      DeltaPhi_L_HotvrTopjets->Fill(deltaPhi(topjet.v4(), muon), weight);

      if (n_jets == 1)
	{
	  Pt_HotvrTopjet1->Fill(pt_topjet, weight);
	  Eta_HotvrTopjet1->Fill(topjet.v4().eta(), weight);
	  M_HotvrTopjet1->Fill(topjet.v4().M(), weight);
	  R_HotvrTopjet1->Fill(R_topjet, weight);
	  NSub_HotvrTopjet1->Fill(subjets.size(), weight);
	  Fpt_HotvrTopjet1->Fill(fpt, weight);
	  Mpair_HotvrTopjet1->Fill(mpair, weight);
	  Pt_HotvrTopjet1_Sub1->Fill(pt_sub1, weight);
	  Pt_HotvrTopjet1_Sub2->Fill(pt_sub2, weight);
	  Pt_HotvrTopjet1_Sub3->Fill(pt_sub3, weight);

	  DeltaR_L_HotvrTopjet1->Fill(deltaR(topjet.v4(), muon), weight);
	  DeltaPhi_L_HotvrTopjet1->Fill(deltaPhi(topjet.v4(), muon), weight);
	  //   }
	}
      
    }
  N_HotvrTopjets->Fill(n_jets, weight);


}

HOTVRGenHists::~HOTVRGenHists(){}
