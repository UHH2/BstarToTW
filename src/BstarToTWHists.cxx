#include "UHH2/core/include/Utils.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/BstarToTW/include/BstarToTWHists.h"


#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>


using namespace std;
using namespace uhh2;


bool operator<(const run_lumi & rl1, const run_lumi & rl2){
    if(rl1.run == rl2.run){
        return rl1.lumiblock < rl2.lumiblock;
    }
    else{
        return rl1.run < rl2.run;
    }
}

BstarToTWHists::BstarToTWHists(Context & ctx, const string & dirname): 
  Hists(ctx, dirname) {

  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");

  // book all histograms here
  MET = book<TH1F>("MET", "#slash{E}_{T} [GeV/c]", 20, 0, 1000);
  HT_lep = book<TH1F>("HT_lep", "HT_{lep} [GeV/c]", 70, 0, 3500);
  HT_jet = book<TH1F>("HT_jet", "HT_{jet} [GeV/c]", 70, 0, 3500);
  ST = book<TH1F>("ST", "S_{T} [GeV/c]", 80, 0, 4000);
  rho = book<TH1F>("rho", "#rho [GeV/c]", 12, 0, 65);
  deltaPhi_blep = book<TH1F>("deltaPhi_blep", "#Delta #phi_{b,l}", 35, 0, 3.5);
  deltaPhi_btop = book<TH1F>("deltaPhi_btop", "#Delta #phi_{b,t}", 35, 0, 3.5);
  deltaPhi_hotvr_ak4 = book<TH1F>("deltaPhi_hotvr_ak4", "#Delta #phi_{hotvr,ak4}", 35, 0, 3.5);
  deltaPhi_lep_ak4 = book<TH1F>("deltaPhi_lep_ak4", "#Delta #phi_{l,ak4}", 35, 0, 3.5);
  deltaPhi_hotvr_lak4 = book<TH1F>("deltaPhi_hotvr_leadingak4", "#Delta #phi_{hotvr,ak4}", 35, 0, 3.5);
  deltaPhi_lep_lak4 = book<TH1F>("deltaPhi_lep_leadingak4", "#Delta #phi_{l,ak4}", 35, 0, 3.5);
  // deltaRmin_vs_pTrel = book <TH2F>("deltaRmin_vs_pTrel", "2D_cut", );
}


void BstarToTWHists::fill(const Event & event){
  double weight = event.weight;
  vector<Jet> jets = *event.jets;
  vector<Electron> electrons = *event.electrons;
  vector<Muon> muons = *event.muons;
  vector<TopJet> topjets = *event.topjets;
  
  
  double ht_lep = 0;
  double ht_jet = 0;

  const Particle &primlep = event.get(h_primlep);

  for (Electron ele : electrons)
    {
      ht_lep += ele.v4().pt();
    }

  for (Muon muo : muons)
    {
      ht_lep += muo.v4().pt();
    }
  MET->Fill(event.met->pt(), weight);
  ht_lep += event.met->pt();
  int i = 0;
  for (Jet jet : jets)
    {
      ht_jet += jet.v4().pt();
      if (topjets.size() > 0)
	deltaPhi_hotvr_ak4->Fill(deltaPhi(jet.v4(),topjets.at(0).v4()), weight);
      deltaPhi_lep_ak4->Fill(deltaPhi(jet.v4(),primlep.v4()), weight);
      if(btag_loose(jet, event))
	{
	  if (i == 0)
	    deltaPhi_blep->Fill(deltaPhi(jet.v4(),primlep.v4()), weight);
	  for (TopJet topjet : topjets)
	    {
	      deltaPhi_btop->Fill(deltaPhi(jet.v4(),topjet.v4()), weight);
	    }
	  ++i;
	}
    }

  if (jets.size() > 0)
    {
      if (topjets.size() > 0)
	deltaPhi_hotvr_lak4->Fill(deltaPhi(jets.at(0).v4(),topjets.at(0).v4()), weight);
      deltaPhi_lep_lak4->Fill(deltaPhi(jets.at(0).v4(),primlep.v4()), weight);
    }

  HT_lep->Fill(ht_lep, weight);
  HT_jet->Fill(ht_jet, weight);
  ST->Fill(ht_lep + ht_jet, weight);
  
  rho->Fill(event.rho, weight);
  
}

BstarToTWHists::~BstarToTWHists(){}

BstarToTWBackgroundHists::BstarToTWBackgroundHists(Context & ctx, const string & dirname, const string & hyps_name, const TString & path):
  Hists(ctx, dirname){
  h_hyps = ctx.get_handle<std::vector<BstarToTWHypothesis>>(hyps_name);
  TFile* f = new TFile(path);

  TH1F* ha = (TH1F*)(f->Get("ratio"))->Clone("ratio");
  // TH1F* ha = (TH1F*)(f->Get("signal"))->Clone("signal");
  TF1* fa = ha->GetFunction("fitfun");
  double* pa = fa->GetParameters();
  int nparam = 3;
  m_pa = new double[nparam];
  for (int i = 0; i<nparam; ++i)
    {
      m_pa[i] = pa[i];
    }
  // TH1F* hb = (TH1F*)(f->Get("control"))->Clone("control");
  // TF1* fb = hb->GetFunction("fitfun");  
  // double* pb = fb->GetParameters();
  // m_pb = new double[nparam];
  // for (int i = 0; i<nparam; ++i)
  //   {
  //     m_pb[i] = pb[i];
  //   }
  
  f->Close();

  double xbins[17] = {0, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1700, 1900, 2100, 2400, 5000};
  double xbins_fitbin[13] = {0,500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100, 2300, 2500, 5000};
  Bstar_reco_M_unbinned = book<TH1F>("Bstar_reco_M_unbinned", "M_{tW} [GeV/c^{2}]", 25, 100, 5100);
  Bstar_reco_M_rebin = book<TH1F>("Bstar_reco_M_rebin", "M_{tW} [GeV/c^{2}]", 12, xbins_fitbin);
  Bstar_reco_M = book<TH1F>("Bstar_reco_M", "M_{tW} [GeV/c^{2}]", 16, xbins);
}

void BstarToTWBackgroundHists::fill(const Event & event){
  
  std::vector<BstarToTWHypothesis> hyps = event.get(h_hyps);
  const BstarToTWHypothesis* hyp = get_best_hypothesis( hyps, "Chi2" );
  if (!hyp)
    {
      //cout << "WARNING: " + m_hyps_name  +": " + m_discriminator_name + " No hypothesis was valid!" << endl;
      return;
    }
  
  double mbstar = 0;
  if((hyp->get_topjet() + hyp->get_w()).isTimelike())
    {    
      mbstar = (hyp->get_topjet() + hyp->get_w()).M();
    }
  else
    {
      mbstar = sqrt(-(hyp->get_topjet()+hyp->get_w()).mass2());
    }

  // double background_weight = m_pa[0] + m_pa[1] * mbstar;
  double background_weight = exp(m_pa[0] + m_pa[1] * mbstar) + m_pa[2];
  
  // double weight_a = (1.+erf((mbstar - m_pa[3])/m_pa[4]))/2. * exp(m_pa[0] + m_pa[1] * mbstar + m_pa[2] * mbstar * mbstar);
  // double weight_b = (1.+erf((mbstar - m_pb[3])/m_pb[4]))/2. * exp(m_pb[0] + m_pb[1] * mbstar + m_pb[2] * mbstar * mbstar);
  // double background_weight = weight_a / weight_b;

  // cout << "weight a: " << weight_a << endl;
  // cout << "weight b: " << weight_b << endl;
  // cout << "Background weight: " << background_weight << endl;
  if ( background_weight < 0) 
    background_weight = 0;
  const double event_weight = event.weight;
  background_weight *= event_weight;
  if (mbstar < 5000.) 
    {
      Bstar_reco_M->Fill(mbstar, background_weight);
      Bstar_reco_M_rebin->Fill(mbstar, background_weight);
      Bstar_reco_M_unbinned->Fill(mbstar, background_weight);
    }
  else 
    {
      Bstar_reco_M->Fill(4999., background_weight);
      Bstar_reco_M_rebin->Fill(4999., background_weight);
      Bstar_reco_M_unbinned->Fill(4999., background_weight);
    }

}

BstarToTWBackgroundHists::~BstarToTWBackgroundHists(){}

BstarToTWAnalysisHists::BstarToTWAnalysisHists(uhh2::Context &ctx, const std::string &dirname):
  Hists(ctx,dirname) {
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_toptag = ctx.get_handle<vector<TopJet> >("toptag");
  h_btag = ctx.get_handle<vector<Jet> >("btag_loose");
  h_ht = ctx.get_handle<double>("HT");

  // --- Book histograms
  // Event variables in nicer binning
  pt_lep = book<TH1F>("lepton_pt","p_{T, lep} [GeV]", 40, 0, 2000);
  event_met = book<TH1F>("MET", "miss. E_{T} [GeV]", 40, 0, 2000);
  event_ht = book<TH1F>("HT", "H_{T} [GeV]", 50, 0, 2500);
  event_htlep = book<TH1F>("HTlep", "H_{T, lep} [GeV]", 50, 0, 2500);
  event_st = book<TH1F>("ST", "S_{T} [GeV]", 50, 0, 5000);
  event_st_half = book<TH1F>("ST_half", "S_{T} [GeV]", 25, 0, 5000);
  event_st_fifth = book<TH1F>("ST_fifth", "S_{T} [GeV]", 10, 0, 5000);
  event_rho = book<TH1F>("rho", "#rho [GeV]", 20, 0, 100);

  ratio_ht = book<TH1F>("HTratio", "ratio H_{T}", 20, 0, 4);
  asymm_st = book<TH1F>("STasymm", "asymm. S_{T}", 20, -1, 1);

  lep_isolation = book<TH1F>("lepton_isolation", "rel. iso", 50, 0, 0.5);

  // Angular distributions
  deltaPhi_blep = book<TH1F>("deltaPhi_blep", "#Delta #phi_{b,l}", 35, 0, 3.5);
  deltaPhi_btop = book<TH1F>("deltaPhi_btop", "#Delta #phi_{b,t}", 35, 0, 3.5);
  deltaPhi_toplep = book<TH1F>("deltaPhi_toplep", "#Delta #phi_{t,l}", 35, 0, 3.5);
  deltaPhi_lepmet = book<TH1F>("deltaPhi_lepmet", "#Delta #phi_{l,MET}", 35, 0, 3.5);
  deltaR_blep = book<TH1F>("deltaR_blep", "#Delta R_{b,l}", 60, 0, 6);
  deltaR_btop = book<TH1F>("deltaR_btop", "#Delta R_{b,t}", 60, 0, 6);
  deltaR_toplep = book<TH1F>("deltaR_toplep", "#Delta R_{t,l}", 60, 0, 6);

}

void BstarToTWAnalysisHists::fill(const Event &event) {
  // check if handles were set, otherwise return
  if (!event.is_valid(h_primlep) || !event.is_valid(h_toptag) || !event.is_valid(h_btag) || !event.is_valid(h_ht))
    return;

  double weight = event.weight;
  FlavorParticle lep = event.get(h_primlep);
  vector<TopJet> topjets = event.get(h_toptag);
  vector<Jet> bjets = event.get(h_btag);
  double met = event.met->pt();
  double ht = event.get(h_ht);
  double htlep = lep.pt() + met;
  double st = ht + lep.pt() + met;
  
  pt_lep->Fill(lep.pt(), weight);
  event_met->Fill(met, weight);
  event_ht->Fill(ht, weight);
  event_htlep->Fill(htlep, weight);
  event_st->Fill(st, weight);
  event_st_half->Fill(st, weight);
  event_st_fifth->Fill(st, weight);
  event_rho->Fill(event.rho, weight);

  ratio_ht->Fill(ht/htlep, weight);
  asymm_st->Fill((ht-htlep)/st, weight);
  
  for (TopJet &topjet : topjets) 
    {
      deltaPhi_toplep->Fill(deltaPhi(topjet.v4(), lep.v4()), weight);
      deltaR_toplep->Fill(deltaR(topjet.v4(), lep.v4()), weight);
      for (Jet &bjet : bjets)
	{
	  deltaPhi_btop->Fill(deltaPhi(topjet.v4(), bjet.v4()), weight);
	  deltaR_btop->Fill(deltaR(topjet.v4(), bjet.v4()), weight);
	}
    }
  for (Jet &bjet : bjets)
    {
      deltaPhi_blep->Fill(deltaPhi(lep.v4(), bjet.v4()), weight);
      deltaR_blep->Fill(deltaR(lep.v4(), bjet.v4()), weight);
    }

  deltaPhi_lepmet->Fill(deltaPhi(lep.v4(),event.met->v4()), weight);

  for (const Electron &ele : *event.electrons)
    {
      lep_isolation->Fill(ele.relIsorho(event.rho), weight);
    }

  for (const Muon &muo : *event.muons)
    {
      lep_isolation->Fill(muo.relIso(), weight);
    }
 
}

BstarToTWAnalysisHists::~BstarToTWAnalysisHists(){}

void TopMatchHists::fill(const Event &event) {
  
  if (event.isRealData) return;
  
  const vector<GenTop> &gentops = event.get(h_tophad);
  const vector<TopJet> &topjets = *event.topjets;

  for (const TopJet &topjet : topjets)
    {
      
    }
}

TopMatchHists::~TopMatchHists(){}
