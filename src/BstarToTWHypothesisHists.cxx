#include "UHH2/BstarToTW/include/BstarToTWHypothesisHists.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace uhh2;
using namespace std;

BstarToTWHypothesisHists::BstarToTWHypothesisHists(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name ): Hists(ctx, dirname){

  TString name = discriminator_name;
    if(discriminator_name=="Chi2"){
      name = "#chi^{2}";
    }
    else{
      name += " discriminator";
    }
    Discriminator =   book<TH1F>("Discriminator", name,   100, 0, 500);
    Discriminator_2 = book<TH1F>("Discriminator_2", name, 50,  0, 500);
    Discriminator_3 = book<TH1F>("Discriminator_3", name, 300, 0,  30);     

    double xbins[17] = {0, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1700, 1900, 2100, 2400, 3500};
    double xbins_fitbin[13] = {0,500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100, 2300, 2500, 3500};
    Bstar_reco_M_unbinned = book<TH1F>("Bstar_reco_M_unbinned", "M_{tW} [GeV]", 25, 100, 3500);
    Bstar_reco_M_fine = book<TH1F>("Bstar_reco_M_fine", "M_{tW} [GeV]", 100, 100, 3500);
    Bstar_reco_M_rebin = book<TH1F>("Bstar_reco_M_rebin", "M_{tW} [GeV]", 12, xbins_fitbin);
    Bstar_reco_M = book<TH1F>("Bstar_reco_M", "M_{tW} [GeV]", 16, xbins);
    
    // Bstar_reco_M = book<TH1F>("Bstar_reco_M", "M_{b*}^{reco} [GeV/c^{2}]", 50, 0, 5000);
    Bstar_reco_Pt = book<TH1F>("Bstar_reco_Pt", "p_{T, b#ast}^{reco}", 40, 0, 400);
    
    W_reco_M = book<TH1F>("W_reco_M", "M_{W}^{reco} [GeV]", 30, 0, 300);
    W_reco_Pt = book<TH1F>("W_reco_Pt", "p_{T, W}^{reco} [GeV]", 100, 0, 2000);

    Top_reco_M = book<TH1F>("Top_reco_M", "M_{t}^{reco} [GeV]", 70, 0, 700);
    Top_reco_Pt = book<TH1F>("Top_reco_Pt", "p_{T, t}^{reco} [GeV]", 100, 0, 2000);
    DeltaR_top_W = book<TH1F>("DeltaR_top_W", "#Delta R_{t,W}", 60, 0, 6);
    DeltaPhi_top_W = book<TH1F>("DeltaPhi_top_W", "#Delta #phi_{t,W}", 50, 2.5, 3.5);
    DeltaPt_top_W = book<TH1F>("DeltaPt_top_W", "#Delta p_{T, t,W} [GeV]", 80, -400, 400);
    DeltaPt_top_W_over_pt = book<TH1F>("DeltaPt_top_W_over_pt", "#Delta p_{T t,W}/p_{T t}", 50, -0.5, 0.5);
    // DeltaEta_top_W = book<TH1F>("DeltaEta_top_W", "#Delta #eta_{t,W}", 20,-5,5);

    PtTop_over_PtW = book<TH1F>("PtTop_over_PtW", "p_{T}^{jet}/p_{T}^{W}", 40, 0, 2);
    
    Discriminator_vs_M_bstar = book<TH2F>("Discriminator_vs_M_bstar", name+" vs M_{b*}^{rec}" , 50, 0, 500, 100, 0, 5000) ;
    
    h_hyps = ctx.get_handle<std::vector<BstarToTWHypothesis>>(hyps_name);
    m_name = hyps_name;
    m_discriminator_name = discriminator_name;
}


void BstarToTWHypothesisHists::fill(const uhh2::Event & e){
  std::vector<BstarToTWHypothesis> hyps = e.get(h_hyps);
  const BstarToTWHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );
  if (!hyp)
    {
      cout << "WARNING: " + m_name  +": " + m_discriminator_name + " No hypothesis was valid!" << endl;
      return;
    }
  double weight = e.weight;

  double mbstar_reco = 0;
  if((hyp->get_topjet() + hyp->get_w()).isTimelike())
    {    
      mbstar_reco = (hyp->get_topjet() + hyp->get_w()).M();
    }
  else
    {
      mbstar_reco = sqrt(-(hyp->get_topjet()+hyp->get_w()).mass2());
    }

  double ptbstar_reco = (hyp->get_topjet()+hyp->get_w()).Pt();
  Bstar_reco_M->Fill(mbstar_reco, weight);
  // if (mbstar_reco < 5000.) 
  //   {
  //     Bstar_reco_M->Fill(mbstar_reco, weight);
  //   }
  // else 
  //   {
  //     Bstar_reco_M->Fill(4990., weight);
  //   }

  Bstar_reco_M_rebin->Fill(mbstar_reco, weight);
  Bstar_reco_M_unbinned->Fill(mbstar_reco, weight);
  Bstar_reco_M_fine->Fill(mbstar_reco, weight);


  Bstar_reco_Pt->Fill (ptbstar_reco, weight);

  Discriminator->Fill(hyp->get_discriminator(m_discriminator_name) ,weight);
  Discriminator_2->Fill(hyp->get_discriminator(m_discriminator_name) ,weight);
  Discriminator_3->Fill(hyp->get_discriminator(m_discriminator_name) ,weight);

  Discriminator_vs_M_bstar->Fill(hyp->get_discriminator(m_discriminator_name), mbstar_reco, weight);

  double mtop = 0;
  if(hyp->get_topjet().isTimelike()) mtop = hyp->get_topjet().M();
  Top_reco_M->Fill(mtop, weight);
  Top_reco_Pt->Fill(hyp->get_topjet().Pt(), weight);

  double mw = 0;
  if(hyp->get_w().isTimelike()) mw = hyp->get_w().M();
  W_reco_M->Fill(mw, weight);
  W_reco_Pt->Fill(hyp->get_w().Pt(), weight);


  double delPhi = (hyp->get_w().Phi() - hyp->get_topjet().Phi());
  if (delPhi < 0) delPhi += 2*M_PI;

  DeltaR_top_W->Fill(deltaR(hyp->get_w(), hyp->get_topjet()),weight);
  DeltaPhi_top_W->Fill(abs(hyp->get_w().Phi() - hyp->get_topjet().Phi()), weight);
  DeltaPt_top_W->Fill((hyp->get_topjet().Pt() - hyp->get_w().Pt()), weight);
  DeltaPt_top_W_over_pt->Fill((hyp->get_topjet().Pt() - hyp->get_w().Pt()) / hyp->get_topjet().Pt(), weight);

  PtTop_over_PtW->Fill(hyp->get_topjet().Pt() / hyp->get_w().Pt(), weight);

}

LeptonicTopHypothesisHists::LeptonicTopHypothesisHists(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name ): Hists(ctx, dirname){

  TString name = discriminator_name;
    if(discriminator_name=="Chi2"){
      name = "#chi^{2}";
    }
    else{
      name += " discriminator";
    }
    Discriminator =   book<TH1F>("Discriminator", name,   100, 0, 500);
    Discriminator_2 = book<TH1F>("Discriminator_2", name, 50,  0, 500);
    Discriminator_3 = book<TH1F>("Discriminator_3", name, 300, 0,  30);     

    Top_reco_M = book<TH1F>("toplep_M", "M_{t}^{reco} [GeV]", 70, 0, 700);
    Top_reco_pt = book<TH1F>("toplep_pt", "p_{T, t}^{reco} [GeV]", 100, 0, 2000);
    
    h_hyps = ctx.get_handle<std::vector<LeptonicTopHypothesis>>(hyps_name);
    m_name = hyps_name;
    m_discriminator_name = discriminator_name;
}


void LeptonicTopHypothesisHists::fill(const uhh2::Event & event){
  std::vector<LeptonicTopHypothesis> hyps = event.get(h_hyps);
  const LeptonicTopHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );
  if (!hyp)
    {
      cout << "WARNING: " + m_name  +": " + m_discriminator_name + " No hypothesis was valid!" << endl;
      return;
    }
  double weight = event.weight;
  
  Discriminator->Fill(hyp->get_discriminator(m_discriminator_name), weight);
  Discriminator_2->Fill(hyp->get_discriminator(m_discriminator_name), weight);
  Discriminator_3->Fill(hyp->get_discriminator(m_discriminator_name), weight);

  double mtop = (hyp->get_toplep().isTimelike()) ? hyp->get_toplep().M() : sqrt(-(hyp->get_toplep()).mass2());
  Top_reco_M->Fill(mtop, weight);
  Top_reco_pt->Fill(hyp->get_toplep().pt(), weight);
}
