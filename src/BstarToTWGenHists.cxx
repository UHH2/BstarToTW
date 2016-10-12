#include "UHH2/BstarToTW/include/BstarToTWGenHists.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace uhh2;
using namespace std;

BstarToTWGenHists::BstarToTWGenHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // Book histograms
  book<TH1F>( "M_bstar", " M^{b*} [GeV/c^{2}]; M^{b*} [GeV/c^{2}]; N_{Events} ", 30, 1000, 4000 );
  book<TH1F>( "pt_bstar", " p_{T}^{b*} [GeV/c]; p_{T}^{b*} [GeV/c]; N_{Events} ", 50, 0, 1000 );

  book<TH1F>( "pt_top", " p_{T}^{top} [GeV/c]; p_{T}^{top} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  book<TH1F>( "eta_top", " eta^{top}; eta^{top}; N_{Events} ", 20, -5, 5 );
  book<TH1F>( "phi_top", " phi^{top}; phi^{top}; N_{Events} ", 25, -M_PI, M_PI );
  book<TH1F>( "M_top", " M^{top} [GeV/c^{2}]; M^{top} [GeV/c^{2}]; N_{Events} ", 100, 150, 200 );

  book<TH1F>( "pt_W1", " p_{T}^{W_{1}} [GeV/c]; p_{T}^{W_{1}} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  book<TH1F>( "eta_W1", " eta^{W_{1}}; eta^{W_{1}}; N_{Events} ", 20, -5, 5 );
  book<TH1F>( "phi_W1", " phi^{W_{1}}; phi^{W_{1}}; N_{Events} ", 25, -M_PI, M_PI );
  book<TH1F>( "M_W1", " M^{W_{1}} [GeV/c^{2}]; M^{W_{1}} [GeV/c^{2}]; N_{Events} ", 100, 50, 100 );

  book<TH1F>( "pt_W2", " p_{T}^{W2} [GeV/c]; p_{T}^{W2} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  book<TH1F>( "eta_W2", " eta^{W2}; eta^{W2}; N_{Events} ", 20, -5, 5 );
  book<TH1F>( "phi_W2", " phi^{W2}; phi^{W2}; N_{Events} ", 25, -M_PI, M_PI );
  book<TH1F>( "M_W2", " M^{W_{2}} [GeV/c^{2}]; M^{W_{2}} [GeV/c^{2}]; N_{Events} ", 100, 50, 100 );

  book<TH1F>( "pt_b", " p_{T}^{b} [GeV/c]; p_{T}^{b} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  book<TH1F>( "eta_b", " eta^{b}; eta^{b}; N_{Events} ", 20, -5, 5 );
  book<TH1F>( "phi_b", " phi^{b}; phi^{b}; N_{Events} ", 25, -M_PI, M_PI );
  book<TH1F>( "M_b", " M^{b} [GeV/c^{2}]; M^{b} [GeV/c^{2}]; N_{Events} ", 20, 0, 10 );

  book<TH1F>( "pt_ele", " p_{T}^{e} [GeV/c]; p_{T}^{e} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  book<TH1F>( "eta_ele", " eta^{e}; eta^{e}; N_{Events} ", 20, -5, 5 );
  book<TH1F>( "phi_ele", " phi^{e}; phi^{e}; N_{Events} ", 25, -M_PI, M_PI );
  book<TH1F>( "M_ele", " M^{e} [GeV/c^{2}]; M^{e} [GeV/c^{2}]; N_{Events} ", 20, 0, 1 );

  book<TH1F>( "pt_muo", " p_{T}^{#mu} [GeV/c]; p_{T}^{#mu} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  book<TH1F>( "eta_muo", " eta^{#mu}; eta^{#mu}; N_{Events} ", 20, -5, 5 );
  book<TH1F>( "phi_muo", " phi^{#mu}; phi^{#mu}; N_{Events} ", 25, -M_PI, M_PI );
  book<TH1F>( "M_muo", " M^{#mu} [GeV/c^{2}]; M^{#mu} [GeV/c^{2}]; N_{Events} ", 20, 0, 1 );

  book<TH1F>( "pt_elenu", " p_{T}^{#nu_{e}} [GeV/c]; p_{T}^{#nu_{e}} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  book<TH1F>( "eta_elenu", " eta^{#nu_{e}}; eta^{#nu_{e}}; N_{Events} ", 20, -5, 5 );
  book<TH1F>( "phi_elenu", " phi^{#nu_{e}}; phi^{#nu_{e}}; N_{Events} ", 25, -M_PI, M_PI );

  book<TH1F>( "pt_muonu", " p_{T}^{#nu_{#mu}} [GeV/c]; p_{T}^{#nu_{#mu}} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  book<TH1F>( "eta_muonu", " eta^{#nu_{#mu}}; eta^{#nu_{#mu}}; N_{Events} ", 20, -5, 5 );
  book<TH1F>( "phi_muonu", " phi^{#nu_{#mu}}; phi^{#nu_{#mu}}; N_{Events} ", 25, -M_PI, M_PI );

  book<TH1F>( "pt_q1", " p_{T}^{q_{1}} [GeV/c]; p_{T}^{q_{1}} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  book<TH1F>( "eta_q1", " eta^{q_{1}}; eta^{q_{1}}; N_{Events} ", 20, -5, 5 );
  book<TH1F>( "phi_q1", " phi^{q_{1}}; phi^{q_{1}}; N_{Events} ", 25, -M_PI, M_PI );
  book<TH1F>( "M_q1", " M^{q_{1}} [GeV/c^{2}]; M^{q_{1}} [GeV/c^{2}]; N_{Events} ", 20, 0, 5 );

  book<TH1F>( "pt_q2", " p_{T}^{q_{2}} [GeV/c]; p_{T}^{q_{2}} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  book<TH1F>( "eta_q2", " eta^{q_{2}}; eta^{q_{2}}; N_{Events} ", 20, -5, 5 );
  book<TH1F>( "phi_q2", " phi^{q_{2}}; phi^{q_{2}}; N_{Events} ", 25, -M_PI, M_PI );
  book<TH1F>( "M_q2", " M^{q_{2}} [GeV/c^{2}]; M^{q_{2}} [GeV/c^{2}]; N_{Events} ", 20, 0, 5 );
 

  // Get container BstarToTWgen
  h_BstarToTWgen = ctx.get_handle<BstarToTWGen>("BstarToTWgen");
}

void BstarToTWGenHists::fill(const uhh2::Event & e){
  // Do not fill histograms if BstarToTWgen information has not been filled
  if(!e.is_valid(h_BstarToTWgen)){
    return;
  }

  double eventweight = e.weight;

  // Get handles from BstarToTWGen
  const auto & BstarToTWgen = e.get(h_BstarToTWgen);

  LorentzVector bstar = BstarToTWgen.bstar().v4();
  LorentzVector top = BstarToTWgen.tbstar().v4(); 
  LorentzVector W1 = BstarToTWgen.Wbstar().v4();
  LorentzVector W2 = BstarToTWgen.tW().v4();
  LorentzVector b = BstarToTWgen.tb().v4();
  LorentzVector Wdecay1 = BstarToTWgen.Wdecay1().v4();
  LorentzVector Wdecay2 = BstarToTWgen.Wdecay2().v4();
  LorentzVector tWdecay1 = BstarToTWgen.tWdecay1().v4();
  LorentzVector tWdecay2 = BstarToTWgen.tWdecay2().v4();

  // Fill histograms
  hist("M_bstar")->Fill(bstar.M(), eventweight);
  hist("pt_bstar")->Fill(bstar.Pt(), eventweight);

  hist("pt_top")->Fill(top.Pt(), eventweight);
  hist("eta_top")->Fill(top.eta(), eventweight);
  hist("phi_top")->Fill(top.phi(), eventweight);
  hist("M_top")->Fill(top.M(), eventweight);

  hist("pt_W1")->Fill(W1.Pt(), eventweight);
  hist("eta_W1")->Fill(W1.eta(), eventweight);
  hist("phi_W1")->Fill(W1.phi(), eventweight);  
  hist("M_W1")->Fill(W1.M(), eventweight);

  hist("pt_W2")->Fill(W2.Pt(), eventweight);
  hist("eta_W2")->Fill(W2.eta(), eventweight);
  hist("phi_W2")->Fill(W2.phi(), eventweight);  
  hist("M_W2")->Fill(W2.M(), eventweight);

  hist("pt_b")->Fill(b.Pt(), eventweight);
  hist("eta_b")->Fill(b.eta(), eventweight);
  hist("phi_b")->Fill(b.phi(), eventweight);  
  hist("M_b")->Fill(b.M(), eventweight);

  // Only lepton+jets
  if(BstarToTWgen.IsSemiLeptonicDecay()){
    LorentzVector lep = BstarToTWgen.ChargedLepton().v4();
    LorentzVector nu = BstarToTWgen.Neutrino().v4();
    if(abs(BstarToTWgen.ChargedLepton().pdgId()) == 11){
      hist("pt_ele")->Fill(lep.Pt(), eventweight);
      hist("eta_ele")->Fill(lep.eta(), eventweight);
      hist("phi_ele")->Fill(lep.phi(), eventweight);
      hist("M_ele")->Fill(lep.M(), eventweight);

      hist("pt_elenu")->Fill(nu.Pt(), eventweight);
      hist("eta_elenu")->Fill(nu.eta(), eventweight);
      hist("phi_elenu")->Fill(nu.phi(), eventweight);

    }
    if(abs(BstarToTWgen.ChargedLepton().pdgId()) == 13){
      hist("pt_muo")->Fill(lep.Pt(), eventweight);
      hist("eta_muo")->Fill(lep.eta(), eventweight);
      hist("phi_muo")->Fill(lep.phi(), eventweight);
      hist("M_muo")->Fill(lep.M(), eventweight);  

      hist("pt_muonu")->Fill(nu.Pt(), eventweight);
      hist("eta_muonu")->Fill(nu.eta(), eventweight);
      hist("phi_muonu")->Fill(nu.phi(), eventweight);

    }
    if(BstarToTWgen.IsPrimaryHadronicDecay()){
      hist("pt_q1")->Fill(Wdecay1.Pt(), eventweight);
      hist("eta_q1")->Fill(Wdecay1.eta(), eventweight);
      hist("phi_q1")->Fill(Wdecay1.phi(), eventweight);
      hist("M_q1")->Fill(Wdecay1.M(), eventweight); 

      hist("pt_q2")->Fill(Wdecay2.Pt(), eventweight);
      hist("eta_q2")->Fill(Wdecay2.eta(), eventweight);
      hist("phi_q2")->Fill(Wdecay2.phi(), eventweight);
      hist("M_q2")->Fill(Wdecay2.M(), eventweight); 

    }
    else{     
      hist("pt_q1")->Fill(tWdecay1.Pt(), eventweight);
      hist("eta_q1")->Fill(tWdecay1.eta(), eventweight);
      hist("phi_q1")->Fill(tWdecay1.phi(), eventweight);
      hist("M_q1")->Fill(tWdecay1.M(), eventweight); 

      hist("pt_q2")->Fill(tWdecay2.Pt(), eventweight);
      hist("eta_q2")->Fill(tWdecay2.eta(), eventweight);
      hist("phi_q2")->Fill(tWdecay2.phi(), eventweight);
      hist("M_q2")->Fill(tWdecay2.M(), eventweight); 
    }
  }
}
