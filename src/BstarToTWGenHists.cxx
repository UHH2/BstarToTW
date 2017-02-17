#include "UHH2/BstarToTW/include/BstarToTWGenHists.h"
#include "TH1F.h"
#include "TH2F.h"
#include <vector>

using namespace uhh2;
using namespace std;

BstarToTWGenHists::BstarToTWGenHists(uhh2::Context & ctx, const std::string & dirname): 
  Hists(ctx, dirname),
  h_BstarToTWgen(ctx.get_handle<BstarToTWGen>("BstarToTWgen"))
{
  // --------------------------------------------------------------
  // Book histograms
  M_bstar  = book<TH1F>( "M_bstar", " M^{b*} [GeV/c^{2}]; M^{b*} [GeV/c^{2}]; N_{Events} ", 100, 0, 4000 );
  pt_bstar = book<TH1F>( "pt_bstar", " p_{T}^{b*} [GeV/c]; p_{T}^{b*} [GeV/c]; N_{Events} ", 50, 0, 1000 );

  pt_top   = book<TH1F>( "pt_top", " p_{T}^{top} [GeV/c]; p_{T}^{top} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  eta_top  = book<TH1F>( "eta_top", " #eta^{top}; #eta^{top}; N_{Events} ", 20, -5, 5 );
  phi_top  = book<TH1F>( "phi_top", " #phi^{top}; #phi^{top}; N_{Events} ", 25, -M_PI, M_PI );
  M_top    = book<TH1F>( "M_top", " M^{top} [GeV/c^{2}]; M^{top} [GeV/c^{2}]; N_{Events} ", 100, 150, 200 );
  M_tophad = book<TH1F>( "M_tophad", " M^{top,had} [GeV/c^{2}]; M^{top,had} [GeV/c^{2}]; N_{Events} ", 100, 150, 200 );
  M_toplep = book<TH1F>( "M_toplep", " M^{top,lep} [GeV/c^{2}]; M^{top,lep} [GeV/c^{2}]; N_{Events} ", 100, 150, 200 );

  pt_W1  = book<TH1F>( "pt_W1", " p_{T}^{W_{1}} [GeV/c]; p_{T}^{W_{1}} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  eta_W1 = book<TH1F>( "eta_W1", " #eta^{W_{1}}; #eta^{W_{1}}; N_{Events} ", 20, -5, 5 );
  phi_W1 = book<TH1F>( "phi_W1", " #phi^{W_{1}}; #phi^{W_{1}}; N_{Events} ", 25, -M_PI, M_PI );
  M_W1   = book<TH1F>( "M_W1", " M^{W_{1}} [GeV/c^{2}]; M^{W_{1}} [GeV/c^{2}]; N_{Events} ", 100, 50, 100 );

  pt_W2  = book<TH1F>( "pt_W2", " p_{T}^{W2} [GeV/c]; p_{T}^{W2} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  eta_W2 = book<TH1F>( "eta_W2", " #eta^{W2}; #eta^{W2}; N_{Events} ", 20, -5, 5 );
  phi_W2 = book<TH1F>( "phi_W2", " #phi^{W2}; #phi^{W2}; N_{Events} ", 25, -M_PI, M_PI );
  M_W2   = book<TH1F>( "M_W2", " M^{W_{2}} [GeV/c^{2}]; M^{W_{2}} [GeV/c^{2}]; N_{Events} ", 100, 50, 100 );

  pt_b  = book<TH1F>( "pt_b", " p_{T}^{b} [GeV/c]; p_{T}^{b} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  eta_b = book<TH1F>( "eta_b", " #eta^{b}; #eta^{b}; N_{Events} ", 20, -5, 5 );
  phi_b = book<TH1F>( "phi_b", " #phi^{b}; #phi^{b}; N_{Events} ", 25, -M_PI, M_PI );
  M_b   = book<TH1F>( "M_b", " M^{b} [GeV/c^{2}]; M^{b} [GeV/c^{2}]; N_{Events} ", 20, 0, 10 );

  pt_ele  = book<TH1F>( "pt_ele", " p_{T}^{e} [GeV/c]; p_{T}^{e} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  eta_ele = book<TH1F>( "eta_ele", " #eta^{e}; #eta^{e}; N_{Events} ", 20, -5, 5 );
  phi_ele = book<TH1F>( "phi_ele", " #phi^{e}; #phi^{e}; N_{Events} ", 25, -M_PI, M_PI );

  pt_muo  = book<TH1F>( "pt_muo", " p_{T}^{#mu} [GeV/c]; p_{T}^{#mu} [GeV/c]; N_{Events} ", 50, 0, 2000 );
  eta_muo = book<TH1F>( "eta_muo", " #eta^{#mu}; #eta^{#mu}; N_{Events} ", 20, -5, 5 );
  phi_muo = book<TH1F>( "phi_muo", " #phi^{#mu}; #phi^{#mu}; N_{Events} ", 25, -M_PI, M_PI );

  cosThetastar_bstar = book<TH1F>( "cosThetastar_bstar", "cos(#theta *_{b*})", 20, -1, 1);
  cosThetastar_t = book<TH1F>( "cosThetastar_t", "cos(#theta *_{t})", 20, -1, 1);
  deltaRMax = book<TH1F>( "deltaRMax", "#Delta R_{max}", 20, 0, 2);
}

void BstarToTWGenHists::fill(const uhh2::Event & e)
{
  // Do not fill histograms if BstarToTWgen information has not been filled
  if(!e.is_valid(h_BstarToTWgen))
    {
      return;
    }

  double eventweight = e.weight;
  const auto & BstarToTWgen = e.get(h_BstarToTWgen);

  // --------------------------------------------------------------
  // Get genparticles
  vector<LorentzVector> genStable = BstarToTWgen.stable();

  LorentzVector bstar = BstarToTWgen.bstar();
  LorentzVector top = BstarToTWgen.tbstar();
  LorentzVector W1 = BstarToTWgen.Wbstar();
  LorentzVector W2 = BstarToTWgen.tW();
  LorentzVector b = BstarToTWgen.tb();
  //LorentzVector Wdecay1 = BstarToTWgen.Wdecay1();
  //LorentzVector Wdecay2 = BstarToTWgen.Wdecay2();
  LorentzVector tWdecay1 = BstarToTWgen.tWdecay1();
  LorentzVector tWdecay2 = BstarToTWgen.tWdecay2();

  TLorentzVector tbstar(bstar.X(), bstar.Y(), bstar.Z(), bstar.T());
  TLorentzVector ttop(top.X(), top.Y(), top.Z(), top.T());
  TLorentzVector tW2(W2.X(), W2.Y(), W2.Z(), W2.T());



  // --------------------------------------------------------------
  // Fill histograms
  // Bstar
  M_bstar->Fill( bstar.M(), eventweight );
  pt_bstar->Fill( bstar.Pt(), eventweight );

  // Top
  pt_top->Fill( top.Pt(), eventweight );
  eta_top->Fill( top.eta(), eventweight );
  phi_top->Fill( top.phi(), eventweight );
  M_top->Fill( top.M(), eventweight );

  // W_bstar
  pt_W1->Fill( W1.Pt(), eventweight );
  eta_W1->Fill( W1.eta(), eventweight );
  phi_W1->Fill( W1.phi(), eventweight );  
  M_W1->Fill( W1.M(), eventweight );

  // W_top
  pt_W2->Fill( W2.Pt(), eventweight );
  eta_W2->Fill( W2.eta(), eventweight );
  phi_W2->Fill( W2.phi(), eventweight );  
  M_W2->Fill( W2.M(), eventweight );

  // b_top
  pt_b->Fill( b.Pt(), eventweight );
  eta_b->Fill( b.eta(), eventweight );
  phi_b->Fill( b.phi(), eventweight );  
  M_b->Fill( b.M(), eventweight );
 
  // Only lepton+jets
  if(BstarToTWgen.IsSemiLeptonicDecay())
    {
      if(BstarToTWgen.IsWHadronicDecay())
	{
      
	  M_toplep->Fill( ( tWdecay1 + tWdecay2 + b ).M(), eventweight );
	}
      else if(BstarToTWgen.IsTopHadronicDecay())
	{     
	  LorentzVector lep = BstarToTWgen.ChargedLepton();
	  if(BstarToTWgen.IsElectronDecay())
	    {
	      pt_ele->Fill( lep.Pt(), eventweight );
	      eta_ele->Fill( lep.eta(), eventweight );
	      phi_ele->Fill( lep.phi(), eventweight );
	    }

	  if(BstarToTWgen.IsMuonDecay())
	    {
	      pt_muo->Fill( lep.Pt(), eventweight );
	      eta_muo->Fill( lep.eta(), eventweight );
	      phi_muo->Fill( lep.phi(), eventweight );
	    }
	  M_tophad->Fill( ( tWdecay1 + tWdecay2 + b ).M(), eventweight );

	  // get theta*_b* and theta*_t
	  tW2.Boost(-ttop.BoostVector());
	  double thetastar_t = tW2.Angle(ttop.Vect());

	  ttop.Boost(-tbstar.BoostVector());
	  double thetastar_bstar = ttop.Angle(tbstar.Vect());

	  cosThetastar_bstar->Fill( cos(thetastar_bstar), eventweight );
	  cosThetastar_t->Fill( cos(thetastar_t), eventweight );

	  double dRMax = max( max(deltaR(b, tWdecay1), deltaR(b, tWdecay2)), deltaR(tWdecay1, tWdecay2));
	  deltaRMax->Fill(dRMax, eventweight);
	}

    }
}
