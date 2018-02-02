#include "UHH2/BstarToTW/include/HOTVRJetCorrectionHists.h"
#include "UHH2/BstarToTW/include/myutils.h"

using namespace std;
using namespace uhh2;

HOTVRJetCorrectionHists::HOTVRJetCorrectionHists(uhh2::Context & ctx, const std::string & dirname): 
  Hists(ctx, dirname) {

  unsigned int no_ptbins = 6; // rows
  unsigned int no_etabins = 12; // columns
  unsigned int no_rhobins = 12; // rho bins
  TH1F* initial_value;
  TH2F* initial_value_2d;
  pt_reso_raw.resize(no_ptbins, std::vector<TH1F*>(no_etabins, initial_value));
  rho_res_raw.resize(no_rhobins);
  pt_reso_L1.resize(no_ptbins, std::vector<TH1F*>(no_etabins, initial_value));
  rho_res_L1.resize(no_rhobins);
  pt_reso.resize(no_ptbins, std::vector<TH1F*>(no_etabins, initial_value));
  rho_res.resize(no_rhobins);
  pt_eta.resize(no_ptbins, std::vector<TH2F*>(no_etabins, initial_value_2d));
  event_count.resize(no_ptbins, std::vector<TH1F*>(no_etabins, initial_value));

  for(unsigned int rho_bin = 0; rho_bin < no_rhobins; ++rho_bin)
    {
      std::string rho_string = std::to_string(rho_bin);
      rho_res_raw[rho_bin] = book<TH1F>("RhoResRaw_" + rho_string, "#rho", 90, -0.5, 2.5);
      rho_res_L1[rho_bin] = book<TH1F>("RhoResL1_" + rho_string, "#rho", 90, -0.5, 2.5);
      rho_res[rho_bin] = book<TH1F>("RhoRes_" + rho_string, "#rho", 90, -0.5, 2.5);
    }

  for(unsigned int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(unsigned int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      std::string pt_string = std::to_string(pt_bin);
      std::string eta_string = std::to_string(eta_bin);
      pt_reso_raw[pt_bin][eta_bin] = book<TH1F>("PtResRaw_" + pt_string + "_" + eta_string, "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet})", 90, -0.5, 2.5);
      pt_reso_L1[pt_bin][eta_bin] = book<TH1F>("PtResL1_" + pt_string + "_" + eta_string, "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet})", 90, -0.5, 2.5);
      pt_reso[pt_bin][eta_bin] = book<TH1F>("PtRes_" + pt_string + "_" + eta_string, "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet})", 90, -0.5, 2.5);
      pt_eta[pt_bin][eta_bin] = book<TH2F>("pt_eta_" + pt_string + "_" + eta_string, "x=pt y=eta",  50, 0, 500, 20, -4, 4);
      event_count[pt_bin][eta_bin] = book<TH1F>("Count_" + pt_string + "_" + eta_string, "a.u.", 1, 0.5, 1.5);
    }
  }

}

void HOTVRJetCorrectionHists::fill(const Event &event) {

  // get weight
  double weight = event.weight;
  double rho = event.rho;
  vector<TopJet> recjets = *event.topjets;
  vector<GenTopJet> genjets = *event.gentopjets;


 /* ******************************************************************************
     matching between gen and reco jets:
     - a rec jet is called isolated if the next jet is not within 2R. An isolated jet should be spherical and more simelar to an ak4 jet
     - for each reco jet, calc distance to all other gen jets
     - gen jet with lowest distance is a match if distance is < 0.2
     - then calculate resolution with reco jet and matched gen jet
     - to do: account for double counting
     ********************************************************************************* */

  for (TopJet &recjet : recjets)
    {
      // match fatjets
      const GenTopJet *matched_genjet = match(recjet, genjets, 0.3);
      if(matched_genjet) // match subjets
	{
	  for (const Jet &recsubj : recjet.subjets())
	    {
	      const Particle *matched_gensubj = match(recsubj, matched_genjet->subjets(), 0.2);
	      if (matched_gensubj)
		{
		  fill_hists(recsubj.v4(), matched_gensubj->v4(), recsubj.JEC_L1factor_raw(), recsubj.JEC_factor_raw(), rho, weight);
		}
	    }
	}
    }
}

void HOTVRJetCorrectionHists::fill_hists(const LorentzVector &rec, const LorentzVector &gen, double raw, double l1raw, double rho, double weight) {

  double pt_rec  = rec.Pt();
  double eta_rec = rec.Eta();
  double pt_gen  = gen.Pt();
  fill_in_bin(pt_rec * raw, pt_gen, eta_rec, rho, 1, weight);   // fill raw
  fill_in_bin(pt_rec * raw / l1raw, pt_gen, eta_rec, rho, 2, weight); // fill L1
  fill_in_bin(pt_rec, pt_gen, eta_rec, rho, 3, weight);         // fill corr

}

void HOTVRJetCorrectionHists::fill_in_bin(double pt_rec, double pt_gen, double eta_rec, double rho, int flag, double weight)
{
  /* Flag: 1 - raw
           2 - L1
	   3 - corr
  */

  int pt_bin = 100;
  int eta_bin = 100;
  int rho_bin = 100;  

  double R = pt_rec/pt_gen; 

  // bins in pt_gen
  if(pt_rec <= 80) pt_bin = 0;
  else if(pt_rec > 80 && pt_rec <= 130) pt_bin = 1;
  else if(pt_rec > 130 && pt_rec <= 180) pt_bin = 2; 
  else if(pt_rec > 180 && pt_rec <= 250) pt_bin = 3;
  else if(pt_rec > 250 && pt_rec <= 350) pt_bin = 4;
  else if(pt_rec > 350) pt_bin = 5;
  // bins in eta_gen
  if(eta_rec <= -1.5) eta_bin = 0;
  else if(eta_rec > -1.5 && eta_rec <= -1.0) eta_bin = 1;
  else if(eta_rec > -1.0 && eta_rec <= -0.7) eta_bin = 2; 
  else if(eta_rec > -0.7 && eta_rec <= -0.4) eta_bin = 3;
  else if(eta_rec > -0.4 && eta_rec <= -0.2) eta_bin = 4;
  else if(eta_rec > -0.2 && eta_rec <= -0.0) eta_bin = 5;
  else if(eta_rec > 0.0 && eta_rec <= 0.2) eta_bin = 6;
  else if(eta_rec > 0.2 && eta_rec <= 0.4) eta_bin = 7;
  else if(eta_rec > 0.4 && eta_rec <= 0.7) eta_bin = 8;
  else if(eta_rec > 0.7 && eta_rec <= 1.0) eta_bin = 9;
  else if(eta_rec > 1.0 && eta_rec <= 1.5) eta_bin = 10;
  else if(eta_rec > 1.5) eta_bin = 11;
  if (rho < 5) rho_bin = 0;
  else if(rho < 10) rho_bin = 1;
  else if(rho < 15) rho_bin = 2;
  else if(rho < 20) rho_bin = 3;
  else if(rho < 25) rho_bin = 4;
  else if(rho < 30) rho_bin = 5;
  else if(rho < 35) rho_bin = 6;
  else if(rho < 40) rho_bin = 7;
  else if(rho < 45) rho_bin = 8;
  else if(rho < 50) rho_bin = 9;
  else if(rho < 55) rho_bin = 10;
  else if(rho >= 55) rho_bin = 11;
  if(pt_bin != 100 && eta_bin != 100)
    {
      if (flag == 1) 
	{
	  pt_reso_raw[pt_bin][eta_bin]->Fill(R, weight);
	  rho_res_raw[rho_bin]->Fill(R, weight);
	}
      else if (flag ==2) 
	{
	  pt_reso_L1[pt_bin][eta_bin]->Fill(R, weight);
	  rho_res_L1[rho_bin]->Fill(R, weight);
	}
      else if (flag ==3)
	{
	  pt_reso[pt_bin][eta_bin]->Fill(R, weight);
	  rho_res[rho_bin]->Fill(R, weight);
	  pt_eta[pt_bin][eta_bin]->Fill(pt_rec, eta_rec, weight);
	  event_count[pt_bin][eta_bin]->Fill(1, weight);
	}
    } 
}
