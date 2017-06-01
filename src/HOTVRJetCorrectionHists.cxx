#include "UHH2/BstarToTW/include/HOTVRJetCorrectionHists.h"

using namespace std;
using namespace uhh2;

HOTVRJetCorrectionHists::HOTVRJetCorrectionHists(uhh2::Context & ctx, const std::string & dirname): 
  Hists(ctx, dirname) {

  unsigned int no_ptbins = 6; // rows
  unsigned int no_etabins = 12; // columns
  TH1F* initial_value;
  TH2F* initial_value_2d;
  pt_reso.resize(no_ptbins, std::vector<TH1F*>(no_etabins, initial_value));
  pt_eta.resize(no_ptbins, std::vector<TH2F*>(no_etabins, initial_value_2d));
  event_count.resize(no_ptbins, std::vector<TH1F*>(no_etabins, initial_value));

  for(unsigned int pt_bin = 0; pt_bin < no_ptbins; pt_bin++){
    for(unsigned int eta_bin = 0; eta_bin < no_etabins; eta_bin++){
      std::string pt_string = std::to_string(pt_bin);
      std::string eta_string = std::to_string(eta_bin);
      pt_reso[pt_bin][eta_bin] = book<TH1F>("PtRes_" + pt_string + "_" + eta_string, "(p^{rec}_{T, jet} - p^{gen}_{T, jet}) / p^{gen}_{T, jet})", 90, -1.5, 1.5);
      pt_eta[pt_bin][eta_bin] = book<TH2F>("pt_eta_" + pt_string + "_" + eta_string, "x=pt y=eta",  50, 0, 500, 20, -4, 4);
      event_count[pt_bin][eta_bin] = book<TH1F>("Count_" + pt_string + "_" + eta_string, "a.u.", 1, 0.5, 1.5);
    }
  }

}

void HOTVRJetCorrectionHists::fill(const Event &event) {

  // get weight
  double weight = event.weight;

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

  double dR;
  double dR_temp;
  int nearest_j;
  double gen_pt;
  double gen_eta;
  double rec_pt;
  double rec_eta;
  double R;
  int pt_bin = 100;
  int eta_bin = 100;
  // do matching
  for(unsigned int i=0; i<recjets.size(); i++){
    dR = 1000;
    nearest_j = 100;
    for(unsigned int j=0; j<genjets.size(); j++){
      dR_temp = uhh2::deltaR(recjets.at(i), genjets.at(j));
      if(dR_temp < dR){
	dR = dR_temp;
	nearest_j = j;
      }
    }
    gen_pt=genjets.at(nearest_j).v4().Pt();
    gen_eta=genjets.at(nearest_j).v4().Eta();
    rec_pt=recjets.at(i).v4().Pt();
    rec_eta=recjets.at(i).v4().Eta();
    R = rec_pt/gen_pt;
    if(nearest_j != 100 && dR <= 0.2){
      // bins in pt_gen
      if(rec_pt <= 80) pt_bin = 0;
      if(rec_pt > 80 && rec_pt <= 130) pt_bin = 1;
      if(rec_pt > 130 && rec_pt <= 180) pt_bin = 2; 
      if(rec_pt > 180 && rec_pt <= 250) pt_bin = 3;
      if(rec_pt > 250 && rec_pt <= 350) pt_bin = 4;
      if(rec_pt > 350) pt_bin = 5;
      // bins in eta_gen
      if(rec_eta <= -1.5) eta_bin = 0;
      if(rec_eta > -1.5 && rec_eta <= -1.0) eta_bin = 1;
      if(rec_eta > -1.0 && rec_eta <= -0.7) eta_bin = 2; 
      if(rec_eta > -0.7 && rec_eta <= -0.4) eta_bin = 3;
      if(rec_eta > -0.4 && rec_eta <= -0.2) eta_bin = 4;
      if(rec_eta > -0.2 && rec_eta <= -0.0) eta_bin = 5;
      if(rec_eta > 0.0 && rec_eta <= 0.2) eta_bin = 6;
      if(rec_eta > 0.2 && rec_eta <= 0.4) eta_bin = 7;
      if(rec_eta > 0.4 && rec_eta <= 0.7) eta_bin = 8;
      if(rec_eta > 0.7 && rec_eta <= 1.0) eta_bin = 9;
      if(rec_eta > 1.0 && rec_eta <= 1.5) eta_bin = 10;
      if(rec_eta > 1.5) eta_bin = 11;
      if(pt_bin != 100 && eta_bin != 100){
	pt_reso[pt_bin][eta_bin]->Fill(R, weight);
	pt_eta[pt_bin][eta_bin]->Fill(rec_pt, rec_eta, weight);
	event_count[pt_bin][eta_bin]->Fill(1, weight);

      }
    }

  }

  //---------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------
 


}
