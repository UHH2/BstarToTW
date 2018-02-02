#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <vector>

class HOTVRJetCorrectionHists: public uhh2::Hists {
 public:
  HOTVRJetCorrectionHists(uhh2::Context & ctx, const std::string & dirname);
    
  virtual void fill(const uhh2::Event & event) override;

 protected:
  void fill_hists(const LorentzVector &rec, const LorentzVector &gen, double raw, double l1raw, double rho, double weight);
  void fill_in_bin(double pt_rec, double pt_gen, double eta_rec,double rho, int flag, double weight);
  std::vector<std::vector<TH1F*>> pt_reso, pt_reso_raw, pt_reso_L1;
  std::vector<TH1F*> rho_res, rho_res_raw, rho_res_L1;
  std::vector<std::vector<TH1F*>> event_count;
  std::vector<std::vector<TH2F*>> pt_eta;


};
