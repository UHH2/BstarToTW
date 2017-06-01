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

  std::vector<std::vector<TH1F*>> pt_reso;
  std::vector<std::vector<TH1F*>> event_count;
  std::vector<std::vector<TH2F*>> pt_eta;


};
