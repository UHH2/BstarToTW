#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"

#include <vector>

namespace uhh2 {

  class EfficiencyHists: public uhh2::Hists {
  public:
    EfficiencyHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~EfficiencyHists();
  protected:
    
    TH1F *hotvr_counter, *cms_counter;

    uhh2::Event::Handle<std::vector<TopJet>> h_AK8Jets;
  };
   
}
