#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/TopJetIds.h"

#include <TH1D.h>
#include <TFile.h>
#include <TGraph.h>

class HOTVRPileUpUncertainty : public uhh2::AnalysisModule {

 public:
  explicit HOTVRPileUpUncertainty(uhh2::Context &ctx, TString path, TString sys_direction);
  virtual bool process(uhh2::Event &event) override;

 private: 
  bool is_mc;
  std::unique_ptr<TGraph> sf_hist;
  TString m_sys_direction;
};
