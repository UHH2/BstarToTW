#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

class HOTVRJetCorrector: public uhh2::AnalysisModule {
 public:
  explicit HOTVRJetCorrector(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;

  virtual ~HOTVRJetCorrector();
 private:
  bool mIsMC;
};
