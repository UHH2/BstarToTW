#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

class MuonScaleFactors2018: public uhh2::AnalysisModule {
 public:
  explicit MuonScaleFactors2018(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;
  
 private:
  const int m_hlt_runnr = 316361;
  std::unique_ptr<AnalysisModule> m_sf_trigger_before, m_sf_trigger_after, m_sf_id, m_sf_iso;
};
