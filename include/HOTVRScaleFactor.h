#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/TopJetIds.h"

#include <TH1D.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>

class HOTVRScaleFactor : public uhh2::AnalysisModule {
 public:
  explicit HOTVRScaleFactor(uhh2::Context &ctx, std::string signal_name, TopJetId id_topjet,  TString path, TString sys_direction);
  virtual bool process(uhh2::Event &event) override;

 private:
  bool m_do_weight;
  TopJetId m_id_topjet;
  std::unique_ptr<TGraphAsymmErrors> sf_hist;
  TString m_sys_direction;
};

class CMSTTScaleFactor : public uhh2::AnalysisModule {

 public:

  explicit CMSTTScaleFactor(uhh2::Context &ctx, std::string signal_name, TString sys_direction);
  virtual bool process(uhh2::Event &event) override;

 private:

  bool m_do_weight;
  TString m_sys_direction;
};
