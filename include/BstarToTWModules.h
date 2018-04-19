#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include <TFile.h>
#include <TGraphAsymmErrors.h>


class ElectronTriggerWeights: public uhh2::AnalysisModule{

 public:
  explicit ElectronTriggerWeights(uhh2::Context & ctx, TString path_, TString SysDirection_);
  virtual bool process(uhh2::Event & event) override;

 private:
  TString path, SysDirection;
  std::unique_ptr<TGraphAsymmErrors> Eff_lowpt_MC, Eff_lowpt_DATA, Eff_highpt_MC, Eff_highpt_DATA;

};

class CMSTTScaleFactor : public uhh2::AnalysisModule {

 public:

  explicit CMSTTScaleFactor(uhh2::Context &ctx, std::string signal_name, TString sys_direction);
  virtual bool process(uhh2::Event &event) override;

 private:

  bool m_do_weight;
  TString m_sys_direction;
};

