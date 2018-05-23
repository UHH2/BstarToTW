#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/BstarToTW/include/BstarToTWHypothesis.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>


class ElectronTriggerWeights: public uhh2::AnalysisModule{

 public:
  explicit ElectronTriggerWeights(uhh2::Context & ctx, TString path_, TString SysDirection_);
  virtual bool process(uhh2::Event & event) override;

 private:
  TString path, SysDirection;
  std::unique_ptr<TGraphAsymmErrors> Eff_lowpt_MC, Eff_lowpt_DATA, Eff_highpt_MC, Eff_highpt_DATA;

};

class BackgroundShapeNormWeights: public uhh2::AnalysisModule{

 public:
  explicit BackgroundShapeNormWeights(uhh2::Context &ctx, const TString &path, const std::string &hyps_name, const std::string &discriminator_name, const std::string &sys_direction="nominal");
  virtual bool process(uhh2::Event & event) override;

 private:
  double m_p0, m_p1;
  std::string m_sys_direction, m_hyps_name, m_discriminator_name;
  uhh2::Event::Handle<std::vector<BstarToTWHypothesis>> h_hyps;
};

class CMSTTScaleFactor : public uhh2::AnalysisModule {

 public:

  explicit CMSTTScaleFactor(uhh2::Context &ctx, std::string signal_name, TString sys_direction);
  virtual bool process(uhh2::Event &event) override;

 private:

  bool m_do_weight;
  TString m_sys_direction;
};

