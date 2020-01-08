#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/YearRunSwitchers.h"

#include <TRandom.h>

class MuonScaleFactors2016: public uhh2::AnalysisModule {
 public:
  explicit MuonScaleFactors2016(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;
  
 private:
  std::unique_ptr<AnalysisModule> m_sf_trigger, m_sf_id, m_sf_iso;
};

class MuonScaleFactors2017: public uhh2::AnalysisModule {
 public:
  explicit MuonScaleFactors2017(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;
  
 private:
  std::unique_ptr<AnalysisModule> m_sf_trigger, m_sf_id, m_sf_iso;
};

class MuonScaleFactors2018: public uhh2::AnalysisModule {
 public:
  explicit MuonScaleFactors2018(uhh2::Context &ctx, long int seed = 123456789);
  virtual bool process(uhh2::Event &event) override;
  
 private:
  TRandom *m_rng;
  long int m_seed;
  const double m_lumi_fraction = 8.95 / 59.74;
  std::unique_ptr<AnalysisModule> m_sf_trigger_before, m_sf_trigger_after, m_sf_id, m_sf_iso;
};

class ElectronScaleFactors2016: public uhh2::AnalysisModule {
 public:
  explicit ElectronScaleFactors2016(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;
  
 private:
  std::unique_ptr<AnalysisModule> m_sf_trigger, m_sf_id, m_sf_reco;
};

class ElectronScaleFactors2017: public uhh2::AnalysisModule {
 public:
  explicit ElectronScaleFactors2017(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;
  
 private:
  std::unique_ptr<AnalysisModule> m_sf_trigger, m_sf_id, m_sf_reco;
};

class ElectronScaleFactors2018: public uhh2::AnalysisModule {
 public:
  explicit ElectronScaleFactors2018(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;
  
 private:
  std::unique_ptr<AnalysisModule> m_sf_trigger, m_sf_id, m_sf_reco;
};

class LeptonScaleFactors: public uhh2::AnalysisModule {
 public:
  explicit LeptonScaleFactors(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event);

 private:
  std::unique_ptr<YearSwitcher> m_sf_lepton;
};
