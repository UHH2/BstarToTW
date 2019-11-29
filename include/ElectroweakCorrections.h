#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include <TFile.h>
#include <TF1.h>

class ElectroweakCorrections: public uhh2::AnalysisModule {
 public:
  explicit ElectroweakCorrections(uhh2::Context &ctx);
  virtual bool process(uhh2::Event &event) override;

 private:
  uhh2::Event::Handle<float> h_ewk_weight;
  GenParticle get_genp(uhh2::Event &event);
  TF1 *m_fun;
  int m_pdg_id = 0;

};
