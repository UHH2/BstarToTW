#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/ObjectIdUtils.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/HOTVRGenIds.h"

class GenTopJetCleaner: public uhh2::AnalysisModule {
 public:
  explicit GenTopJetCleaner(uhh2::Context & ctx, GenTopJetId gentopjetid);
    virtual bool process(uhh2::Event & event) override;

 private:
    GenTopJetId m_gentopjetid;
};

class GenJetCleaner: public uhh2::AnalysisModule {
 public:
  explicit GenJetCleaner(uhh2::Context & ctx, double pt, double eta);
    virtual bool process(uhh2::Event & event) override;

 private:
    double m_pt;
    double m_eta;
};
