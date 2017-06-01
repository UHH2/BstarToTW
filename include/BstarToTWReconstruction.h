#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/PrimaryLepton.h"

#include "UHH2/BstarToTW/include/BstarToTWHypothesis.h"
#include "UHH2/BstarToTW/include/HOTVRIds.h"


typedef std::function< std::vector<LorentzVector> (const LorentzVector & lepton, const LorentzVector & met)> NeutrinoReconstructionMethod;


class BstarToTWReconstruction: public uhh2::AnalysisModule {
 public:
  explicit BstarToTWReconstruction(uhh2::Context&, const NeutrinoReconstructionMethod&, const std::string& label="BstarToTWReco", TopJetId id=HOTVRTopTag());
  virtual bool process(uhh2::Event&) override;

 private:
  NeutrinoReconstructionMethod m_neutrinofunction;
  uhh2::Event::Handle<std::vector<BstarToTWHypothesis>> h_recohyps;
  uhh2::Event::Handle<FlavorParticle> h_primlep;

  TopJetId m_topjetid;
};
