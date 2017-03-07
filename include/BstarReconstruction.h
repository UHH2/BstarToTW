#pragma once

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/LorentzVector.h"

#include "UHH2/common/include/TTbarReconstruction.h"

#include "UHH2/BstarToTW/include/BstarReconstructionHypothesis.h"

#include <vector>



class BstarReconstruction: public uhh2::AnalysisModule {
typedef std::function< std::vector<LorentzVector> (const LorentzVector & lepton, const LorentzVector & met)> NeutrinoReconstructionMethod;
 public:
  explicit BstarReconstruction(uhh2::Context&, const NeutrinoReconstructionMethod&,const std::string &name = "BstarToTWReconstruction");
  virtual bool process(uhh2::Event &event) override;

 private:
  NeutrinoReconstructionMethod m_NeutrinoReco;

  uhh2::Event::Handle<std::vector<BstarReconstructionHypothesis>> h_BstarHyps;
  
};
