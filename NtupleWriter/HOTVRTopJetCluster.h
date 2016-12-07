#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/HOTVR.hh"
#include "fastjet/contrib/HOTVRinfo.hh"

#include "vector"

class UniversalJetCluster
{
 public:
  // TODO Overload Constructor for GenParticles and PFParticles and not other mehtods
  UniversalJetCluster(const std::vector<GenParticle>);
  UniversalJetCluster(const std::vector<PFParticle>);

 private:
  std::vector<fastjet::PseudoJet> _psj;
  std::vector<TopJet> _hotvrTop;

  fastjet::PseudoJet ConvertGenToPsj(const GenParticle & genp);
  fastjet::PseudoJet ConvertPFToPsj(const PFParticle & pfp);
  Jet ConvertPsjToJet(const fastjet::PseudoJet & psj);
  TopJet ConvertPsjToTopJet(const fastjet::PseudoJet & psj, const std::vector<fastjet::PseudoJet> subpsj);

  void ClusterHOTVR();
};
