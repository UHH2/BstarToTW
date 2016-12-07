#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/HOTVR.hh"
#include "fastjet/contrib/HOTVRinfo.hh"

#include "vector"

class GenJetCluster
{
 public:
  GenJetCluster(const uhh2::Event & event);
  std::vector<fastjet::PseudoJet> getPseudojets() const;  
  std::vector<fastjet::PseudoJet> getAK4Jets() const;
  std::vector<fastjet::PseudoJet> getAK8Jets() const;  
  std::vector<fastjet::PseudoJet> getCMSTopTagged() const;
  std::vector<fastjet::PseudoJet> getHOTVRTopTagged() const;

 private:
  fastjet::PseudoJet convertGenParticle(const GenParticle & genp);
  std::vector<fastjet::PseudoJet> _stablePsj;
  std::vector<fastjet::PseudoJet> _ak4jets;
  std::vector<fastjet::PseudoJet> _ak8jets;
  std::vector<fastjet::PseudoJet> _cmstoptagged;
  std::vector<fastjet::PseudoJet> _hotvrtoptagged;
};

class GenJetProducer: public uhh2::AnalysisModule
{
 public:
  explicit GenJetProducer(uhh2::Context & ctx, const std::string & name);
  virtual bool process(uhh2::Event & event) override;

 private:
  uhh2::Event::Handle<GenJetCluster> h_cluster;
  uhh2::Event::Handle<GenJetCluster> h_genjetcluster;
};
