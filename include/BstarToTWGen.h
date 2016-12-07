#pragma once

#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include <vector>


/** \brief bstar generator truth information
 * 
 * Interpret a vector of GenParticle as bstar event, providing easier access to
 * the various components in the bstar decay chain.
 * 
 * The class can either be used directly by passing the genparticles vector;
 * another option is to use the BstarToTWGenProducer (see below), which writes the
 * result to the event container, where it can be accessed later.
 */
class BstarToTWGen {
 public:
  
  /** construct from genparticles.
   * 
   * The event should be an actual bstar event, i.e.:
   *   - there are a top and two Ws
   *   - the top must have exactl two daughters
   *   - one of the top daughters must be a W
   *   - the Ws must have exactly two daughters
   * 
   * Note that it is allowed that more particles than those belonging to the ttbar system
   * are in genparts; those are ignored.
   * 
   * In case one of these conditions is not fulfilled, the behavior
   * depends on the throw_on_failure parameter: if it is true (the default), a runtime_error
   * is thrown with an explanation what failed. If it is false, no exception is thrown, but
   * not all contents are valid; most will return a 0 vector. The one thing guaranteed is that the
   * decaychannel will be e_notfound. If using throw_on_failure = false, it is thus a good idea
   * to check the decaychannel.
   */
  explicit BstarToTWGen(const uhh2::Event & event);

  LorentzVector bstar() const;
  LorentzVector Wbstar() const;
  LorentzVector tbstar() const;
  LorentzVector Wdecay1() const;
  LorentzVector Wdecay2() const;
  LorentzVector tW() const;
  LorentzVector tb() const;
  LorentzVector tWdecay1() const;
  LorentzVector tWdecay2() const;
  LorentzVector ChargedLepton() const;
  LorentzVector Neutrino() const;

  std::vector<LorentzVector> stable() const;
  std::vector<LorentzVector> electrons() const;
  std::vector<LorentzVector> muons() const;

  bool IsSemiLeptonicDecay() const;
  bool IsTopHadronicDecay() const;
  bool IsWHadronicDecay() const;
  bool IsMuonDecay() const;
  bool IsElectronDecay() const;

 private:
  // void MatchHardProcess(const std::vector<GenParticle> & genparticles);

  // Particles from hard process
  LorentzVector m_bstar; 
  LorentzVector m_Wbstar;
  LorentzVector m_tbstar; 
  LorentzVector m_Wdecay1;
  LorentzVector m_Wdecay2; 
  LorentzVector m_tW;
  LorentzVector m_tb;
  LorentzVector m_tWdecay1;
  LorentzVector m_tWdecay2;
  LorentzVector m_cl;
  LorentzVector m_nl;

  // Stable particle collections
  std::vector<LorentzVector> m_stable;
  std::vector<LorentzVector> m_ele;
  std::vector<LorentzVector> m_muo;

  bool is_tophad, is_whad;
  bool is_singlemuo, is_singleele;

};


class BstarToTWGenProducer: public uhh2::AnalysisModule {
 public:
  explicit BstarToTWGenProducer(uhh2::Context & ctx, const std::string & name = "BstarToTWgen");
  virtual bool process(uhh2::Event & event) override;
    
 private:
  uhh2::Event::Handle<BstarToTWGen> h_BstarToTWgen;
  bool throw_on_failure;
};


