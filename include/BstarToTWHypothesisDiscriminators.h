#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/BstarToTW/include/BstarToTWHypothesis.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"



const BstarToTWHypothesis * get_best_hypothesis(const std::vector<BstarToTWHypothesis> & hyps, const std::string & label);
const LeptonicTopHypothesis * get_best_hypothesis(const std::vector<LeptonicTopHypothesis> & hyps, const std::string & label);


class BstarToTWChi2Discriminator: public uhh2::AnalysisModule {

 public:

  BstarToTWChi2Discriminator(uhh2::Context&, const std::string&);
  virtual bool process(uhh2::Event&) override;

  virtual void set_mtop_mean (const float m){ m_mtop_mean  = m; }
  virtual void set_mtop_sigma(const float s){
    m_mtop_sigma = s; 
    if(s <= 0.) throw std::runtime_error("Chi2Discriminator::set_mtop_sigma -- logic error: non-positive input value: "+std::to_string(s));
  }

  virtual float get_mtop_mean () const { return m_mtop_mean; }
  virtual float get_mtop_sigma() const { return m_mtop_sigma; }

  virtual void set_deltaPhi_mean (const float m){ m_deltaPhi_mean  = m; }
  virtual void set_deltaPhi_sigma(const float s){
    m_deltaPhi_sigma = s; 
    if(s <= 0.) throw std::runtime_error("Chi2Discriminator::set_deltaPhi_sigma -- logic error: non-positive input value: "+std::to_string(s));
  }

  virtual float get_deltaPhi_mean () const { return m_deltaPhi_mean; }
  virtual float get_deltaPhi_sigma() const { return m_deltaPhi_sigma; }

  virtual void set_deltaPt_mean (const float m){ m_deltaPt_mean  = m; }
  virtual void set_deltaPt_sigma(const float s){
    m_deltaPt_sigma = s; 
    if(s <= 0.) throw std::runtime_error("Chi2Discriminator::set_deltaPt_sigma -- logic error: non-positive input value: "+std::to_string(s));
  }

  virtual float get_deltaPt_mean () const { return m_deltaPt_mean; }
  virtual float get_deltaPt_sigma() const { return m_deltaPt_sigma; }

 private:
  uhh2::Event::Handle<std::vector<BstarToTWHypothesis>> h_hyps;
  uhh2::Event::Handle<std::vector<Jet>> h_bjets;
  float m_mtop_mean, m_mtop_sigma;
  float m_deltaPhi_mean, m_deltaPhi_sigma;
  float m_deltaPt_mean, m_deltaPt_sigma;

};

class LeptonicTopChi2Discriminator: public uhh2::AnalysisModule {

 public:

  LeptonicTopChi2Discriminator(uhh2::Context&, const std::string&);
  virtual bool process(uhh2::Event&) override;

  virtual void set_mtop_mean (const float m){ m_mtop_mean  = m; }
  virtual void set_mtop_sigma(const float s){
    m_mtop_sigma = s; 
    if(s <= 0.) throw std::runtime_error("Chi2Discriminator::set_mtop_sigma -- logic error: non-positive input value: "+std::to_string(s));
  }

  virtual float get_mtop_mean () const { return m_mtop_mean; }
  virtual float get_mtop_sigma() const { return m_mtop_sigma; }

 private:
  uhh2::Event::Handle<std::vector<LeptonicTopHypothesis>> h_hyps;
  float m_mtop_mean, m_mtop_sigma;

};


/** \brief Try to match the reconstruction hypotheses to Monte-Carlo truth, jet-by-jet
 * 
 * Requires a TTbarGen object in the event (see TTbarGen.h).
 * 
 * Writes a "BstarToTWMatch" quality flags to the reconstruction hypotheses, which is the sum of Delta R values
 * between the four generated and reconstructed matrix-element final-state partons and the DeltaR between the
 * true neutrino and the reconstructed neutrino. The discriminator is set to infinity
 * if one of the final-state partons could not be matched to a jet within Delta R < 0.3 (note that no
 * such match is done for the neutrino).
 * 
 * NOTE: This class only works for events which are (on gen-level) electron+jets or muon+jets. Otherwise, all discriminator
 * values are set to +infinity. The reconstructed lepton is ignored in this discriminator criterion.
 */
class BstarToTWMatchDiscriminator: public uhh2::AnalysisModule {
public:
    struct cfg {
        std::string bstartotwgen_name;
        std::string discriminator_label;
        cfg(): bstartotwgen_name("BstarToTWgen"), discriminator_label("Match"){}
    };
    
    BstarToTWMatchDiscriminator(uhh2::Context & ctx, const std::string & rechyps_name, const cfg & config = cfg());
    virtual bool process(uhh2::Event & event) override;
    
private:
    
    uhh2::Event::Handle<std::vector<BstarToTWHypothesis>> h_hyps;
    uhh2::Event::Handle<BstarToTWGen> h_bstartotwgen;
    cfg config;
};
