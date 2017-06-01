#pragma once

#include "UHH2/core/include/LorentzVector.h"
#include <map>

class BstarToTWHypothesis {
public:
  explicit BstarToTWHypothesis(){};

  LorentzVector get_topjet() const{return m_topjet;}
  LorentzVector get_neutrino() const{return m_neutrino;} 
  LorentzVector get_lepton() const{return m_lepton;}
  LorentzVector get_w() const{return (m_lepton + m_neutrino);}

  /// get the discriminator value for this hypothesis; thows a runtime_error if it does not exist.
  float get_discriminator(const std::string & label) const {
      auto it = m_discriminators.find(label);
      if(it == m_discriminators.end()){
          throw std::runtime_error("ReconstructionHypothesis::discriminator: discriminator with label '" + label + "' not set");
      }
      return it->second;
  }
  
  /// test if a discriminator value with a certian label has already been added
  bool has_discriminator(const std::string & label) const {
      return m_discriminators.find(label) != m_discriminators.end();
  }

  void set_neutrino(LorentzVector v4){m_neutrino = v4;}
  void set_lepton(const LorentzVector v4){m_lepton = v4;}
  void set_topjet(const LorentzVector v4){m_topjet = v4;}
  void set_discriminator(const std::string & label, float discr){m_discriminators[label] = discr;}
  
private:

  LorentzVector m_neutrino;
  LorentzVector m_lepton;
  LorentzVector m_topjet;

  std::map<std::string, float> m_discriminators;
};
