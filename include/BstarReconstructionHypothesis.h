#pragma once

#include "UHH2/core/include/LorentzVector.h"

class BstarReconstructionHypothesis {
 public:
  LorentzVector get_Bstar() const{return m_bstar;}
  LorentzVector get_Top() const{return m_top;}
  LorentzVector get_Muon() const{return m_muon;}
  LorentzVector get_Neutrino() const{return m_neutrino;}
  LorentzVector get_W() const{return m_muon + m_neutrino;}

  void set_Bstar(LorentzVector v4) {m_bstar = v4;}
  void set_Top(LorentzVector v4) {m_top = v4;}
  void set_Muon(LorentzVector v4) {m_muon = v4;}
  void set_Neutrino(LorentzVector v4) {m_neutrino = v4;}

 private:
  LorentzVector m_bstar;
  LorentzVector m_top;
  LorentzVector m_muon;
  LorentzVector m_neutrino;

};
