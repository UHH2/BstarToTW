#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Utils.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/HOTVR/include/HOTVRJetCorrector.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

class JEC2016Module: public uhh2::AnalysisModule {
  /* 
     Jet Energy Corrections Module for 2016 era
   */
 public:
  explicit JEC2016Module(uhh2::Context & context);
  virtual bool process(uhh2::Event &event) override;
  void switch_jet_lepton_cleaning(bool b = true) {b_jetlep = b;}

 private:
  std::unique_ptr<JetCorrector> jec_ak4_mc, jec_ak4_B, jec_ak4_C, jec_ak4_D, jec_ak4_E, jec_ak4_F, jec_ak4_G, jec_ak4_H;
  std::unique_ptr<JetLeptonCleaner_by_KEYmatching> jlc_ak4_mc, jlc_ak4_B, jlc_ak4_C, jlc_ak4_D, jlc_ak4_E, jlc_ak4_F, jlc_ak4_G, jlc_ak4_H;
  std::unique_ptr<HOTVRJetCorrector> jec_hotvr_mc, jec_hotvr_B, jec_hotvr_C, jec_hotvr_D, jec_hotvr_E, jec_hotvr_F, jec_hotvr_G, jec_hotvr_H;
  std::unique_ptr<AnalysisModule> jlc_hotvr;
  std::unique_ptr<JetResolutionSmearer> jet_resolution_smearer;

  bool is_mc;
  bool b_jetlep = true;

  // Run Numbers taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2016Analysis
  const int runnr_B = 275376;
  const int runnr_C = 276283;
  const int runnr_D = 276811;
  const int runnr_E = 277420;
  const int runnr_F = 278808;
  const int runnr_G = 280385;
  const int runnr_H = 284044;
};

class JEC2017Module: public uhh2::AnalysisModule {
  /* 
     Jet Energy Corrections Module for 2016 era
   */
 public:
  explicit JEC2017Module(uhh2::Context & context);
  virtual bool process(uhh2::Event &event) override;
  void switch_jet_lepton_cleaning(bool b = true) {b_jetlep = b;}

 private:
  std::unique_ptr<JetCorrector> jec_ak4_mc, jec_ak4_B, jec_ak4_C, jec_ak4_D, jec_ak4_E, jec_ak4_F;
  std::unique_ptr<JetLeptonCleaner_by_KEYmatching> jlc_ak4_mc, jlc_ak4_B, jlc_ak4_C, jlc_ak4_D, jlc_ak4_E, jlc_ak4_F;
  std::unique_ptr<HOTVRJetCorrector> jec_hotvr_mc, jec_hotvr_B, jec_hotvr_C, jec_hotvr_D, jec_hotvr_E, jec_hotvr_F;
  std::unique_ptr<AnalysisModule> jlc_hotvr;
  std::unique_ptr<JetResolutionSmearer> jet_resolution_smearer;

  bool is_mc;
  bool b_jetlep = true;

  // Run Numbers taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2017Analysis
  const int runnr_B = 299329;
  const int runnr_C = 302029;
  const int runnr_D = 303434;
  const int runnr_E = 304826;
  const int runnr_F = 306462;
};

class JEC2018Module: public uhh2::AnalysisModule {
  /* 
     Jet Energy Corrections Module for 2018 era
   */
 public:
  explicit JEC2018Module(uhh2::Context & context);
  virtual bool process(uhh2::Event &event) override;
  void switch_jet_lepton_cleaning(bool b = true) {b_jetlep = b;}

 private:
  std::unique_ptr<AnalysisModule> cl_jetmet;
  std::unique_ptr<JetCorrector> jec_ak4_mc, jec_ak4_A, jec_ak4_B, jec_ak4_C, jec_ak4_D;
  std::unique_ptr<JetLeptonCleaner_by_KEYmatching> jlc_ak4_mc, jlc_ak4_A, jlc_ak4_B, jlc_ak4_C, jlc_ak4_D;
  std::unique_ptr<AnalysisModule> jec_hotvr_mc, jec_hotvr_A, jec_hotvr_B, jec_hotvr_C, jec_hotvr_D;
  std::unique_ptr<AnalysisModule> jlc_hotvr;
  std::unique_ptr<JetResolutionSmearer> jet_resolution_smearer;

  bool is_mc;
  bool b_jetlep = true;

  // Run Numbers taken from https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2018Analysis
  const int runnr_A = 316995;
  const int runnr_B = 319312;
  const int runnr_C = 320393;
  const int runnr_D = 325273;
};

class ObjectSetup: public uhh2::AnalysisModule {
 public:
  explicit ObjectSetup(uhh2::Context & context);
  virtual bool process(uhh2::Event &event) override;

  void switch_lepton_cleaning(bool b = true) {b_lepclean = b;}


 private:
  bool is_mc;
  bool b_lepclean = true;
  Year year;

  double lep_eta_max = 2.4;
  double lepveto_pt_min  = 30.0;
  double lep_pt_min  = 50.0; 
  double muo_iso_max = 0.15;
  double jet_pt_min  = 30.0;
  double jet_eta_max = 2.4;
  double top_pt_min = 200.0;
  double top_eta_max = 2.5;

  std::unique_ptr<AnalysisModule> jec_module;
  std::unique_ptr<AnalysisModule> cl_pv, cl_muo, cl_ele, cl_jetpfid, cl_jet, cl_topjet;
  std::unique_ptr<AndSelection> metfilters_selection;

};
