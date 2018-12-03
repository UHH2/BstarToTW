#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

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

#include "UHH2/BstarToTW/include/BstarToTWHypothesis.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>

class ObjectCleaner: public uhh2::AnalysisModule {
 public:
  explicit ObjectCleaner(uhh2::Context & context);
  virtual bool process(uhh2::Event &event) override;

 private:
  bool is_mc;
  const int runnr_BCD = 276811;
  const int runnr_EFearly = 278802;
  const int runnr_FlateG = 280385;

  double lep_eta_max = 2.4;
  double lepveto_pt_min  = 30.0;
  double lep_pt_min  = 50.0; 
  double muo_iso_max = 0.15;
  double jet_pt_min  = 30.0;
  double jet_eta_max = 2.4;
  double top_pt_min = 200.0;
  double top_eta_max = 2.5;

  std::unique_ptr<JetCorrector> jec_ak4_mc, jec_ak4_BCD, jec_ak4_EFearly, jec_ak4_FlateG, jec_ak4_H;
  std::unique_ptr<JetLeptonCleaner> jlc_mc, jlc_BCD, jlc_EFearly, jlc_FlateG, jlc_H;
  std::unique_ptr<HOTVRJetCorrector> jec_hotvr_mc, jec_hotvr_BCD, jec_hotvr_EFearly, jec_hotvr_FlateG, jec_hotvr_H;
  std::unique_ptr<AnalysisModule> cl_pv, cl_muo, cl_ele, cl_jetpfid, cl_jet, cl_topjet;
  std::unique_ptr<JetResolutionSmearer> jet_resolution_smearer;
  std::unique_ptr<AndSelection> metfilters;
};

class SystematicsModule: public uhh2::AnalysisModule {
 public:
  explicit SystematicsModule();
  virtual bool process(uhh2::Event &event) override;
  void add_module(std::unique_ptr<uhh2::AnalysisModule> &module){modules.push_back(std::move(module));}
  
 private:
  bool is_mc;
  std::vector<std::unique_ptr<uhh2::AnalysisModule>> modules;
};

class ElectronTriggerWeights: public uhh2::AnalysisModule{

 public:
  explicit ElectronTriggerWeights(uhh2::Context & ctx, TString path_, TString SysDirection_);
  virtual bool process(uhh2::Event & event) override;

 private:
  TString path, SysDirection;
  std::unique_ptr<TGraphAsymmErrors> Eff_lowpt_MC, Eff_lowpt_DATA, Eff_highpt_MC, Eff_highpt_DATA;
};

class CMSTTScaleFactor : public uhh2::AnalysisModule {

 public:

  explicit CMSTTScaleFactor(uhh2::Context &ctx, std::string signal_name, TString sys_direction);
  virtual bool process(uhh2::Event &event) override;

 private:

  bool m_do_weight;
  TString m_sys_direction;
};
