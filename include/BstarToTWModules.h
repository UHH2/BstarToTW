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
#include "UHH2/HOTVR/include/HOTVRIds.h"

#include "UHH2/BstarToTW/include/BstarToTWHypothesis.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"

#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>

class ObjectCleaner: public uhh2::AnalysisModule {
 public:
  explicit ObjectCleaner(uhh2::Context & context);
  virtual bool process(uhh2::Event &event) override;

  void switch_lepton_cleaning(bool b = true) {b_lepclean = b;}
  void switch_jet_lepton_cleaning(bool b = true) {b_jetlep = b;}

 private:
  bool is_mc;
  bool b_lepclean = true;
  bool b_jetlep = true;
  
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
  std::unique_ptr<AnalysisModule> cl_pv, cl_muo, cl_ele, cl_jetpfid, cl_jet, cl_topjet, jlc_hotvr;
  std::unique_ptr<JetResolutionSmearer> jet_resolution_smearer;
  std::unique_ptr<AndSelection> metfilters;
};

class ObjectTagger: public uhh2::AnalysisModule {
 public:
  explicit ObjectTagger();
  virtual bool process(uhh2::Event &event) override;
  
  void set_toptag_name(std::string name) {fail_if_init(); m_toptag_name = name;}
  void set_toptag_id(TopJetId id) {fail_if_init(); m_toptag_id = id;}
  std::string get_toptag_name() {return m_toptag_name;}

  void set_btag_loose_name(std::string name) {fail_if_init(); m_btag_loose_name = name;}
  void set_btag_loose_id(JetId id) {fail_if_init(); m_btag_loose_id = id;}
  std::string get_btag_loose_name() {return m_btag_loose_name;}

  void set_btag_medium_name(std::string name) {fail_if_init(); m_btag_medium_name = name;}
  void set_btag_medium_id(JetId id) {fail_if_init(); m_btag_medium_id = id;}
  std::string get_btag_medium_name() {return m_btag_medium_name;}

  void set_btag_tight_name(std::string name) {fail_if_init(); m_btag_tight_name = name;}
  void set_btag_tight_id(JetId id) {fail_if_init(); m_btag_tight_id = id;}
  std::string get_btag_tight_name() {return m_btag_tight_name;}  
  
  void init(Context &ctx);

 private:
  bool b_init = false;
  void fail_if_init();

  std::string m_toptag_name = "toptag";
  std::string m_btag_loose_name = "btag_loose";
  std::string m_btag_medium_name = "btag_medium";
  std::string m_btag_tight_name = "btag_tight";

  uhh2::Event::Handle<std::vector<TopJet> > h_toptag;
  uhh2::Event::Handle<std::vector<Jet> > h_btag_loose;
  uhh2::Event::Handle<std::vector<Jet> > h_btag_medium;
  uhh2::Event::Handle<std::vector<Jet> > h_btag_tight;

  TopJetId m_toptag_id = AndId<TopJet>(HOTVRTopTag(0.8, 140.0, 220.0, 50.0), Tau32Groomed(0.56));
  JetId m_btag_tight_id = CSVBTag(CSVBTag::WP_TIGHT);
  JetId m_btag_medium_id = CSVBTag(CSVBTag::WP_MEDIUM);
  JetId m_btag_loose_id = CSVBTag(CSVBTag::WP_LOOSE);

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

class BstarToTWOutputModule:public uhh2::AnalysisModule {
 public:
  explicit BstarToTWOutputModule(uhh2::Context &ctx, const std::string hyps_name);
  virtual bool process(uhh2::Event &event) override;

 private:
  uhh2::Event::Handle<std::vector<std::vector<double> > > h_topjets;
  uhh2::Event::Handle<float> h_weight;
  uhh2::Event::Handle<float> h_recomass;
  uhh2::Event::Handle<std::vector<BstarToTWHypothesis> > h_hyps;
};
