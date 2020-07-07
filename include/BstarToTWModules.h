#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TTbarGen.h" 
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
  JetId m_btag_tight_id = DeepCSVBTag(DeepCSVBTag::WP_TIGHT);
  JetId m_btag_medium_id = DeepCSVBTag(DeepCSVBTag::WP_MEDIUM);
  JetId m_btag_loose_id = DeepCSVBTag(DeepCSVBTag::WP_LOOSE);

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

class ElectronTriggerWeights: public uhh2::AnalysisModule{

 public:
  explicit ElectronTriggerWeights(uhh2::Context & ctx, TString path_, TString syst_direction_);
  virtual bool process(uhh2::Event & event) override;

 private:
  TString path, syst_direction;
  Year year;
  std::vector<float> pt_bins;
  std::unique_ptr<TGraphAsymmErrors> g_sf_pt1, g_sf_pt2;
  uhh2::Event::Handle<float> h_ele_weight, h_ele_weight_up, h_ele_weight_down;

};

class MuonTriggerWeights: public uhh2::AnalysisModule{

 public:
  explicit MuonTriggerWeights(uhh2::Context & ctx, TString path_, TString syst_direction_);
  virtual bool process(uhh2::Event & event) override;

 private:
  TString path, syst_direction;
  Year year;
  std::unique_ptr<TGraphAsymmErrors> g_sf_0to0p9, g_sf_0p9to1p2, g_sf_1p2to2p1, g_sf_2p1to2p4;
  uhh2::Event::Handle<float> h_muo_weight, h_muo_weight_up, h_muo_weight_down;

};

class TopPtReweighting : public uhh2::AnalysisModule {
 public:
  explicit TopPtReweighting(uhh2::Context& ctx,
			 float a, float b,
			 const std::string& syst_a,
			 const std::string& syst_b,
			 const std::string& ttgen_name ="",
			 const std::string& weight_name="weight_ttbar");


  virtual bool process(uhh2::Event& event) override;
 private:
  uhh2::Event::Handle<TTbarGen> h_ttbargen_;
  uhh2::Event::Handle<float> h_weight_;
  float a_, b_;
  std::string version_;
  std::string ttgen_name_;
  std::string weight_name_;
};
