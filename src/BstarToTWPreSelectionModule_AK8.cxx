#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/BstarToTWHists.h"
#include "UHH2/BstarToTW/include/BstarToTWSelections.h"

#include "UHH2/HOTVR/include/HOTVRHists.h"
#include "UHH2/HOTVR/include/HOTVRJetCorrectionHists.h"
#include "UHH2/HOTVR/include/HOTVRJetCorrector.h"

using namespace std;
using namespace uhh2;
namespace JERFiles {
  extern const std::vector<std::string> Fall17_17Nov2017_V6_L123_AK8PFchs_MC;

  extern const std::vector<std::string> Fall17_17Nov2017_V6_B_L123_AK8PFchs_DATA;
  
  extern const std::vector<std::string> Fall17_17Nov2017_V6_C_L123_AK8PFchs_DATA;
  
  extern const std::vector<std::string> Fall17_17Nov2017_V6_D_L123_AK8PFchs_DATA;
  
  extern const std::vector<std::string> Fall17_17Nov2017_V6_E_L123_AK8PFchs_DATA;
  
  extern const std::vector<std::string> Fall17_17Nov2017_V6_F_L123_AK8PFchs_DATA;
}

namespace uhh2 {

  class BstarToTWPreSelectionModule_AK8: public AnalysisModule {
  public:
    
    explicit BstarToTWPreSelectionModule_AK8(Context & ctx);
    virtual bool process(Event & event) override;

  private:  

    // Common Modules
    std::unique_ptr<CommonModules> common;
    std::unique_ptr<AnalysisModule> jec_topj_mc, jec_topj_B, jec_topj_C, jec_topj_D, jec_topj_E, jec_topj_F;

    // Cleaner
    std::unique_ptr<AnalysisModule> cl_muon;
    std::unique_ptr<AnalysisModule> cl_topjet;

    // Selections
    bool is_data;
    std::unique_ptr<Selection> sel_lumi;
    std::unique_ptr<Selection> trig_IsoMu27;//, trig_IsoTkMu24;
    std::unique_ptr<Selection> sel_nmuo;
    std::unique_ptr<Selection> sel_met;
    std::unique_ptr<Selection> sel_st;
    std::unique_ptr<Selection> sel_ntop;

    // Hists
    std::unique_ptr<Hists> hist_nocuts, hist_trigger, hist_cleaner, hist_nmuo, hist_met, hist_st, hist_topcleaner, hist_ntop;

    bool is_mc;

    const int runnr_B = 299329; 
    const int runnr_C = 302029; 
    const int runnr_D = 303434; 
    const int runnr_E = 304797; 

  };

  BstarToTWPreSelectionModule_AK8::BstarToTWPreSelectionModule_AK8(Context & ctx) {

    is_mc = ctx.get("dataset_type") == "MC";

    // Kinematic Variables
    double met_min = 50.;
    double st_min = 400.;

    double lep_eta_max = 2.4;
    double muo_pt_min  = 50.; 
    double muo_iso_max = 0.15;

    double ele_pt_min  = 30.0;

    double jet_pt_min  = 30.0;
    double jet_eta_max = 2.4;

    double top_pt_min = 200.0;
    double top_eta_max = 2.5;


    // IDs
    // MuonId id_muo_loose = AndId<Muon>(MuonIDLoose(), PtEtaCut(muo_pt_min, lep_eta_max), MuonIso(muo_iso_max));
    MuonId id_muo_loose = AndId<Muon>(MuonID(Muon::Selector::CutBasedIdLoose),MuonID(Muon::Selector::PFIsoTight));
    // MuonId id_muo_tight = AndId<Muon>(MuonIDTight(), PtEtaCut(muo_pt_min, lep_eta_max), MuonIso(muo_iso_max));
    MuonId id_muo_tight = AndId<Muon>(MuonID(Muon::Selector::CutBasedIdTight), MuonID(Muon::Selector::PFIsoTight));
    ElectronId id_ele = AndId<Electron>(ElectronID_Spring16_veto_noIso, PtEtaCut(ele_pt_min, lep_eta_max));
    JetId id_jet = AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_LEPVETO), PtEtaCut(jet_pt_min, jet_eta_max));
    TopJetId id_topjet =  PtEtaCut(top_pt_min, top_eta_max);

    // Additional Modules
    common.reset(new CommonModules());
    common->switch_jetlepcleaner(true);
    common->set_muon_id(id_muo_loose);
    common->set_electron_id(id_ele);
    common->set_jet_id(id_jet);
    if(is_mc) common->disable_metfilters();
    common->init(ctx);

    if(is_mc)
      {	
	jec_topj_mc.reset(new TopJetCorrector(ctx, JERFiles::Fall17_17Nov2017_V6_L123_AK8PFchs_MC));
      }
    else
      { 
 	jec_topj_B.reset(new TopJetCorrector(ctx, JERFiles::Fall17_17Nov2017_V6_B_L123_AK8PFchs_DATA));
    	jec_topj_C.reset(new TopJetCorrector(ctx, JERFiles::Fall17_17Nov2017_V6_C_L123_AK8PFchs_DATA));
    	jec_topj_D.reset(new TopJetCorrector(ctx, JERFiles::Fall17_17Nov2017_V6_D_L123_AK8PFchs_DATA));
    	jec_topj_E.reset(new TopJetCorrector(ctx, JERFiles::Fall17_17Nov2017_V6_E_L123_AK8PFchs_DATA));
    	jec_topj_F.reset(new TopJetCorrector(ctx, JERFiles::Fall17_17Nov2017_V6_F_L123_AK8PFchs_DATA));
      }

    // Cleaner
    cl_muon.reset(new MuonCleaner(id_muo_tight));
    cl_topjet.reset(new TopJetCleaner(ctx, id_topjet));

    // Selections
    sel_lumi.reset(new LumiSelection(ctx));
    trig_IsoMu27.reset(new TriggerSelection("HLT_IsoMu27_v*"));
    // trig_IsoTkMu24.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    sel_nmuo.reset(new NMuonSelection(1, 1));
    sel_met.reset(new METSelection(met_min));
    sel_st.reset(new STSelection(st_min));
    sel_ntop.reset(new NTopJetSelection(1, -1));

    // Hists
    hist_nocuts.reset(new AndHists(ctx, "NoCuts"));
    hist_trigger.reset(new AndHists(ctx, "Trigger"));
    hist_cleaner.reset(new AndHists(ctx, "Cleaning"));
    hist_nmuo.reset(new AndHists(ctx, "1MuonCut"));
    hist_met.reset(new AndHists(ctx, "50METCut"));
    hist_st.reset(new AndHists(ctx, "350STCut"));
    hist_topcleaner.reset(new AndHists(ctx, "Topjet_Cleaning"));
    hist_ntop.reset(new AndHists(ctx, "1TopJetCut"));

  }

  bool BstarToTWPreSelectionModule_AK8::process(Event & event) {

    if(!is_mc)
      {
	if(!sel_lumi->passes(event)) return false;
      }

    hist_nocuts->fill(event);

    // Trigger
    if(!trig_IsoMu27->passes(event)) return false;
    hist_trigger->fill(event);

    // Cleaner
    if(!common->process(event)) return false;
    hist_cleaner->fill(event);

    // Muon selection
    if(!sel_nmuo->passes(event)) return false;
    cl_muon->process(event);
    if(!sel_nmuo->passes(event)) return false;
    hist_nmuo->fill(event);

    // MET selection
    if(!sel_met->passes(event)) return false;
    hist_met->fill(event);

    // ST selection
    if(!sel_st->passes(event)) return false;
    hist_st->fill(event);

    // TopJEC
    if (is_mc) jec_topj_mc->process(event);
    else
      {
	if(event.run <= runnr_B)      jec_topj_B->process(event);
	else if(event.run <= runnr_C) jec_topj_C->process(event);
	else if(event.run <= runnr_D) jec_topj_D->process(event);
	else if(event.run <= runnr_E) jec_topj_E->process(event);
	else if(event.run >  runnr_E) jec_topj_F->process(event);
      }
    cl_topjet->process(event);
    hist_topcleaner->fill(event);

    // TopJet Selection
    if(!sel_ntop->passes(event)) return false;
    hist_ntop->fill(event);

    // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWPreSelectionModule_AK8)

}

const std::vector<std::string> JERFiles::Fall17_17Nov2017_V6_L123_AK8PFchs_MC = {
  "JECDatabase/textFiles/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L1FastJet_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_L3Absolute_AK8PFchs.txt",
};
const std::vector<std::string> JERFiles::Fall17_17Nov2017_V6_B_L123_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L1FastJet_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017B_V6_DATA/Fall17_17Nov2017B_V6_DATA_L2L3Residual_AK8PFchs.txt",
};
const std::vector<std::string> JERFiles::Fall17_17Nov2017_V6_C_L123_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L1FastJet_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017C_V6_DATA/Fall17_17Nov2017C_V6_DATA_L2L3Residual_AK8PFchs.txt",
};
const std::vector<std::string> JERFiles::Fall17_17Nov2017_V6_D_L123_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017D_V6_DATA/Fall17_17Nov2017D_V6_DATA_L1FastJet_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V6_DATA/Fall17_17Nov2017D_V6_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V6_DATA/Fall17_17Nov2017D_V6_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017D_V6_DATA/Fall17_17Nov2017D_V6_DATA_L2L3Residual_AK8PFchs.txt",
};
const std::vector<std::string> JERFiles::Fall17_17Nov2017_V6_E_L123_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L1FastJet_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017E_V6_DATA/Fall17_17Nov2017E_V6_DATA_L2L3Residual_AK8PFchs.txt",
};
const std::vector<std::string> JERFiles::Fall17_17Nov2017_V6_F_L123_AK8PFchs_DATA = {
  "JECDatabase/textFiles/Fall17_17Nov2017F_V6_DATA/Fall17_17Nov2017F_V6_DATA_L1FastJet_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V6_DATA/Fall17_17Nov2017F_V6_DATA_L2Relative_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V6_DATA/Fall17_17Nov2017F_V6_DATA_L3Absolute_AK8PFchs.txt",
  "JECDatabase/textFiles/Fall17_17Nov2017F_V6_DATA/Fall17_17Nov2017F_V6_DATA_L2L3Residual_AK8PFchs.txt",
};
