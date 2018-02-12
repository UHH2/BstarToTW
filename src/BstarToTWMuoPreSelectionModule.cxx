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
#include "UHH2/BstarToTW/include/HOTVRHists.h"
#include "UHH2/BstarToTW/include/HOTVRJetCorrectionHists.h"
#include "UHH2/BstarToTW/include/HOTVRJetCorrector.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWMuoPreSelectionModule: public AnalysisModule {
  public:
    
    explicit BstarToTWMuoPreSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  

    // Common Modules
    std::unique_ptr<CommonModules> common;
    // JEC 
    std::unique_ptr<AnalysisModule> jec_subj_mc, jec_subj_BCD, jec_subj_EFearly, jec_subj_FlateG, jec_subj_H;

    // Cleaner
    std::unique_ptr<AnalysisModule> cl_muon;
    std::unique_ptr<AnalysisModule> cl_topjet;

    // Selections
    bool is_data;
    std::unique_ptr<Selection> sel_lumi;
    std::unique_ptr<Selection> trig_IsoMu24, trig_IsoTkMu24;
    std::unique_ptr<Selection> sel_nmuo;
    std::unique_ptr<Selection> sel_nele;
    std::unique_ptr<Selection> sel_met;
    std::unique_ptr<Selection> sel_st;
    std::unique_ptr<Selection> sel_ntop;

    // Hists
    std::unique_ptr<Hists> hist_nocuts, hist_trigger, hist_cleaner, hist_nmuo, hist_met, hist_st, hist_subj_jec, hist_jec_corr, hist_topcleaner, hist_ntop;
    std::unique_ptr<Hists> hist_pileup;

    bool is_mc;

    const int runnr_BCD = 276811;
    const int runnr_EFearly = 278802;
    const int runnr_FlateG = 280385;

  };

  BstarToTWMuoPreSelectionModule::BstarToTWMuoPreSelectionModule(Context & ctx) {

    is_mc = ctx.get("dataset_type") == "MC";

    // Kinematic Variables
    double met_min = 50.;
    double st_min = 400.;

    double lep_eta_max = 2.4;
    double lep_pt_min  = 30.0;
    double muo_pt_min  = 50.0; 
    double muo_iso_max = 0.15;

    double jet_pt_min  = 30.0;
    double jet_eta_max = 2.4;

    double top_pt_min = 200.0;
    double top_eta_max = 2.5;


    // IDs
    MuonId id_muo_loose = AndId<Muon>(MuonIDLoose(), PtEtaCut(muo_pt_min, lep_eta_max), MuonIso(muo_iso_max));
    MuonId id_muo_tight = AndId<Muon>(MuonIDTight(), PtEtaCut(muo_pt_min, lep_eta_max), MuonIso(muo_iso_max));
    ElectronId id_ele = AndId<Electron>(ElectronID_Spring16_veto_noIso, PtEtaCut(lep_pt_min, lep_eta_max));
    JetId id_jet = AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(jet_pt_min, jet_eta_max));
    TopJetId id_topjet =  PtEtaCut(top_pt_min, top_eta_max);

    // Additional Modules
    common.reset(new CommonModules());
    common->switch_jetlepcleaner(true);
    common->set_muon_id(id_muo_loose);
    common->set_electron_id(id_ele);
    common->set_jet_id(id_jet);
    common->init(ctx);

    if(is_mc)
      {
    	jec_subj_mc.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L23_AK4PFchs_MC));
      }
    else
      { 
    	jec_subj_BCD.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_BCD_L23_AK4PFchs_DATA));
    	jec_subj_EFearly.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_EF_L23_AK4PFchs_DATA));
    	jec_subj_FlateG.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_G_L23_AK4PFchs_DATA));
    	jec_subj_H.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_H_L23_AK4PFchs_DATA));

      }
    hist_jec_corr.reset(new HOTVRJetCorrectionHists(ctx, "JEC_corrections"));

    // Cleaner
    cl_muon.reset(new MuonCleaner(id_muo_tight));
    cl_topjet.reset(new TopJetCleaner(ctx, id_topjet));

    // Selections
    sel_lumi.reset(new LumiSelection(ctx));
    trig_IsoMu24.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    trig_IsoTkMu24.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    sel_nmuo.reset(new NMuonSelection(1, 1));
    sel_nele.reset(new NElectronSelection(0, 0));
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
    hist_subj_jec.reset(new AndHists(ctx, "Subjet_Corrections"));
    hist_topcleaner.reset(new AndHists(ctx, "Topjet_Cleaning"));
    hist_ntop.reset(new AndHists(ctx, "1TopJetCut"));

    hist_pileup.reset(new HOTVRPileUpHists(ctx, "HOTVR_PileUp"));
  }

  bool BstarToTWMuoPreSelectionModule::process(Event & event) {

    if(!is_mc)
      {
	if(!sel_lumi->passes(event)) return false;
      }

    hist_nocuts->fill(event);

    // Trigger
    if(!(trig_IsoMu24->passes(event) || trig_IsoTkMu24->passes(event))) return false;
    hist_trigger->fill(event);

    // Cleaner
    if(!common->process(event)) return false;
    hist_cleaner->fill(event);

    // Muon selection & additional lepton veto
    if(!sel_nmuo->passes(event) || !sel_nele->passes(event)) return false;
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
    hist_pileup->fill(event);
    if (is_mc) jec_subj_mc->process(event);
    else
      {
	if(event.run <= runnr_BCD)         jec_subj_BCD->process(event);
	else if(event.run < runnr_EFearly) jec_subj_EFearly->process(event); //< is correct, not <= 
	else if(event.run <= runnr_FlateG) jec_subj_FlateG->process(event);
	else if(event.run > runnr_FlateG)  jec_subj_H->process(event);
      }
    // if (is_mc)
    //   {
    // 	hist_jec_corr->fill(event);
    //   }

    cl_topjet->process(event);
    hist_topcleaner->fill(event);

    // TopJet Selection
    if(!sel_ntop->passes(event)) return false;
    hist_ntop->fill(event);

    // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWMuoPreSelectionModule)

}
