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

namespace uhh2 {

  class BstarToTWPreSelectionModule: public AnalysisModule {
  public:
    
    explicit BstarToTWPreSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  

    // Common Modules
    std::unique_ptr<CommonModules> common;
    // JEC 
    std::unique_ptr<AnalysisModule> jec_subj_mc, jec_subj_BCD, jec_subj_EFearly, jec_subj_FlateG, jec_subj_H;

    // Cleaner
    std::unique_ptr<AnalysisModule> cl_muo, cl_ele;
    std::unique_ptr<AnalysisModule> cl_topjet;

    // Selections
    bool is_data;
    std::unique_ptr<Selection> sel_lumi;
    std::unique_ptr<Selection> trig_IsoMu24, trig_IsoTkMu24;
    std::unique_ptr<Selection> trig_Ele27, trig_Ele115;
    std::unique_ptr<Selection> veto_muo;
    std::unique_ptr<Selection> veto_ele;
    std::unique_ptr<Selection> sel_1muo;
    std::unique_ptr<Selection> sel_1ele;
    std::unique_ptr<Selection> sel_met;
    std::unique_ptr<Selection> sel_st;
    std::unique_ptr<Selection> sel_ntop;

    // Hists
    std::unique_ptr<Hists> hist_nocuts, hist_trigger, hist_cleaner, hist_1lep, hist_met, hist_st, hist_subj_jec, hist_topcleaner, hist_ntop;
    std::unique_ptr<Hists> hist_pileup;

    bool is_mc, is_ele, is_muo;

    const int runnr_BCD = 276811;
    const int runnr_EFearly = 278802;
    const int runnr_FlateG = 280385;

  };

  BstarToTWPreSelectionModule::BstarToTWPreSelectionModule(Context & ctx) {

    // --- Config ---
    is_mc = ctx.get("dataset_type") == "MC";
    is_ele = ctx.get("analysis_channel") == "ELECTRON";
    is_muo = ctx.get("analysis_channel") == "MUON";

    // -- Kinematic Variables --
    double met_min = 50.;
    double st_min = 400.;

    double lep_eta_max = 2.4;
    double lepveto_pt_min  = 30.0;
    double lep_pt_min  = 50.0; 
    double muo_iso_max = 0.15;

    double jet_pt_min  = 30.0;
    double jet_eta_max = 2.4;

    double top_pt_min = 200.0;
    double top_eta_max = 2.5;


    // --- IDs ---
    MuonId id_muo_veto = AndId<Muon>(MuonIDLoose(), PtEtaCut(lepveto_pt_min, lep_eta_max), MuonIso(muo_iso_max)); // muon veto ID
    MuonId id_muo_tight = AndId<Muon>(MuonIDTight(), PtEtaCut(lep_pt_min, lep_eta_max), MuonIso(muo_iso_max)); // muon ID

    ElectronId id_ele_veto = AndId<Electron>(ElectronID_Spring16_veto, PtEtaCut(lepveto_pt_min, lep_eta_max)); // electron veto ID
    ElectronId id_ele_tight = AndId<Electron>(ElectronID_Spring16_tight, PtEtaCut(lep_pt_min, lep_eta_max)); // electron ID

    JetId id_jet = AndId<Jet>(JetPFID(JetPFID::WP_TIGHT), PtEtaCut(jet_pt_min, jet_eta_max)); // jet ID
    TopJetId id_topjet =  PtEtaCut(top_pt_min, top_eta_max); // maybe also deltaPhiCut ??

    // --- Additional Modules ---
    common.reset(new CommonModules());
    common->switch_jetlepcleaner(true);
    common->set_muon_id(id_muo_veto);
    common->set_electron_id(id_ele_veto);
    common->set_jet_id(id_jet);
    common->init(ctx);

    if(is_mc)
      {
    	jec_subj_mc.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK4PFPuppi_MC));
      }
    else
      { 
    	jec_subj_BCD.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFPuppi_DATA));
    	jec_subj_EFearly.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFPuppi_DATA));
    	jec_subj_FlateG.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFPuppi_DATA));
    	jec_subj_H.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFPuppi_DATA));
      }

    // Cleaner
    cl_muo.reset(new MuonCleaner(id_muo_tight));
    cl_ele.reset(new ElectronCleaner(id_ele_tight));
    cl_topjet.reset(new TopJetCleaner(ctx, id_topjet));

    // --- Selections ---
    // - Lumi -
    sel_lumi.reset(new LumiSelection(ctx));
    // - Trigger -
    // Muon
    trig_IsoMu24.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    trig_IsoTkMu24.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    // Electron
    trig_Ele27.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));
    trig_Ele115.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*")); // 105
    // - Lepton vetos
    veto_muo.reset(new NMuonSelection(0,0));
    veto_ele.reset(new NElectronSelection(0,0));
    // - Lepton
    sel_1muo.reset(new NMuonSelection(1,1));
    sel_1ele.reset(new NElectronSelection(1,1));
    // - MET
    sel_met.reset(new METSelection(met_min));
    // - ST
    sel_st.reset(new STSelection(st_min));
    // - TopJet
    sel_ntop.reset(new NTopJetSelection(1, -1));

    // --- Hists ---
    // - Initial -
    hist_nocuts.reset(new AndHists(ctx, "NoCuts"));
    // - Trigger -
    hist_trigger.reset(new AndHists(ctx, "Trigger"));
    // - Cleaning -
    hist_cleaner.reset(new AndHists(ctx, "Cleaning"));
    // - Lepton -
    hist_1lep.reset(new AndHists(ctx, "1LepCut"));
    // - MET -
    hist_met.reset(new AndHists(ctx, "50METCut"));
    // - ST -
    hist_st.reset(new AndHists(ctx, "350STCut"));
    // - HOTVR JEC -
    hist_subj_jec.reset(new AndHists(ctx, "Subjet_Corrections"));
    // - TopJet Cleaning -
    hist_topcleaner.reset(new AndHists(ctx, "Topjet_Cleaning"));
    // - TopJet -
    hist_ntop.reset(new AndHists(ctx, "1TopJetCut"));

    // hist_pileup.reset(new HOTVRPileUpHists(ctx, "HOTVR_PileUp"));
  }

  bool BstarToTWPreSelectionModule::process(Event & event) {

    if(!is_mc)
      {
	if(!sel_lumi->passes(event)) return false;
      }

    hist_nocuts->fill(event);

    // --- Cleaner
    if(!common->process(event)) return false;
    hist_cleaner->fill(event);

    // --- Lepton selection & additional lepton veto
    // - Muon Channel
    if (is_muo)
      {
	if(!(sel_1muo->passes(event) && veto_ele->passes(event))) return false;
	cl_muo->process(event);
	if(!sel_1muo->passes(event)) return false;
	hist_1lep->fill(event);
      }
    // - Electron Channel
    else if (is_ele)
      {
	if(!(sel_1ele->passes(event) && veto_muo->passes(event))) return false;
	cl_ele->process(event);
	if(!sel_1ele->passes(event)) return false;
	hist_1lep->fill(event);
      }
    
    // --- Trigger
    // - Muon
    if (is_muo)
      {
	if(!(trig_IsoMu24->passes(event) || trig_IsoTkMu24->passes(event))) return false;
	hist_trigger->fill(event);
      }
    else if (is_ele)
      {
	if(!(trig_Ele27->passes(event) || trig_Ele115->passes(event))) return false;
	hist_trigger->fill(event);
      }
    
    // --- MET selection
    if(!sel_met->passes(event)) return false;
    hist_met->fill(event);

    // --- ST selection
    if(!sel_st->passes(event)) return false;
    hist_st->fill(event);

    // --- HOTVR JEC
    // hist_pileup->fill(event);
    if (is_mc) 
      jec_subj_mc->process(event);
    else
      {
	if(event.run <= runnr_BCD)         jec_subj_BCD->process(event);
	else if(event.run < runnr_EFearly) jec_subj_EFearly->process(event); // < is correct, not <= 
	else if(event.run <= runnr_FlateG) jec_subj_FlateG->process(event);
	else if(event.run > runnr_FlateG)  jec_subj_H->process(event);
      }
    hist_subj_jec->fill(event);

    // --- TopJet Cleaning
    cl_topjet->process(event);
    hist_topcleaner->fill(event);

    // --- TopJet Selection
    if(!sel_ntop->passes(event)) return false;
    hist_ntop->fill(event);

    // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWPreSelectionModule)

}
