#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/PrimaryLepton.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/BstarToTWModules.h"
#include "UHH2/BstarToTW/include/BstarToTWHists.h"
#include "UHH2/BstarToTW/include/BstarToTWSelections.h"

#include "UHH2/HOTVR/include/HOTVRHists.h"
#include "UHH2/HOTVR/include/HOTVRJetCorrectionModule.h"
#include "UHH2/HOTVR/include/HOTVRJetCorrector.h"
#include "UHH2/HOTVR/include/HOTVRJetCorrectionHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWPreSelectionModule: public AnalysisModule {
  public:
    
    explicit BstarToTWPreSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  

    Year year;
    std::unique_ptr<CommonModules> common_modules;
    std::unique_ptr<AnalysisModule> hotvr_jlc_module, hotvr_jec_module;
    std::unique_ptr<AnalysisModule> primary_lepton, ht_calculator;
    // Cleaner
    std::unique_ptr<AnalysisModule> cl_muo_tight, cl_ele_tight, cl_hotvr;

    // Selections
    bool is_data;
    std::unique_ptr<Selection> sel_lumi;
    std::unique_ptr<Selection> veto_muo;
    std::unique_ptr<Selection> veto_ele;
    std::unique_ptr<Selection> sel_1muo;
    std::unique_ptr<Selection> sel_1ele;
    std::unique_ptr<Selection> sel_njet; 
    std::unique_ptr<Selection> sel_ht;
    std::unique_ptr<Selection> sel_met;
    std::unique_ptr<Selection> sel_st;
    std::unique_ptr<Selection> sel_ntop;

    std::unique_ptr<Selection> sel_deltaPhi_metlep; // delta phi between primary lepton and MET
    std::unique_ptr<AndHists> hist_sel_deltaPhi_metlep;
    std::unique_ptr<Selection> sel_deltaR_jetlep; // delta R between primary lepton and ak4-jets
    std::unique_ptr<AndHists> hist_sel_deltaR_jetlep;


    // Hists
    std::unique_ptr<Hists> hist_cleaner, hist_1lep, hist_njet, hist_ht, hist_met, hist_st, hist_ntop;

    bool is_mc, is_ele, is_muo;

    const bool b_debug = false;

  };

  BstarToTWPreSelectionModule::BstarToTWPreSelectionModule(Context & ctx) {
    if (b_debug) cout << "setting up pre-selection module ..." << endl;

    // --- Config ---
    year = extract_year(ctx);
    is_mc = ctx.get("dataset_type") == "MC";
    is_ele = ctx.get("analysis_channel") == "ELECTRON";
    is_muo = ctx.get("analysis_channel") == "MUON";

    // -- Kinematic Variables -- 
    double jet_pt_min     = 30.0;
    double jet_eta_max    = 2.4;
    double lep_eta_max    = 2.4;
    double lepveto_pt_min = 30.0;
    double lep_pt_min     = 50.0; 
    double muo_iso_max    = 0.15;
    double ht_min         = 200.; 
    double met_min        = 50.;
    double st_min         = 400.;
    double top_pt_min     = 200.;
    double top_eta_max    = 2.5;
    double deltaPhi_metlep_max = M_PI/2; // maximum deltaPhi between missign ET and lepton

    // --- IDs ---
    JetId id_jet = PtEtaCut(jet_pt_min, jet_eta_max);
    MuonId id_muo_veto = AndId<Muon>(MuonID(Muon::Selector::CutBasedIdLoose), PtEtaCut(lepveto_pt_min, lep_eta_max), MuonIso(muo_iso_max)); // muon veto ID
    MuonId id_muo_tight = AndId<Muon>(MuonID(Muon::Selector::CutBasedIdTight), PtEtaCut(lep_pt_min, lep_eta_max), MuonIso(muo_iso_max)); // muon ID
    ElectronId id_ele_veto = AndId<Electron>(ElectronID_Fall17_loose, PtEtaCut(lepveto_pt_min, lep_eta_max)); // electron veto ID
    ElectronId id_ele_tight = AndId<Electron>(ElectronID_Fall17_tight, PtEtaCut(lep_pt_min, lep_eta_max)); // electron ID
    // if (year == Year::is2016v2 || year == Year::is2016v3)
    //   {
    if (year == Year::is2016v2)
      {
	id_muo_veto = AndId<Muon>(MuonID(Muon::Selector::Loose), PtEtaCut(lepveto_pt_min, lep_eta_max), MuonIso(muo_iso_max));
	id_muo_tight = AndId<Muon>(MuonID(Muon::Selector::Tight), PtEtaCut(lep_pt_min, lep_eta_max), MuonIso(muo_iso_max));
      }
    // 	id_ele_veto = AndId<Electron>(ElectronID_Summer16_veto, PtEtaCut(lepveto_pt_min, lep_eta_max));
    // 	id_ele_tight = AndId<Electron>(ElectronID_Summer16_tight, PtEtaCut(lep_pt_min, lep_eta_max)); // electron ID
    //   }
    // else
    //   {
    // 	id_ele_veto = AndId<Electron>(ElectronID_Fall17_veto, PtEtaCut(lepveto_pt_min, lep_eta_max));
    // 	id_ele_tight = AndId<Electron>(ElectronID_Fall17_tight, PtEtaCut(lep_pt_min, lep_eta_max)); // electron ID
    //   }
    TopJetId id_hotvr = PtEtaCut(top_pt_min, top_eta_max);

    // --- Additional Modules ---
    common_modules.reset(new CommonModules());
    common_modules->switch_jetlepcleaner(true);
    common_modules->switch_jetPtSorter(true);
    common_modules->switch_metcorrection(true);
    common_modules->set_jet_id(id_jet);
    common_modules->set_muon_id(id_muo_veto);
    common_modules->set_electron_id(id_ele_veto);
    common_modules->init(ctx);

    hotvr_jlc_module.reset(new HOTVRJetLeptonCleaner());
    hotvr_jec_module.reset(new HOTVRJetCorrectionModule(ctx));

    primary_lepton.reset(new PrimaryLepton(ctx));
    ht_calculator.reset(new HTCalculator(ctx));
    
    // Cleaner
    cl_muo_tight.reset(new MuonCleaner(id_muo_tight));
    cl_ele_tight.reset(new ElectronCleaner(id_ele_tight));
    cl_hotvr.reset(new TopJetCleaner(ctx, id_hotvr));
    // --- Selections ---
    // - Lumi -
    sel_lumi.reset(new LumiSelection(ctx));
    // - Lepton vetos
    veto_muo.reset(new NMuonSelection(0, 0));
    veto_ele.reset(new NElectronSelection(0, 0));
    // - Lepton
    sel_1muo.reset(new NMuonSelection(1, 1));
    sel_1ele.reset(new NElectronSelection(1, 1));
    // - JET
    sel_njet.reset(new NJetSelection(1, -1));
    // - HT
    sel_ht.reset(new HTSelection(ctx, ht_min));
    // - MET
    sel_met.reset(new METSelection(met_min));
    // - ST
    sel_st.reset(new STSelection(ctx, st_min));
    // - deltaPhi met-lep
    sel_deltaPhi_metlep.reset(new METDeltaPhiSelection(ctx, deltaPhi_metlep_max));
    hist_sel_deltaPhi_metlep.reset(new AndHists(ctx, "deltaPhiMETLep"));
    // - deltaR jet-lep
    sel_deltaR_jetlep.reset(new JetDeltaRSelection(ctx, 0.4));
    hist_sel_deltaR_jetlep.reset(new AndHists(ctx, "deltaR_jetlep"));
    // - TopJet
    sel_ntop.reset(new NTopJetSelection(1, -1));

    // --- Hists ---
    // - Cleaning -
    hist_cleaner.reset(new AndHists(ctx, "Cleaning"));
    // - Lepton -
    hist_1lep.reset(new AndHists(ctx, "1LepCut"));
    // - Jet -
    hist_njet.reset(new AndHists(ctx, "1JetCut"));
    // - HT -
    hist_ht.reset(new AndHists(ctx, "100HTCut"));
    // - MET -
    hist_met.reset(new AndHists(ctx, "50METCut"));
    // - ST -
    hist_st.reset(new AndHists(ctx, "400STCut"));
    // - TopJet -
    hist_ntop.reset(new AndHists(ctx, "1TopJetCut"));

  }

  bool BstarToTWPreSelectionModule::process(Event & event) {

    
    // --- Object Setup and Cleaning
    if (b_debug) cout << "processing common moduels ..." << endl;
    if (!common_modules->process(event)) return false;
    if (b_debug) cout << "processing primary lepton ..." << endl;
    primary_lepton->process(event);
    if (b_debug) cout << "processing ht calculator ..." << endl;
    ht_calculator->process(event);
    if (b_debug) cout << "done, fill hists ..." << endl;
    hist_cleaner->fill(event);

    // --- Lepton selection & additional lepton veto
    if (b_debug) cout << "processing lepton selection ..." << endl;
    // - Muon Channel
    if (is_muo)
      {
	if(!(sel_1muo->passes(event) && veto_ele->passes(event))) return false;
	cl_muo_tight->process(event);
	if(!sel_1muo->passes(event)) return false;
	if (b_debug) cout << "done, fill hists ..." << endl;
	hist_1lep->fill(event);
      }    
    // - Electron Channel
    else if (is_ele)
      {
	if(!(sel_1ele->passes(event) && veto_muo->passes(event))) return false;
	cl_ele_tight->process(event);
	if(!sel_1ele->passes(event)) return false;
	if (b_debug) cout << "done, fill hists ..." << endl;
	hist_1lep->fill(event);
      }

    // --- Jet selection
    if (b_debug) cout << "processing ak4-jet selection ..." << endl;
    if (!sel_njet->passes(event)) return false;
    if (b_debug) cout << "done, fill hists ..." << endl;
    hist_njet->fill(event);

    // --- HT selection
    if (b_debug) cout << "processing HT selection ..." << endl;
    if(!sel_ht->passes(event)) return false;
    if (b_debug) cout << "done, fill hists ..." << endl;
    hist_ht->fill(event);

    // --- MET selection
    if (b_debug) cout << "processing MET selection ..." << endl;
    if(!sel_met->passes(event)) return false;
    if (b_debug) cout << "done, fill hists ..." << endl;
    hist_met->fill(event);

    // --- ST selection
    if (b_debug) cout << "processing ST selection ..." << endl;
    if(!sel_st->passes(event)) return false;
    if (b_debug) cout << "done, fill hists ..." << endl;
    hist_st->fill(event);

    // deltaR jet-lep
    if (b_debug) cout << "processing deltaR jet-lepton selection ..." << endl;
    if(!sel_deltaR_jetlep->passes(event)) return false;
    if (b_debug) cout << "done, fill hists ..." << endl;
    hist_sel_deltaR_jetlep->fill(event);

    // deltaPhi met-lep
    if (b_debug) cout << "processing deltaPhi MET-lepton selection ..." << endl;
    if(!sel_deltaPhi_metlep->passes(event)) return false;
    if (b_debug) cout << "done, fill hists ..." << endl;
    hist_sel_deltaPhi_metlep->fill(event);

    // --- TopJet Selection
    if (b_debug) cout << "processing HOTVR JEC/JER ..." << endl;
    hotvr_jlc_module->process(event);
    hotvr_jec_module->process(event);
    cl_hotvr->process(event);
    if (b_debug) cout << "processing HOTVR-jet selection ..." << endl;
    if(!sel_ntop->passes(event)) return false;
    if (b_debug) cout << "done, fill hists ..." << endl;
    hist_ntop->fill(event);
     
    // done
    if (b_debug) cout << "done!" << endl;
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWPreSelectionModule)

}
