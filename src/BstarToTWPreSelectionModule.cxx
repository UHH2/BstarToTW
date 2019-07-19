#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

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
#include "UHH2/BstarToTW/include/BstarToTWCommonModules.h"
#include "UHH2/BstarToTW/include/BstarToTWHists.h"
#include "UHH2/BstarToTW/include/BstarToTWSelections.h"

#include "UHH2/HOTVR/include/HOTVRHists.h"
#include "UHH2/HOTVR/include/HOTVRJetCorrectionHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWPreSelectionModule: public AnalysisModule {
  public:
    
    explicit BstarToTWPreSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  

    std::unique_ptr<AnalysisModule> primary_lepton, ht_calculator;
    std::unique_ptr<AnalysisModule> object_setup;
    std::unique_ptr<AnalysisModule> sf_mc_lumiweight, sf_mc_pu_reweight;
    // Cleaner
    std::unique_ptr<AnalysisModule> cl_muo_tight, cl_ele_tight;

    // Selections
    bool is_data;
    std::unique_ptr<Selection> sel_lumi;
    std::unique_ptr<Selection> sel_trigger;
    std::unique_ptr<Selection> veto_muo;
    std::unique_ptr<Selection> veto_ele;
    std::unique_ptr<Selection> sel_1muo;
    std::unique_ptr<Selection> sel_1ele;
    std::unique_ptr<Selection> sel_njet;
    std::unique_ptr<Selection> sel_met;
    std::unique_ptr<Selection> sel_st;
    std::unique_ptr<Selection> sel_ntop;

    // Hists
    std::unique_ptr<Hists> hist_trigger, hist_cleaner, hist_1lep, hist_njet, hist_met, hist_st, hist_ntop;

    bool is_mc, is_ele, is_muo;
  };

  BstarToTWPreSelectionModule::BstarToTWPreSelectionModule(Context & ctx) {

    // --- Config ---
    is_mc = ctx.get("dataset_type") == "MC";
    is_ele = ctx.get("analysis_channel") == "ELECTRON";
    is_muo = ctx.get("analysis_channel") == "MUON";

    // --- Additional Modules ---
    primary_lepton.reset(new PrimaryLepton(ctx));
    ht_calculator.reset(new HTCalculator(ctx));
    object_setup.reset(new ObjectSetup(ctx));
    if (is_mc)
      {
	sf_mc_lumiweight.reset(new MCLumiWeight(ctx));
	std::string syst_pu = "central";
	sf_mc_pu_reweight.reset(new MCPileupReweight(ctx, syst_pu));
      }

    // -- Kinematic Variables --
    double met_min = 50.;
    double st_min = 400.;
    double lep_eta_max = 2.4;
    double lep_pt_min  = 50.0; 
    double muo_iso_max = 0.15;

    // --- IDs ---
    MuonId id_muo_tight = AndId<Muon>(MuonID(Muon::Selector::CutBasedIdTight), PtEtaCut(lep_pt_min, lep_eta_max), MuonIso(muo_iso_max)); // muon ID
    ElectronId id_ele_tight = AndId<Electron>(ElectronID_Summer16_tight, PtEtaCut(lep_pt_min, lep_eta_max)); // electron ID

    // Cleaner
    cl_muo_tight.reset(new MuonCleaner(id_muo_tight));
    cl_ele_tight.reset(new ElectronCleaner(id_ele_tight));
    
    // --- Selections ---
    // - Lumi -
    sel_lumi.reset(new LumiSelection(ctx));
    // - Trigger -
    sel_trigger.reset(new BstarToTWTriggerSelection(ctx));
    // - Lepton vetos
    veto_muo.reset(new NMuonSelection(0, 0));
    veto_ele.reset(new NElectronSelection(0, 0));
    // - Lepton
    sel_1muo.reset(new NMuonSelection(1, 1));
    sel_1ele.reset(new NElectronSelection(1, 1));
    // - JET
    sel_njet.reset(new NJetSelection(1, -1));
    // - MET
    sel_met.reset(new METSelection(met_min));
    // - ST
    sel_st.reset(new STSelection(st_min));
    // - TopJet
    sel_ntop.reset(new NTopJetSelection(1, -1));

    // --- Hists ---
    // - Cleaning -
    hist_cleaner.reset(new AndHists(ctx, "Cleaning"));
    // - Lepton -
    hist_1lep.reset(new AndHists(ctx, "1LepCut"));
    // - Jet -
    hist_njet.reset(new AndHists(ctx, "1JetCut"));
    // - Trigger -
    hist_trigger.reset(new AndHists(ctx, "Trigger"));
    // - MET -
    hist_met.reset(new AndHists(ctx, "50METCut"));
    // - ST -
    hist_st.reset(new AndHists(ctx, "400STCut"));
    // - TopJet -
    hist_ntop.reset(new AndHists(ctx, "1TopJetCut"));

  }

  bool BstarToTWPreSelectionModule::process(Event & event) {
    
    // --- Object Setup and Cleaning
    if(!is_mc)
      {
	if(!sel_lumi->passes(event)) return false;
      }
    else if(is_mc)
      {
	sf_mc_lumiweight->process(event);
	sf_mc_pu_reweight->process(event);
      }
    
    if(!object_setup->process(event)) return false;
    primary_lepton->process(event);
    ht_calculator->process(event);

    hist_cleaner->fill(event);

    // --- Lepton selection & additional lepton veto
    // - Muon Channel
    if (is_muo)
      {
	if(!(sel_1muo->passes(event) && veto_ele->passes(event))) return false;
	cl_muo_tight->process(event);
	if(!sel_1muo->passes(event)) return false;
	hist_1lep->fill(event);
      }
    
    // - Electron Channel
    else if (is_ele)
      {
	if(!(sel_1ele->passes(event) && veto_muo->passes(event))) return false;
	cl_ele_tight->process(event);
	if(!sel_1ele->passes(event)) return false;
	hist_1lep->fill(event);
      }

    // --- Jet selection
    if (!sel_njet->passes(event)) return false;
    hist_njet->fill(event);

    // --- MET selection
    if(!sel_met->passes(event)) return false;
    hist_met->fill(event);

    // --- ST selection
    if(!sel_st->passes(event)) return false;
    hist_st->fill(event);

    // --- TopJet Selection
    if(!sel_ntop->passes(event)) return false;
    hist_ntop->fill(event);

    if (!sel_trigger->passes(event)) return false;
    hist_trigger->fill(event);
     
    // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWPreSelectionModule)

}
