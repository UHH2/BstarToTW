#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"

#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/BstarToTW/include/HOTVRIds.h"
#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/HOTVRHists.h"
#include "UHH2/BstarToTW/include/CutflowHists.h"

/* ToDo:
 * - Add b-tag selection
 *
 */


using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWAnalysisModule: public AnalysisModule {
  public:
    
    explicit BstarToTWAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  

    // additional Modules
    std::unique_ptr<CommonModules> common;

    // Scale factors
    std::unique_ptr<AnalysisModule> sf_muo_id, sf_muo_iso, sf_muo_trigger;

    // Object Ids
    MuonId id_muo;
    TopJetId id_hotvr, id_hotvr_deltaphi;
    TopJetId id_hotvr_only_fpt, id_hotvr_only_mpair, id_hotvr_only_mjet; // Ids for every single HOTVR parameter
    TopJetId id_hotvr_fpt_mpair, id_hotvr_fpt_mjet, id_hotvr_mpair_mjet; // Id for combinations of two HOTVR parameters

    // Selections
    std::unique_ptr<Selection> sel_met, sel_ntop;
    std::unique_ptr<Selection> sel_hotvr_only_fpt, sel_hotvr_only_mpair, sel_hotvr_only_mjet; // selections for every single HOTVR parameter
    std::unique_ptr<Selection> sel_hotvr_fpt_mpair, sel_hotvr_fpt_mjet, sel_hotvr_mpair_mjet; // selections for combinations of two HOTVR parameters
    std::unique_ptr<Selection> sel_hotvr, sel_hotvr_deltaphi;

    // Hists
    std::unique_ptr<AndHists> hist_presel, hist_metcut, hist_ntop;
    std::unique_ptr<AndHists> hist_only_fpt, hist_only_mpair, hist_only_mjet;
    std::unique_ptr<AndHists> hist_fpt_mpair, hist_fpt_mjet, hist_mpair_mjet; 
    std::unique_ptr<AndHists> hist_hotvr, hist_hotvr_deltaphi;
    std::unique_ptr<CutflowHists> cutflow;

  };

  BstarToTWAnalysisModule::BstarToTWAnalysisModule(Context & ctx) {

    // --- Setup Additional Modules
    common.reset(new CommonModules());
    common->disable_jec();
    common->disable_jersmear();
    common->init(ctx);

    // --- Scale Factors
    sf_muo_id.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1, "tightID"));
    sf_muo_iso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1, "tightID"));
    sf_muo_trigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger"));

    // --- Variables for Selections
    double muo_pt_min  = 130.;     // minimum muon pt
    double muo_eta_max = 2.4;      // maximum muon eta
    id_muo = AndId<Muon>(MuonIDTight(), PtEtaCut(muo_pt_min, muo_eta_max));

    double met_min = 60.;          // minimum missing transvere energy

    double deltaPhi_min = M_PI/2;  // minimum delta phi between muon and top

    double top_pt_min    = 200.;   // minimum top jet pt
    double top_eta_max   = 2.4;    // maximum top jet eta
    double top_fpt_max   = 0.8;    // maximum pt fraction of leading subjet
    double top_m_min     = 140.;   // minimum jet mass
    double top_m_max     = 220.;   // maximum jet mass
    double top_mpair_min = 50.;    // minimum pairwise mass of first three subjets

    id_hotvr = AndId<TopJet>(HOTVRTopTag(top_fpt_max,
					 top_m_min, 
					 top_m_max, 
					 top_mpair_min),
			     PtEtaCut(top_pt_min, top_eta_max));

    id_hotvr_deltaphi  = AndId<TopJet>(HOTVRTopTag(top_fpt_max,
						   top_m_min, 
						   top_m_max, 
						   top_mpair_min),
				       DeltaPhiCut(deltaPhi_min),
				       PtEtaCut(top_pt_min, top_eta_max));

    // only fpt
    id_hotvr_only_fpt   = AndId<TopJet>(HOTVRTopTag(top_fpt_max,
						    0, 
						    FLT_MAX, 
						    0),
					PtEtaCut(top_pt_min, top_eta_max)); 

    // only mpair
    id_hotvr_only_mpair = AndId<TopJet>(HOTVRTopTag(FLT_MAX,
						    0, 
						    FLT_MAX, 
						    top_mpair_min),
					PtEtaCut(top_pt_min, top_eta_max)); 

    // only mjet
    id_hotvr_only_mjet  = AndId<TopJet>(HOTVRTopTag(FLT_MAX,
						    top_m_min, 
						    top_m_max, 
						    0),
					PtEtaCut(top_pt_min, top_eta_max)); 

    // fpt && mpair
    id_hotvr_fpt_mpair  = AndId<TopJet>(HOTVRTopTag(top_fpt_max,
						    0, 
						    FLT_MAX, 
						    top_mpair_min),
					PtEtaCut(top_pt_min, top_eta_max)); 

    // fpt && mjet
    id_hotvr_fpt_mjet   = AndId<TopJet>(HOTVRTopTag(top_fpt_max,
						    top_m_min, 
						    top_m_max, 
						    0),
					PtEtaCut(top_pt_min, top_eta_max));

    // mpair && mjet
    id_hotvr_mpair_mjet = AndId<TopJet>(HOTVRTopTag(FLT_MAX,
						    top_m_min, 
						    top_m_max, 
						    top_mpair_min),
					PtEtaCut(top_pt_min, top_eta_max)); 
    //*/

    // --- Selections
    // MET selection
    sel_met.reset(new METSelection(met_min));

    // Kinematic selections

    // HOTVR selections
    sel_ntop.reset(new NTopJetSelection(1, 1, id_hotvr));

    sel_hotvr_only_fpt.reset(new NTopJetSelection(1, -1, id_hotvr_only_fpt));
    sel_hotvr_only_mpair.reset(new NTopJetSelection(1, -1, id_hotvr_only_mpair));
    sel_hotvr_only_mjet.reset(new NTopJetSelection(1, -1, id_hotvr_only_mjet));

    sel_hotvr_fpt_mpair.reset(new NTopJetSelection(1, -1, id_hotvr_fpt_mpair));
    sel_hotvr_fpt_mjet.reset(new NTopJetSelection(1, -1, id_hotvr_fpt_mjet));
    sel_hotvr_mpair_mjet.reset(new NTopJetSelection(1, -1, id_hotvr_mpair_mjet));

    sel_hotvr.reset(new NTopJetSelection(1, -1, id_hotvr));
    sel_hotvr_deltaphi.reset(new NTopJetSelection(1, -1, id_hotvr_deltaphi));


    // --- Hists
    hist_presel.reset(new AndHists(ctx, "PreSel"));
    hist_metcut.reset(new AndHists(ctx, "METCut"));
    hist_ntop.reset(new AndHists(ctx, "NTopCut"));

    hist_only_fpt.reset(new AndHists(ctx, "only_fpt"));
    hist_only_fpt->add_hist(new HOTVRHists(ctx, "only_fpt_HOTVR_tagged", id_hotvr_only_fpt));
    hist_only_mpair.reset(new AndHists(ctx, "only_mpai"));
    hist_only_mpair->add_hist(new HOTVRHists(ctx, "only_mpair_HOTVR_tagged", id_hotvr_only_mpair));
    hist_only_mjet.reset(new AndHists(ctx, "only_mjet"));
    hist_only_mjet->add_hist(new HOTVRHists(ctx, "only_mjet_HOTVR_tagged", id_hotvr_only_mjet));

    hist_fpt_mpair.reset(new AndHists(ctx, "fpt_mpair"));
    hist_fpt_mpair->add_hist(new HOTVRHists(ctx, "fpt_mapir_HOTVR_tagged", id_hotvr_fpt_mpair));
    hist_fpt_mjet.reset(new AndHists(ctx, "fpt_mjet"));
    hist_fpt_mjet->add_hist(new HOTVRHists(ctx, "fpt_mjet_HOTVR_tagged", id_hotvr_fpt_mjet));
    hist_mpair_mjet.reset(new AndHists(ctx, "mpair_mjet"));
    hist_mpair_mjet->add_hist(new HOTVRHists(ctx, "mpair_mjet_HOTVR_tagged", id_hotvr_mpair_mjet));

    hist_hotvr.reset(new AndHists(ctx, "HOTVRCuts"));
    hist_hotvr->add_hist(new HOTVRHists(ctx, "HOTVRCuts_HOTVR_tagged", id_hotvr));
    hist_hotvr_deltaphi.reset(new AndHists(ctx, "DeltaPhhiCuts"));
    hist_hotvr_deltaphi->add_hist(new HOTVRHists(ctx, "DeltaPhiCuts_HOTVR_tagged", id_hotvr));
		     
    cutflow.reset(new CutflowHists(ctx, "Cutflow"));

  }

  bool BstarToTWAnalysisModule::process(Event & event) {

    // after PreSel
    sf_muo_id->process(event);
    sf_muo_iso->process(event);
    sf_muo_trigger->process(event);

    // if(!cl_pteta->process(event)) return false; // all events should pass this.
    if(!common->process(event)) return false;
    hist_presel->fill(event);
    cutflow->fill(event, 0);


    // MET Cut
    if (!sel_met->passes(event)) return false;
    hist_metcut->fill(event);
    cutflow->fill(event, 1);

    // HOTVR Cuts
    if(sel_hotvr_only_fpt->passes(event))
      {
	hist_only_fpt->fill(event);
	cutflow->fill(event, 2);
      }

    if(sel_hotvr_only_mpair->passes(event))
      {
	hist_only_mpair->fill(event);
	cutflow->fill(event, 3);
      }

    if(sel_hotvr_only_mjet->passes(event))
      {
	hist_only_mjet->fill(event);
	cutflow->fill(event, 4);
      }

    if(sel_hotvr_fpt_mpair->passes(event))
      {
	hist_fpt_mpair->fill(event);
	cutflow->fill(event, 5);
      }

    if(sel_hotvr_fpt_mjet->passes(event))
      {
	hist_fpt_mjet->fill(event);
	cutflow->fill(event, 6);
      }

    if(sel_hotvr_mpair_mjet->passes(event))
      {
	hist_mpair_mjet->fill(event);
	cutflow->fill(event, 7);
      }

    if(sel_hotvr->passes(event))
      {
	hist_hotvr->fill(event);
	cutflow->fill(event, 8);
      }

    if(sel_hotvr_deltaphi->passes(event))
      {
	hist_hotvr_deltaphi->fill(event);
	cutflow->fill(event, 8);
      }

    // N_Top Cut
    if(sel_ntop->passes(event))
      {
	hist_ntop->fill(event);
	cutflow->fill(event, 9);
      }



     // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWAnalysisModule)

}
