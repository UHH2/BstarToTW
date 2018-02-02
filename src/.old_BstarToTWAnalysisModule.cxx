#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/TTbarReconstruction.h"


#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/BstarToTW/include/HOTVRIds.h"
#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/HOTVRHists.h"
#include "UHH2/BstarToTW/include/CutflowHists.h"
#include "UHH2/BstarToTW/include/BstarToTWReconstruction.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisHists.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"

/* ToDo:
 * - Add b-tag selection
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
    std::unique_ptr<AnalysisModule> BstarToTWgenprod;
    std::unique_ptr<AnalysisModule> primary_lep, bstar_reco, bstar_disc;
    std::string dataset_name;

    // Scale factors
    std::unique_ptr<AnalysisModule> sf_muo_id, sf_muo_iso, sf_muo_trigger;

    // Object Ids
    MuonId id_muo;
    TopJetId id_hotvr, id_hotvr_deltaphi, id_hotvr_lepoverlap;
    TopJetId id_hotvr_only_fpt, id_hotvr_only_mpair, id_hotvr_only_mjet; // Ids for every single HOTVR parameter
    TopJetId id_hotvr_fpt_mpair, id_hotvr_fpt_mjet, id_hotvr_mpair_mjet; // Id for combinations of two HOTVR parameters

    // Cleaner
    std::unique_ptr<AnalysisModule> cl_hotvr_deltaphi, cl_hotvr;

    // Selections
    std::unique_ptr<Selection> sel_ntop, sel_nmuo, sel_vetoele;
    std::unique_ptr<Selection> sel_hotvr_only_fpt, sel_hotvr_only_mpair, sel_hotvr_only_mjet; // selections for every single HOTVR parameter
    std::unique_ptr<Selection> sel_hotvr_fpt_mpair, sel_hotvr_fpt_mjet, sel_hotvr_mpair_mjet; // selections for combinations of two HOTVR parameters
    std::unique_ptr<Selection> sel_hotvr, sel_hotvr_deltaphi, sel_hotvr_lepoverlap;

    // Hists
    std::unique_ptr<AndHists> hist_presel, hist_nmuo, hist_ntop;
    std::unique_ptr<AndHists> hist_only_fpt, hist_only_mpair, hist_only_mjet;
    std::unique_ptr<AndHists> hist_fpt_mpair, hist_fpt_mjet, hist_mpair_mjet; 
    std::unique_ptr<AndHists> hist_hotvr, hist_hotvr_deltaphi, hist_hotvr_lepoverlap;
    std::unique_ptr<CutflowHists> cutflow;
    std::unique_ptr<Hists> hist_bstar_reco;

  };

  BstarToTWAnalysisModule::BstarToTWAnalysisModule(Context & ctx) {

    // --- Setup Additional Modules
    // common modules
    common.reset(new CommonModules());
    common->disable_lumisel();
    common->disable_jersmear();
    common->disable_jec();
    common->init(ctx);

    BstarToTWgenprod.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    dataset_name = ctx.get("dataset_version");

    // --- Scale Factors
    sf_muo_id.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1, "tightID"));
    sf_muo_iso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1, "iso"));
    sf_muo_trigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger"));


    // --- Variables for Selections
    double deltaPhi_min = M_PI/2;  // minimum delta phi between muon and top

    double top_fpt_max   = 0.8;    // maximum pt fraction of leading subjet
    double top_m_min     = 140.;   // minimum jet mass
    double top_m_max     = 220.;   // maximum jet mass
    double top_mpair_min = 50.;    // minimum pairwise mass of first three subjets


    id_hotvr_deltaphi =  AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), DeltaPhiCut(deltaPhi_min));
    id_hotvr_lepoverlap =   AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), LepOverlap(true));

    id_hotvr = HOTVRTopTag(top_fpt_max,
			   top_m_min, 
			   top_m_max, 
			   top_mpair_min);


    // only fpt
    id_hotvr_only_fpt   = HOTVRTopTag(top_fpt_max,
				      0, 
				      FLT_MAX, 
				      0); 

    // only mpair
    id_hotvr_only_mpair = HOTVRTopTag(FLT_MAX,
				      0, 
				      FLT_MAX, 
				      top_mpair_min); 

    // only mjet
    id_hotvr_only_mjet  = HOTVRTopTag(FLT_MAX,
				      top_m_min, 
				      top_m_max, 
				      0); 

    // fpt && mpair
    id_hotvr_fpt_mpair  = HOTVRTopTag(top_fpt_max,
				      0, 
				      FLT_MAX, 
				      top_mpair_min); 

    // fpt && mjet
    id_hotvr_fpt_mjet   =HOTVRTopTag(top_fpt_max,
				     top_m_min, 
				     top_m_max, 
				     0);

    // mpair && mjet
    id_hotvr_mpair_mjet = HOTVRTopTag(FLT_MAX,
				      top_m_min, 
				      top_m_max, 
				      top_mpair_min); 
    //*/

    // --- Cleaner
    cl_hotvr_deltaphi.reset(new TopJetCleaner(ctx, id_hotvr_deltaphi));
    cl_hotvr.reset(new TopJetCleaner(ctx, id_hotvr));

    // --- Selections

    // Kinematic selections

    // HOTVR selections
    sel_ntop.reset(new NTopJetSelection(1, -1, id_hotvr));

    sel_hotvr_only_fpt.reset(new NTopJetSelection(1, -1, id_hotvr_only_fpt));
    sel_hotvr_only_mpair.reset(new NTopJetSelection(1, -1, id_hotvr_only_mpair));
    sel_hotvr_only_mjet.reset(new NTopJetSelection(1, -1, id_hotvr_only_mjet));

    sel_hotvr_fpt_mpair.reset(new NTopJetSelection(1, -1, id_hotvr_fpt_mpair));
    sel_hotvr_fpt_mjet.reset(new NTopJetSelection(1, -1, id_hotvr_fpt_mjet));
    sel_hotvr_mpair_mjet.reset(new NTopJetSelection(1, -1, id_hotvr_mpair_mjet));

    sel_hotvr.reset(new NTopJetSelection(1, -1, id_hotvr));
    sel_hotvr_deltaphi.reset(new NTopJetSelection(1, -1, id_hotvr_deltaphi));
    sel_hotvr_lepoverlap.reset(new NTopJetSelection(1, -1, id_hotvr_lepoverlap));

    // Muon Selections
    sel_nmuo.reset(new NMuonSelection(1, 1));

    // Electron Selection
    sel_vetoele.reset(new NElectronSelection(0, 0));


    // --- Bstar reconstruction
    if (dataset_name.find("BstarToTW") == 0)
      {
	primary_lep.reset(new PrimaryLepton(ctx));
	bstar_reco.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "BstarToTWReconstruction", id_hotvr_fpt_mpair));
	bstar_disc.reset(new BstarToTWMatchDiscriminator(ctx, "BstarToTWReconstruction"));
      }
    hist_bstar_reco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWReco", "BstarToTWReconstruction", "Match"));

    // --- Hists
    hist_presel.reset(new AndHists(ctx, "PreSel"));
    hist_nmuo.reset(new AndHists(ctx, "NMuonCut"));
    hist_ntop.reset(new AndHists(ctx, "NTopCut"));

    hist_only_fpt.reset(new AndHists(ctx, "only_fpt"));
    hist_only_fpt->add_hist(new HOTVRHists(ctx, "only_fpt_HOTVR_tagged", id_hotvr_only_fpt));
    hist_only_mpair.reset(new AndHists(ctx, "only_mpair"));
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
    hist_hotvr_deltaphi.reset(new AndHists(ctx, "DeltaPhiCuts"));
    hist_hotvr_deltaphi->add_hist(new HOTVRHists(ctx, "DeltaPhiCuts_HOTVR_tagged", id_hotvr_deltaphi));
    hist_hotvr_lepoverlap.reset(new AndHists(ctx, "LepOverlap"));
    hist_hotvr_lepoverlap->add_hist(new HOTVRHists(ctx, "LepOverlap_HOTVR_tagged", id_hotvr_lepoverlap));

		     
    cutflow.reset(new CutflowHists(ctx, "Cutflow"));

  }

  bool BstarToTWAnalysisModule::process(Event & event) {

    // after PreSel
    sf_muo_id->process(event);
    sf_muo_iso->process(event);
    sf_muo_trigger->process(event);

    
    // Bstar Reconstructinon
    if (dataset_name.find("BstarToTW") == 0 && sel_hotvr_fpt_mpair->passes(event))
      {
	BstarToTWgenprod->process(event);
	primary_lep->process(event);
	bstar_reco->process(event);
	bstar_disc->process(event);
	hist_bstar_reco->fill(event);
      }

    if(!common->process(event)) return false;
    hist_presel->fill(event);
    cutflow->fill(event, 0);

    if(!sel_nmuo->passes(event)) return false;
    hist_nmuo->fill(event);

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


    // N_Top Cut
    cl_hotvr->process(event);

    if(sel_hotvr_deltaphi->passes(event))
      {
	hist_hotvr_deltaphi->fill(event);
      }

    if(sel_hotvr_lepoverlap->passes(event))
      {
	hist_hotvr_lepoverlap->fill(event);
      }

    cl_hotvr_deltaphi->process(event);
    if(!sel_ntop->passes(event)) return false;

    hist_ntop->fill(event);
    cutflow->fill(event, 9);


     // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWAnalysisModule)

}
