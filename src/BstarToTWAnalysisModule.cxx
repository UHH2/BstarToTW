#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/EventHists.h"

#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/BstarToTW/include/TopTagIndexer.h"
#include "UHH2/BstarToTW/include/IndexerTools.h"
#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/HOTVRHists.h"


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
    std::unique_ptr<AnalysisModule> indexerprod;

    // Scale factors
    std::unique_ptr<AnalysisModule> sf_muo_id, sf_muo_iso, sf_muo_trigger;

    // Index Cleaner
    std::unique_ptr<AnalysisModule> cl_pteta, cl_deltaphi, cl_fpt, cl_nsub, cl_mtop, cl_mpair;

    // Selections
    std::unique_ptr<Selection> trig_IsoMu24, trig_IsoTkMu24, sel_met, sel_ptmuo, sel_ntop;

    // Hists
    std::unique_ptr<AndHists> hist_common; //for testing
    std::unique_ptr<AndHists> hist_presel, hist_trigger, hist_ptmuo, hist_metcut, hist_deltaphi;
    std::unique_ptr<AndHists> hist_fpt, hist_nsub, hist_mtop, hist_mpair, hist_ntop;

    MuonId id_muo;
  };

  BstarToTWAnalysisModule::BstarToTWAnalysisModule(Context & ctx) {

    // setup additional modules
    common.reset(new CommonModules());
    common->disable_jec();
    common->disable_jersmear();
    common->init(ctx);
    indexerprod.reset(new TopTagIndexerProducer(ctx, "TopTagIndexer"));

    // Scale factors
    sf_muo_id.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1, "tightID"));
sf_muo_iso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1, "tightID"));
sf_muo_trigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger"));

    // variables for selections
    double muo_pt_min  = 130.; 
    double muo_eta_max = 2.4;
    id_muo = AndId<Muon>(MuonIDTight(), PtEtaCut(muo_pt_min, muo_eta_max));

    double met_min            = 80.;

    double deltaPhi_min       = M_PI/2; // minimum deltaPhi between Muon and Top

    double top_fpt_max        = 0.8;
    unsigned int top_nsub_min = 3;
    double top_m_min          = 140.;
    double top_m_max          = 220.;
    double top_mpair_min      = 50.;

    // --- (index) Cleaner
    cl_pteta.reset(new PtEtaTopIndexCleaner(ctx, 200., 2.4));

    // --- Selections
    // muon selections
    trig_IsoMu24.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    trig_IsoTkMu24.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    sel_ptmuo.reset(new NMuonSelection(1, -1, id_muo));

    // MET selection (TODO: Check if this one can go to PreSel)
    sel_met.reset(new METSelection(met_min));

    // Kinematic selections
    cl_deltaphi.reset(new DeltaPhiTopIndexCleaner(ctx, deltaPhi_min));

    // HOTVR selections
    cl_fpt.reset(new FptTopIndexCleaner(ctx, top_fpt_max));
    cl_nsub.reset(new NSubTopIndexCleaner(ctx, top_nsub_min));
    cl_mtop.reset(new MTopIndexCleaner(ctx, top_m_min, top_m_max));
    cl_mpair.reset(new MpairTopIndexCleaner(ctx, top_mpair_min));
    sel_ntop.reset(new NTopIndexSelection(ctx, 1, 1));

    // --- Hists
    hist_common.reset(new AndHists(ctx, "CommonModule"));

    hist_presel.reset(new AndHists(ctx, "PreSel"));
    hist_trigger.reset(new AndHists(ctx, "Trigger"));
    hist_ptmuo.reset(new AndHists(ctx, "PtMuoCut"));
    hist_metcut.reset(new AndHists(ctx, "METCut"));
    hist_deltaphi.reset(new AndHists(ctx, "DeltaPhiCut"));
    hist_fpt.reset(new AndHists(ctx, "FptCut"));
    hist_nsub.reset(new AndHists(ctx, "NSubCut"));
    hist_mtop.reset(new AndHists(ctx, "MTopCut"));
    hist_mpair.reset(new AndHists(ctx, "MPairwiseCut"));
    hist_ntop.reset(new AndHists(ctx, "NTopCut"));

  }

  bool BstarToTWAnalysisModule::process(Event & event) {

    // after PreSel
    sf_muo_id->process(event);
    sf_muo_iso->process(event);

    indexerprod->process(event);
    if(!cl_pteta->process(event)) return false; // all events should pass this.
    hist_presel->fill(event);

    if(!common->process(event)) return false;
    hist_common->fill(event);

    // Muon Trigger
    if(!trig_IsoMu24->passes(event) || !trig_IsoTkMu24->passes(event)) return false;
    sf_muo_trigger->process(event);
    hist_trigger->fill(event);

    // Muon Pt Cut
    if(!sel_ptmuo->passes(event)) return false;
    hist_ptmuo->fill(event);

    // MET Cut
    if (!sel_met->passes(event)) return false;
    hist_metcut->fill(event);

    //Delta Phi Cut
    if (!cl_deltaphi->process(event)) return false;
    hist_deltaphi->fill(event);

    // HOTVR recommended Cuts
    if (!cl_fpt->process(event)) return false;
    hist_fpt->fill(event);

    if (!cl_nsub->process(event)) return false;
    hist_nsub->fill(event);

    if (!cl_mtop->process(event)) return false;
    hist_mtop->fill(event);

    if (!cl_mpair->process(event)) return false;
    hist_mpair->fill(event);

    // N Top cut
    if(!sel_ntop->passes(event)) return false;
    hist_ntop->fill(event);

     // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWAnalysisModule)

}
