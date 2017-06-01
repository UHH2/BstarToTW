#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/JetHists.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/BstarToTW/include/BstarToTWReconstruction.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisHists.h"
#include "UHH2/BstarToTW/include/AndHists.h"

/* ToDo:
 * - Add chi2 selection
 */


using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWAnalysisModule_AK8: public AnalysisModule {
  public:
    
    explicit BstarToTWAnalysisModule_AK8(Context & ctx);
    virtual bool process(Event & event) override;

  private:  
    std::string dataset_name;

    // --- Additional Modules
    // Common Modules
    std::unique_ptr<CommonModules> common;
    // Bstar Reconsturction
    std::unique_ptr<AnalysisModule> BstarToTWgenprod;
    std::unique_ptr<AnalysisModule> primary_lep, bstar_reco, bstar_disc, bstar_matchdisc;

    // Scale factors
    std::unique_ptr<AnalysisModule> sf_muo_id, sf_muo_iso, sf_muo_trigger, sf_btag;

    // --- Cleaner
    std::unique_ptr<AnalysisModule> cl_toptag, cl_toplep;

    // --- Selections
    std::unique_ptr<Selection> sel_toptag, veto_ele, sel_muo, sel_btag;

    // --- Histograms
    std::unique_ptr<Hists> hist_presel, hist_sel_toptag, hist_veto_ele, hist_sel_muo, hist_sel_btag;

    std::unique_ptr<Hists> hist_bstar_reco, hist_bstar_matchreco;
    std::unique_ptr<Hists> hist_BTagMCEfficiency;

    // Scans
    std::vector<std::unique_ptr<Selection>> sel_chi2;
    std::vector<std::unique_ptr<AndHists>> hist_chi2;



  };

  BstarToTWAnalysisModule_AK8::BstarToTWAnalysisModule_AK8(Context & ctx) {
    dataset_name = ctx.get("dataset_version");

    // --- Kinematic Variables
    double deltaPhi_min = M_PI/2;  // minimum delta phi between muon and top

    double top_m_min = 140.;   // minimum jet mass
    double top_m_max = 220.;   // maximum jet mass
    double top_tau32 = 0.68;   // maximum nsubjetiness

    CSVBTag::wp btag_wp = CSVBTag::WP_MEDIUM;

    // --- Object IDs
   


    // --- Setup Additional Modules
    // Common Modules    
    common.reset(new CommonModules());
    common->disable_jec();
    common->disable_jersmear();
    common->init(ctx);

    // // Bstar reconstruction

    // primary_lep.reset(new PrimaryLepton(ctx));
    // bstar_reco.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "BstarToTWReconstruction", id_toptag));
    // bstar_disc.reset(new BstarToTWChi2Discriminator(ctx, "BstarToTWReconstruction"));
    // if (dataset_name.find("BstarToTW") == 0)
    //   {
    // 	BstarToTWgenprod.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    // 	bstar_matchdisc.reset(new BstarToTWMatchDiscriminator(ctx, "BstarToTWReconstruction"));
    //   }

    // // Scale Factors
    // sf_muo_id.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1, "tightID"));
    // sf_muo_iso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1, "iso"));
    // sf_muo_trigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger"));
    // // sf_btag.reset(new MCBTagScaleFactor(ctx, btag_wp));
   
    // // --- Cleaner
 
    // // --- Selections
    // // TopTag
    // sel_toptag.reset(new NTopJetSelection(1, 1));

    // veto_ele.reset(new NElectronSelection(0, 0));
    // sel_muo.reset(new NMuonSelection(1, 1));
    // sel_btag.reset(new NJetSelection(1, 1, id_btag));

    // // --- Histograms
    // hist_presel.reset(new AndHists(ctx, "PreSel"));
    // hist_sel_toptag.reset(new AndHists(ctx, "1TopTag"));
    // hist_veto_ele.reset(new AndHists(ctx, "ElectronVeto"));
    // hist_sel_muo.reset(new AndHists(ctx, "1Muon"));
    // hist_sel_btag.reset(new AndHists(ctx, "1BtagLoose"));

    // hist_BTagMCEfficiency.reset(new BTagMCEfficiencyHists(ctx, "BTagMCEfficiency", btag_wp));

    // hist_bstar_reco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWReco", "BstarToTWReconstruction", "Chi2"));
    // hist_bstar_matchreco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWMatchedReco", "BstarToTWReconstruction", "Match"));

  }

  bool BstarToTWAnalysisModule_AK8::process(Event & event) {

    // // after PreSel
    // // apply scale factors for muons and trigger
 
    // sf_muo_id->process(event);
    // sf_muo_iso->process(event);
    // sf_muo_trigger->process(event);
    // common->process(event);
    
    // hist_presel->fill(event);

    // // TopTag Cut
    // cl_toplep->process(event);
    // cl_toptag->process(event);
    
    // if(!sel_toptag->passes(event)) return false;
    // hist_sel_toptag->fill(event);

    // //Bstar Reconstructinon   	
    // primary_lep->process(event);
    // bstar_reco->process(event);
    // bstar_disc->process(event);
    // hist_bstar_reco->fill(event);
    // if (dataset_name.find("BstarToTW") == 0)
    //   {
    // 	BstarToTWgenprod->process(event);

    // 	bstar_matchdisc->process(event);
    // 	hist_bstar_matchreco->fill(event);
    //   }

    // for (int i = 0; i < 6; ++i)
    //   {
    // 	if (sel_chi2.at(i)->passes(event))
    // 	  {
    // 	    hist_chi2.at(i)->fill(event);
    // 	  }
    //   }

    // // Electron Veto
    // if (!veto_ele->passes(event)) return false;
    // hist_veto_ele->fill(event);

    // // Muon Cut
    // if (!sel_muo->passes(event)) return false;
    // hist_sel_muo->fill(event);

    // // btag Cut
    // hist_BTagMCEfficiency->fill(event);
    // if (!sel_btag->passes(event)) return false;
    // // sf_btag->process(event);
    // hist_sel_btag->fill(event);
   
    // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWAnalysisModule_AK8)

}
