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

#include "UHH2/BstarToTW/include/HOTVRIds.h"
#include "UHH2/BstarToTW/include/HOTVRJetCorrectionHists.h"

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

  class BstarToTWAnalysisModule: public AnalysisModule {
  public:
    
    explicit BstarToTWAnalysisModule(Context & ctx);
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
    std::unique_ptr<Selection> sel_tau32_cr;

    // --- Histograms
    std::unique_ptr<Hists> hist_hotvr_jec;
    std::unique_ptr<Hists> hist_presel, hist_sel_toptag, hist_veto_ele, hist_sel_muo, hist_sel_btag;

    std::unique_ptr<Hists> hist_bstar_reco, hist_bstar_matchreco;
    std::unique_ptr<Hists> hist_BTagMCEfficiency, hist_tau32MCEfficiency;

    // Scans
    std::vector<std::unique_ptr<Selection>> sel_tau32, sel_chi2;
    std::vector<std::unique_ptr<Hists>> hist_tau32, hist_chi2;
    int n_tau32, n_chi2;
    

  };

  BstarToTWAnalysisModule::BstarToTWAnalysisModule(Context & ctx) {
    dataset_name = ctx.get("dataset_version");

    // --- Kinematic Variables --- //
    double deltaPhi_min = M_PI/2;  // minimum delta phi between muon and top

    double top_fpt_max   = 0.8;    // maximum pt fraction of leading subjet
    double top_m_min     = 140.;   // minimum topjet mass
    double top_m_max     = 220.;   // maximum topjet mass
    double top_mpair_min = 50.;    // minimum pairwise mass of first three subjets
    double top_tau32_max = 0.58;   // maximum nsubjetiness tau_3/2

    CSVBTag::wp btag_wp = CSVBTag::WP_MEDIUM;

    // --- Object IDs --- //
    // TopJetId id_toptag = HOTVRTopTag(top_fpt_max, 0, FLT_MAX, top_mpair_min); // hotvr top tag without mass cut
    TopJetId id_toptag = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), Tau32Groomed(top_tau32_max), DeltaPhiCut(deltaPhi_min)); // hotvr top tag with tau_3/2 and delta phi
    JetId id_btag = CSVBTag(btag_wp);


    // --- Setup Additional Modules --- //
    // - Common Modules - //
    common.reset(new CommonModules());
    common->disable_jec();
    common->disable_jersmear();
    common->init(ctx);

    // - Bstar reconstruction - //
    primary_lep.reset(new PrimaryLepton(ctx));
    bstar_reco.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "BstarToTWReconstruction", id_toptag));
    bstar_disc.reset(new BstarToTWChi2Discriminator(ctx, "BstarToTWReconstruction"));
    if (dataset_name.find("BstarToTW") == 0)
      {
	BstarToTWgenprod.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
	bstar_matchdisc.reset(new BstarToTWMatchDiscriminator(ctx, "BstarToTWReconstruction"));
      }

    // - Scale Factors - //
    sf_muo_id.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1, "tightID"));
    sf_muo_iso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1, "iso"));
    sf_muo_trigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger"));
    // sf_btag.reset(new MCBTagScaleFactor(ctx, btag_wp));
   
    // --- Cleaner ---

    // --- Selections --- //
    // - TopTag - //
    sel_toptag.reset(new NTopJetSelection(1, 1, id_toptag));
    sel_tau32_cr.reset(new NTopJetSelection(2, 2, id_toptag));

    // - Leptons - //
    veto_ele.reset(new NElectronSelection(0, 0));
    sel_muo.reset(new NMuonSelection(1, 1));

    // - BTag - //
    sel_btag.reset(new NJetSelection(1, 1, id_btag));

    // --- Histograms --- //
    // - HOTVR JEC - //
    hist_hotvr_jec.reset(new HOTVRJetCorrectionHists(ctx, "HOTVR_JEC"));

    // - Cuts - //
    hist_presel.reset(new AndHists(ctx, "PreSel"));
    hist_sel_toptag.reset(new AndHists(ctx, "1TopTag", id_toptag));
    hist_veto_ele.reset(new AndHists(ctx, "ElectronVeto"));
    hist_sel_muo.reset(new AndHists(ctx, "1Muon"));
    hist_sel_btag.reset(new AndHists(ctx, "1BtagLoose"));

    // - Scale Factors - //
    hist_BTagMCEfficiency.reset(new BTagMCEfficiencyHists(ctx, "BTagMCEfficiency", btag_wp));
    hist_tau32MCEfficiency.reset(new AndHists(ctx, "Tau32_CR", id_toptag));

    // - Reconstruction - //
    hist_bstar_reco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWReco", "BstarToTWReconstruction", "Chi2"));
    hist_bstar_matchreco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWMatchedReco", "BstarToTWReconstruction", "Match"));
   
    // --- Scans --- //

    // - tau_3/2 Scan - //
    // n_tau32 = 11;
    // for (int i = 0; i < n_chi2; ++i)
    //   {
    // 	std::unique_ptr<Selection> sel;
    // 	sel.reset(new NTopJetSelectionSelection(1, 1, AndId<TopJet>(id_toptag, Tau32Groomed(0.65 - i*0.01))));
    // 	sel_chi2.push_back(std::move(sel));
    // 	std::unique_ptr<Hists> hist;
    // 	hist.reset(new BstarToTWHypothesisHists(ctx, "Tau32_Reco" + to_string(i), "BstarToTWReconstruction", "Chi2"));
    // 	hist_chi2.push_back(std::move(hist));
    //   }

    // - Chi2 Scan - //
    // n_chi2 = 6;
    // double chi2[n_chi2] = {1, 5, 10, 15, 20, 30};
    // for (int i = 0; i < n_chi2; ++i)
    //   {
    // 	std::unique_ptr<Selection> sel;
    // 	sel.reset(new Chi2Selection(ctx, "BstarToTWReconstruction", chi2[i]));
    // 	sel_chi2.push_back(std::move(sel));
    // 	std::unique_ptr<Hists> hist;
    // 	hist.reset(new BstarToTWHypothesisHists(ctx, "Chi2_Reco" + to_string(i), "BstarToTWReconstruction", "Chi2"));
    // 	hist_chi2.push_back(std::move(hist));
    //   }

  }

  bool BstarToTWAnalysisModule::process(Event & event) {

    // -- after PreSel -- //
    // apply scale factors for muons and trigger
    sf_muo_id->process(event);
    sf_muo_iso->process(event);
    sf_muo_trigger->process(event);
    common->process(event);
    
    hist_presel->fill(event);

    // if(!event.isRealData)
    //   {
    // 	hist_hotvr_jec->fill(event);
    //   }

    // -- TopTag Cut -- //    
    // Fill tau_3/2 controll region
    if(!sel_tau32_cr->passes(event)) hist_tau32MCEfficiency->fill(event);

    // Cut on TopTag
    if(!sel_toptag->passes(event)) return false;
    hist_sel_toptag->fill(event);

    // -- Bstar Reconstructinon -- //	
    primary_lep->process(event);
    bstar_reco->process(event);
    bstar_disc->process(event);
    hist_bstar_reco->fill(event);
    if (dataset_name.find("BstarToTW") == 0)
      {
    	BstarToTWgenprod->process(event);

    	bstar_matchdisc->process(event);
    	hist_bstar_matchreco->fill(event);
      }

    // -- Tau32 Scan -- //
    // for (int i = 0; i < n_tau32; ++i)
    //   {
    // 	if (sel_tau32.at(i)->passes(event))
    // 	  {
    // 	    hist_tau32.at(i)->fill(event);
    // 	  }
    //   }

    // -- Chi2 Scan -- //
    // for (int i = 0; i < n_chi2; ++i)
    //   {
    // 	if (sel_chi2.at(i)->passes(event))
    // 	  {
    // 	    hist_chi2.at(i)->fill(event);
    // 	  }
    //   }

    // -- Electron Veto -- //
    if (!veto_ele->passes(event)) return false;
    hist_veto_ele->fill(event);

    // -- Muon Cut -- //
    if (!sel_muo->passes(event)) return false;
    hist_sel_muo->fill(event);

    // -- BTag Cut -- //
    hist_BTagMCEfficiency->fill(event);
    if (!sel_btag->passes(event)) return false;
    // sf_btag->process(event);
    hist_sel_btag->fill(event);
   
    // -- Done -- //
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWAnalysisModule)

}
