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
#include "UHH2/common/include/TopPtReweight.h"

#include "UHH2/HOTVR/include/HOTVRHists.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"
#include "UHH2/HOTVR/include/HOTVRScaleFactor.h"
#include "UHH2/HOTVR/include/HOTVRPileUpUncertainty.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWModules.h"
#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/BstarToTW/include/BstarToTWReconstruction.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisHists.h"
#include "UHH2/BstarToTW/include/BstarToTWPDFHists.h"
#include "UHH2/BstarToTW/include/AndHists.h"


using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWAnalysisModule: public AnalysisModule {
  public:
    
    explicit BstarToTWAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  
    std::string dataset_name;

    bool is_ele, is_muo;

    // --- Additional Modules --- //
    // Common Modules
    std::unique_ptr<CommonModules> common;
    // Bstar Reconsturction
    std::unique_ptr<AnalysisModule> BstarToTWgenprod;
    std::unique_ptr<AnalysisModule> primary_lep;
    std::unique_ptr<AnalysisModule> bstar_reco, bstar_disc, bstar_matchdisc;
    std::unique_ptr<AnalysisModule> w_reco, w_disc;
    std::unique_ptr<AnalysisModule> bstar_reco_anti, bstar_disc_anti;
    std::unique_ptr<AndHists> hist_jec_validation;
    std::unique_ptr<Hists> hist_bstar_reco, hist_bstar_reco_anti, hist_bstar_matchreco, hist_w_reco;

    // Scale factors & Uncertainties
    std::unique_ptr<AnalysisModule> sf_muo_id, sf_muo_iso, sf_muo_trigger, sf_muo_trk;
    std::unique_ptr<AnalysisModule> sf_ele_id, sf_ele_rec, sf_ele_trigger;
    std::unique_ptr<AnalysisModule> sf_btag_tight;
    std::unique_ptr<AnalysisModule> sf_top_pt_reweight;
    bool do_top_pt_reweight;
    std::unique_ptr<AnalysisModule> sf_toptag;
    std::unique_ptr<AnalysisModule> scale_variation;
    std::unique_ptr<AnalysisModule> L1_variation;
    bool do_scale_variation, do_pdf_variations;
    std::unique_ptr<Hists> pdf_variations;

    // --- Selections and Histogramms --- //
    std::unique_ptr<AndHists> hist_presel; 

    // - toptag - //
    std::unique_ptr<Selection> sel_1toptag, sel_1antitoptag;
    std::unique_ptr<AndHists> hist_sel_1toptag, hist_sel_1antitoptag;
    // ttbar CR //
    std::unique_ptr<Selection> sel_2btag_tight;
    std::unique_ptr<Selection> sel_antibtag_medium;
    std::unique_ptr<AndHists> hist_sel_2btag, hist_sel_not2btag, hist_toptag_efficiency, hist_toptag_sf_closure;

    // - chi2 - //
    std::unique_ptr<Selection> sel_20chi2, sel_20chi2_anti;
    std::unique_ptr<AndHists> hist_sel_20chi2, hist_sel_20chi2_anti;

    // - additional hists - //
    std::unique_ptr<Hists> hist_BTagMCEfficiency;

   

  };

  BstarToTWAnalysisModule::BstarToTWAnalysisModule(Context & ctx) {


    dataset_name = ctx.get("dataset_version");

    is_ele = ctx.get("analysis_channel") == "ELECTRON";
    is_muo = ctx.get("analysis_channel") == "MUON";

    // --- Selection Variables --- //
    double deltaPhi_min = M_PI/2;  // minimum delta phi between muon and top

    double top_fpt_max   = 0.8;    // maximum pt fraction of leading subjet
    double top_m_min     = 140.;   // minimum topjet mass
    double top_m_max     = 220.;   // maximum topjet mass
    double top_mpair_min = 50.;    // minimum pairwise mass of first three subjets
    double top_tau32_max = 0.56;   // maximum nsubjetiness tau_3/2
    double top_pt_min    = 200.0;
    double top_eta_max   = 2.5;

    double chi2_max      = 20;     // maximum chi2 of the reconstructed bstar hypothesis

    // -- Object IDs -- //
    TopJetId id_topjet =  AndId<TopJet>(PtEtaCut(top_pt_min, top_eta_max), DeltaPhiCut(ctx, deltaPhi_min));
    TopJetId id_toptag_without_tau32 = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), DeltaPhiCut(ctx, deltaPhi_min));
    TopJetId id_toptag = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), Tau32Groomed(top_tau32_max), DeltaPhiCut(ctx, deltaPhi_min)); // hotvr top tag with tau_3/2 and delta phi
    TopJetId id_antitoptag = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), AntiTau32Groomed(top_tau32_max), DeltaPhiCut(ctx, deltaPhi_min)); // hotvr top tag with tau_3/2 and delta phi

    CSVBTag::wp btag_wp_tight = CSVBTag::WP_TIGHT;
    CSVBTag::wp btag_wp_medium = CSVBTag::WP_MEDIUM;
    JetId id_btag_tight  = CSVBTag(btag_wp_tight);
    JetId id_btag_medium  = CSVBTag(btag_wp_medium);


    // --- Setup Additional Modules --- //

    // - Common Modules - //
    string sys_pu = ctx.get("Systematic_PU");
    common.reset(new CommonModules());
    common->disable_jec();
    common->disable_jersmear();
    common->init(ctx, sys_pu);

    // -- Scale Factors & Uncertainties -- //
    do_scale_variation = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down") || (ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");
    scale_variation.reset(new MCScaleVariation(ctx));

    do_pdf_variations = ctx.get("b_PDFUncertainties") == "true";
    pdf_variations.reset(new BstarToTWPDFHists(ctx, "PDF_variations", true, do_pdf_variations));

    do_top_pt_reweight = ctx.get("b_TopPtReweight") == "true";
    sf_top_pt_reweight.reset(new TopPtReweight(ctx, 0.0615, -0.0005, "", "", do_top_pt_reweight, 1.0)); // https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting

    // - Muon Scale Factors - //
    if (is_muo)
      {
	string sys_muo_id      = ctx.get("Systematic_MuonID");
	string sys_muo_trigger = ctx.get("Systematic_MuonTrigger");
	string sys_muo_iso     = ctx.get("Systematic_MuonIso");
	string sys_muo_trk     = ctx.get("Systematic_MuonTrk");
	string path_sf_muo_id      = "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root";
	string path_sf_muo_iso     = "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root";
	string path_sf_muo_trigger = "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root";
	string path_sf_muo_trk     = "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/Tracking_EfficienciesAndSF_BCDEFGH.root";

	sf_muo_id.reset(new MCMuonScaleFactor(ctx, path_sf_muo_id, "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1, "tightID", true, sys_muo_id));
	sf_muo_iso.reset(new MCMuonScaleFactor(ctx, path_sf_muo_iso, "TightISO_TightID_pt_eta", 1, "iso", true, sys_muo_iso));
	sf_muo_trigger.reset(new MCMuonScaleFactor(ctx, path_sf_muo_trigger, "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger", true, sys_muo_trigger));
	sf_muo_trk.reset(new MCMuonTrkScaleFactor(ctx, path_sf_muo_trk, 1, "track", sys_muo_trk, "muons"));
      }

    // - Electron Scale Factors - //
    if (is_ele)
      {
	string sys_ele_id      = ctx.get("Systematic_ElectronID");
	string sys_ele_rec     = ctx.get("Systematic_ElectronRec");
	string sys_ele_trigger = ctx.get("Systematic_ElectronTrigger");
	string path_sf_ele_id      = "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Tight_ID.root";
	string path_sf_ele_rec     = "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root";
	string path_sf_ele_trigger = "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/BstarToTW/data/ElectronEfficiencies.root";
	sf_ele_id.reset(new MCElecScaleFactor(ctx, path_sf_ele_id, 1, "TightID", sys_ele_id));
	sf_ele_rec.reset(new MCElecScaleFactor(ctx, path_sf_ele_rec, 1, "RecEff", sys_ele_rec));
	sf_ele_trigger.reset(new ElectronTriggerWeights(ctx, path_sf_ele_trigger, sys_ele_trigger));
      }

    string sys_btag_tight   = ctx.get("Systematic_BTag");
    sf_btag_tight.reset(new MCBTagScaleFactor(ctx, btag_wp_tight, "jets", sys_btag_tight, "mujets", "incl", "MCBtagEfficienciesTight"));

    string sys_toptag       = ctx.get("Systematic_TopTag");
    string path_sf_toptag = "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/BstarToTW/data/TopTagMCEfficiency.root";
    sf_toptag.reset(new HOTVRScaleFactor(ctx, "BstarToTW", id_toptag, path_sf_toptag, sys_toptag));

    string sys_L1           = ctx.get("Systematic_L1");
    string path_L1_variation = "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/BstarToTW/data/HOTVR_L1Uncertainty.root";
    L1_variation.reset(new HOTVRPileUpUncertainty(ctx, path_L1_variation, sys_L1 ));

    // -- Bstar reconstruction -- //
    primary_lep.reset(new PrimaryLepton(ctx));
    bstar_reco.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "BstarToTWReconstruction", id_toptag));
    bstar_reco_anti.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "BstarToTWReconstruction_Anti", id_antitoptag));
    w_reco.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "WReconstruction", id_topjet));
    bstar_disc.reset(new BstarToTWChi2Discriminator(ctx, "BstarToTWReconstruction"));
    bstar_disc_anti.reset(new BstarToTWChi2Discriminator(ctx, "BstarToTWReconstruction_Anti"));
    w_disc.reset(new BstarToTWChi2Discriminator(ctx, "WReconstruction"));
    if (dataset_name.find("BstarToTW") == 0)
      {
    	BstarToTWgenprod.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    	bstar_matchdisc.reset(new BstarToTWMatchDiscriminator(ctx, "BstarToTWReconstruction"));
      }
    hist_bstar_reco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWReco", "BstarToTWReconstruction", "Chi2"));
    hist_bstar_reco_anti.reset(new BstarToTWHypothesisHists(ctx, "AntiBstarToTWReco", "BstarToTWReconstruction_Anti", "Chi2"));
    hist_w_reco.reset(new BstarToTWHypothesisHists(ctx, "WReco", "WReconstruction", "Chi2"));
    hist_bstar_matchreco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWMatchedReco", "BstarToTWReconstruction", "Match"));

    // --- Selections and Histogramms --- //
    hist_presel.reset(new AndHists(ctx, "PreSel"));
    hist_presel->add_hist(new HOTVRHists(ctx, "PreSel_HOTVRtagged", id_toptag));

    // - TopTag - //
    sel_1toptag.reset(new NTopJetSelection(1, 1, id_toptag));
    hist_sel_1toptag.reset(new AndHists(ctx, "1TopTag", id_toptag));

    sel_1antitoptag.reset(new NTopJetSelection(1, 1, id_antitoptag));
    hist_sel_1antitoptag.reset(new AndHists(ctx, "1AntiTopTag", id_antitoptag));
    // TopTag CR //

    sel_2btag_tight.reset(new NJetSelection(2, -1, id_btag_tight));
    sel_antibtag_medium.reset(new NJetSelection(0, 0, id_btag_tight));
    hist_jec_validation.reset(new AndHists(ctx, "JEC_Validation"));
    hist_sel_2btag.reset(new AndHists(ctx, "2btag"));
    hist_sel_not2btag.reset(new AndHists(ctx, "less2btag"));
    // - Chi2 - //
    sel_20chi2.reset(new Chi2Selection(ctx, "BstarToTWReconstruction", chi2_max));
    sel_20chi2_anti.reset(new Chi2Selection(ctx, "BstarToTWReconstruction_Anti", chi2_max));
    hist_sel_20chi2.reset(new AndHists(ctx, "20Chi2", id_toptag));
    hist_sel_20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "20Chi2_reco", "BstarToTWReconstruction", "Chi2"));
    hist_sel_20chi2_anti.reset(new AndHists(ctx, "20Chi2_Anti", id_antitoptag));
    hist_sel_20chi2_anti->add_hist(new BstarToTWHypothesisHists(ctx, "20Chi2_reco_Anti", "BstarToTWReconstruction_Anti", "Chi2"));

    // - Scale Factor Hists - //
    hist_BTagMCEfficiency.reset(new BTagMCEfficiencyHists(ctx, "BTagMCEfficiency", btag_wp_tight));    
    hist_toptag_efficiency.reset(new AndHists(ctx, "TopTagEfficiency"));   
    hist_toptag_efficiency->add_hist(new HOTVRHists(ctx, "TopTagEfficiency_HOTVRtagged", id_toptag));
    hist_toptag_sf_closure.reset(new AndHists(ctx, "TopTagSFClosure"));
    hist_toptag_sf_closure->add_hist(new HOTVRHists(ctx, "TopTagSFClosure_HOTVRtagged", id_toptag));

  }

  bool BstarToTWAnalysisModule::process(Event & event) {

    // -- after PreSel -- //
    common->process(event);
    primary_lep->process(event);

    // apply scale factors for muons and trigger
    if (is_muo)
      {
	sf_muo_id->process(event);
	sf_muo_iso->process(event);
	sf_muo_trigger->process(event);
	sf_muo_trk->process(event);
      }
    if (is_ele)
      {
	sf_ele_id->process(event);
	sf_ele_rec->process(event);
	// sf_ele_trigger->process(event);
      }

    hist_presel->fill(event);

    L1_variation->process(event);
    if(do_scale_variation) scale_variation->process(event);
    if(do_top_pt_reweight) sf_top_pt_reweight->process(event);   

    // -- TopTag Cut -- //    
    // Fill ttbar controll region
    hist_BTagMCEfficiency->fill(event);
    sf_btag_tight->process(event);
    if(sel_2btag_tight->passes(event)) 
      {
	hist_sel_2btag->fill(event);
	hist_toptag_efficiency->fill(event);

	sf_toptag->process(event);
	hist_toptag_sf_closure->fill(event);
	return false;
      }
    hist_sel_not2btag->fill(event);


    // Cut on TopTag
    if(sel_1toptag->passes(event))
      {
	sf_toptag->process(event);
	hist_sel_1toptag->fill(event);

	// // -- Bstar Reconstructinon -- //	
	bstar_reco->process(event);
	bstar_disc->process(event);
	hist_bstar_reco->fill(event);
	if (dataset_name.find("BstarToTW") == 0)
	  {
	    BstarToTWgenprod->process(event);
	    bstar_matchdisc->process(event);
	    hist_bstar_matchreco->fill(event);
	  }

	// -- Chi2 Cut -- //
	if (!sel_20chi2->passes(event)) return false;
	if (do_pdf_variations) pdf_variations->fill(event);
	hist_sel_20chi2->fill(event);
      }
    else if (sel_1antitoptag->passes(event))
      {
	hist_sel_1antitoptag->fill(event);
	bstar_reco_anti->process(event);
	bstar_disc_anti->process(event);
	hist_bstar_reco_anti->fill(event);
	if (!sel_20chi2_anti->passes(event)) return false;
	hist_sel_20chi2_anti->fill(event);
	return false;
      }

    // -- Done -- //
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWAnalysisModule)

}
