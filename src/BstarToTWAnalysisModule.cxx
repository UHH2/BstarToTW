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

#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWModules.h"
#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/BstarToTW/include/BstarToTWReconstruction.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisHists.h"
#include "UHH2/BstarToTW/include/BstarToTWPDFHists.h"
#include "UHH2/BstarToTW/include/BstarToTWHists.h"


using namespace std;
using namespace uhh2;

namespace uhh2 {
  /* To do:
   * - top reconstruction performance (using PUPPI)
   * - background estimation
   * - Check top mistag (had W instead of had t) in signal
   * - Try using BDT instead of chi2 method after reconstruction
   */

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
    std::unique_ptr<AnalysisModule> cl_topjet;
    std::unique_ptr<AnalysisModule> background_reweight;
    
    // Bstar Reconsturction
    std::unique_ptr<AnalysisModule> BstarToTWgenprod;
    std::unique_ptr<AnalysisModule> primary_lep;
    std::unique_ptr<AnalysisModule> reco_1toptag, disc_1toptag, bstar_matchdisc;
    std::unique_ptr<AnalysisModule> reco_0toptag, disc_0toptag;
    std::unique_ptr<AndHists> hist_jec_validation;
    std::unique_ptr<Hists> hist_bstar_matchreco;

    // --- Scale factors & Uncertainties --- //
    std::unique_ptr<AnalysisModule> sf_muo_id, sf_muo_iso, sf_muo_trigger, sf_muo_trk;
    std::unique_ptr<AnalysisModule> sf_ele_id, sf_ele_rec, sf_ele_trigger;
    std::unique_ptr<AnalysisModule> sf_btag;
    std::unique_ptr<AnalysisModule> sf_top_pt_reweight;
    bool do_top_pt_reweight;
    std::unique_ptr<AnalysisModule> sf_toptag;
    std::unique_ptr<AnalysisModule> scale_variation;
    // std::unique_ptr<AnalysisModule> L1_variation;
    bool do_scale_variation, do_pdf_variations;
    std::unique_ptr<Hists> pdf_variations;

    // --- Selections and Histogramms --- //
    std::unique_ptr<AndHists> hist_presel; 
    /* Regions
     * Signal:    1 b-tag, 1 top-tag
     * ttbar:     2 b-tag, 1 top-tag
     * non-ttbar: 0 b-tag, 0 top-tag
     * CR_A:      0 b-tag, 1 top-tag
     * CR_B:      1 b-tag, 0 top-tag
     * CR_C       2 b-tag, 0 top-tag
     */

    // - Region Hists - //
    std::unique_ptr<AndHists> hist_sel_1btag, hist_sel_2btag, hist_sel_0btag;

    std::unique_ptr<AndHists> hist_sel_1btag1toptag, hist_sel_1btag0toptag;
    std::unique_ptr<AndHists> hist_sel_2btag1toptag, hist_sel_2btag0toptag;
    std::unique_ptr<AndHists> hist_sel_0btag1toptag, hist_sel_0btag0toptag;

    std::unique_ptr<AndHists> hist_sel_1btag1toptag20chi2, hist_sel_1btag0toptag20chi2;
    std::unique_ptr<AndHists> hist_sel_2btag1toptag20chi2, hist_sel_2btag0toptag20chi2;
    std::unique_ptr<AndHists> hist_sel_0btag1toptag20chi2, hist_sel_0btag0toptag20chi2;

    std::unique_ptr<Hists> hist_1btag0toptagreweighted, hist_0btag1toptagreweighted;
    std::unique_ptr<Hists> hist_1btag1toptag20chi2reweighted, hist_2btag1toptag20chi2reweighted;

    std::unique_ptr<AndHists> hist_sel_0btag0toptagreweighted;
    // - btag - //
    std::unique_ptr<Selection> sel_1btag, sel_2btag, veto_btag;

    // - toptag - //
    std::unique_ptr<Selection> sel_1toptag, veto_toptag, sel_1antitoptag, sel_1top;

    // - chi2 - //
    std::unique_ptr<Selection> sel_1toptag20chi2, sel_0toptag20chi2;
    std::unique_ptr<Selection> sel_1toprecomass, sel_0toprecomass;
    std::unique_ptr<AndHists> hist_sel_20chi2;

    // - additional hists - //
    std::unique_ptr<Hists> hist_BTagMCEfficiency;

  };

  BstarToTWAnalysisModule::BstarToTWAnalysisModule(Context & ctx) {

    // --- Flags --- //
    dataset_name = ctx.get("dataset_version");
    is_ele = ctx.get("analysis_channel") == "ELECTRON";
    is_muo = ctx.get("analysis_channel") == "MUON";


    // --- Selection Variables --- //
    double deltaPhi_min  = M_PI/2; // minimum delta phi between muon and top
    double top_fpt_max   = 0.8;    // maximum pt fraction of leading subjet
    double top_m_min     = 140.;   // minimum topjet mass
    double top_m_max     = 220.;   // maximum topjet mass
    double top_mpair_min = 50.;    // minimum pairwise mass of first three subjets
    double top_tau32_max = 0.56;   // maximum nsubjetiness tau_3/2
    // double top_pt_min    = 200.0;  // minimum topjet pt
    // double top_eta_max   = 2.5;    // maximum topjet eta
    double chi2_max      = 20;     // maximum chi2 of the reconstructed bstar hypothesis
    CSVBTag::wp btag_wp  = CSVBTag::WP_LOOSE;  // b-tag workingpoint


    // --- Object IDs --- //
    TopJetId id_topjet = DeltaPhiCut(ctx, deltaPhi_min);
    // Toptags:
    TopJetId id_toptag = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), Tau32Groomed(top_tau32_max));
    TopJetId id_toptag_inv_tau32 = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), VetoId<TopJet>(Tau32Groomed(top_tau32_max)));
    TopJetId id_toptag_veto = VetoId<TopJet>(id_toptag);

    // btag:
    JetId id_btag  = CSVBTag(btag_wp);


    // --- Setup Additional Modules --- //
    primary_lep.reset(new PrimaryLepton(ctx));
    cl_topjet.reset(new TopJetCleaner(ctx, id_topjet));
    // - Common Modules - //
    string sys_pu = ctx.get("Systematic_PU");
    common.reset(new CommonModules());
    common->disable_jec();
    common->disable_jersmear();
    common->switch_jetPtSorter();
    common->init(ctx, sys_pu);

    // --- Scale Factors & Uncertainties --- //
    do_scale_variation = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down") || (ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");
    scale_variation.reset(new MCScaleVariation(ctx));

    do_pdf_variations = ctx.get("b_PDFUncertainties") == "true";
    pdf_variations.reset(new BstarToTWPDFHists(ctx, "PDF_variations", true, do_pdf_variations));

    do_top_pt_reweight = ctx.get("b_TopPtReweight") == "true";
    // sf_top_pt_reweight.reset(new TopPtReweight(ctx, 0.159, -0.00141, "", "", do_top_pt_reweight, 1.0)); // https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting 8TeV recommendations
    sf_top_pt_reweight.reset(new TopPtReweight(ctx, 0.0615, -0.0005, "", "", do_top_pt_reweight, 1.0)); // https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting

    // - Muon Scale Factors - //
    if (is_muo) // for muon channel only
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
    if (is_ele) // for electron channel only
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

    string sys_btag   = ctx.get("Systematic_BTag");
    sf_btag.reset(new MCBTagScaleFactor(ctx, btag_wp, "jets", sys_btag, "mujets", "incl", "MCBtagEfficienciesLoose"));

    string sys_toptag       = ctx.get("Systematic_TopTag");
    string path_sf_toptag = "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/BstarToTW/data/TopTagMCEfficiency.root";
    sf_toptag.reset(new HOTVRScaleFactor(ctx, "BstarToTW", id_toptag, path_sf_toptag, sys_toptag));

    // string sys_L1           = ctx.get("Systematic_L1");
    // string path_L1_variation = "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/BstarToTW/data/HOTVR_L1Uncertainty.root";
    // L1_variation.reset(new HOTVRPileUpUncertainty(ctx, path_L1_variation, sys_L1 ));

    // --- Bstar reconstruction --- //
    reco_1toptag.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "1TopTagReconstruction", id_toptag));
    disc_1toptag.reset(new BstarToTWChi2Discriminator(ctx, "1TopTagReconstruction"));
    reco_0toptag.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "0TopTagReconstruction", id_toptag_veto));
    disc_0toptag.reset(new BstarToTWChi2Discriminator(ctx, "0TopTagReconstruction"));
    if (dataset_name.find("BstarToTW") == 0)
      {
    	BstarToTWgenprod.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    	bstar_matchdisc.reset(new BstarToTWMatchDiscriminator(ctx, "1TopTagReconstruction"));
      }
    // hist_bstar_reco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWReco", "BstarToTWReconstruction", "Chi2"));
    // hist_reco_cr.reset(new BstarToTWHypothesisHists(ctx, "CRReco", "CR_Reconstruction", "Chi2"));
    hist_bstar_matchreco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWMatchedReco", "1TopTagReconstruction", "Match"));

    string sys_background_reweight = "";//ctx.get("Systematic_BG");
    string path_background_reweight = ctx.get("BackgroundShapeNorm");
    background_reweight.reset(new BackgroundShapeNormWeights(ctx,path_background_reweight,"0TopTagReconstruction","Chi2",sys_background_reweight));
    // --- Selections and Histogramms --- //
    sel_1top.reset(new NTopJetSelection(1,-1, id_topjet));
    hist_presel.reset(new AndHists(ctx, "PreSel"));
    hist_presel->add_hist(new HOTVRPerformanceHists(ctx, "PreSel_HOTVRPerformance"));


    // - bTag - //
    sel_1btag.reset(new NJetSelection(1, 1, id_btag));
    sel_2btag.reset(new NJetSelection(2,-1, id_btag));
    veto_btag.reset(new NJetSelection(0, 0, id_btag));

    // - TopTag - //
    veto_toptag.reset(new NTopJetSelection(0, 0, id_toptag));
    sel_1toptag.reset(new NTopJetSelection(1, 1, id_toptag));
    sel_1antitoptag.reset(new NTopJetSelection(1, -1, id_toptag_veto));

    // - Regions - //
    hist_sel_1btag.reset(new AndHists(ctx, "1btag"));
    hist_sel_2btag.reset(new AndHists(ctx, "2btag"));
    hist_sel_0btag.reset(new AndHists(ctx, "0btag"));

    hist_sel_1btag1toptag.reset(new AndHists(ctx, "1btag1toptag"));
    hist_sel_1btag1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "1btag1toptag_reco", "1TopTagReconstruction", "Chi2"));
    hist_sel_1btag0toptag.reset(new AndHists(ctx, "1btag0toptag"));
    hist_sel_1btag0toptag->add_hist(new BstarToTWHypothesisHists(ctx, "1btag0toptag_reco", "0TopTagReconstruction", "Chi2"));
    hist_sel_2btag1toptag.reset(new AndHists(ctx, "2btag1toptag"));
    hist_sel_2btag1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "2btag1toptag_reco", "1TopTagReconstruction", "Chi2"));
    hist_sel_2btag0toptag.reset(new AndHists(ctx, "2btag0toptag"));
    hist_sel_2btag0toptag->add_hist(new BstarToTWHypothesisHists(ctx, "2btag0toptag_reco", "0TopTagReconstruction", "Chi2"));
    hist_sel_0btag1toptag.reset(new AndHists(ctx, "0btag1toptag"));
    hist_sel_0btag1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "0btag1toptag_reco", "1TopTagReconstruction", "Chi2"));
    hist_sel_0btag0toptag.reset(new AndHists(ctx, "0btag0toptag"));
    hist_sel_0btag0toptag->add_hist(new BstarToTWHypothesisHists(ctx, "0btag0toptag_reco", "0TopTagReconstruction", "Chi2"));

    hist_sel_1btag1toptag20chi2.reset(new AndHists(ctx, "1btag1toptag20chi2"));
    hist_sel_1btag1toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "1btag1toptag20chi2_reco", "1TopTagReconstruction", "Chi2"));
    hist_sel_1btag0toptag20chi2.reset(new AndHists(ctx, "1btag0toptag20chi2"));
    hist_sel_1btag0toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "1btag0toptag20chi2_reco", "0TopTagReconstruction", "Chi2"));
    hist_sel_2btag1toptag20chi2.reset(new AndHists(ctx, "2btag1toptag20chi2"));
    hist_sel_2btag1toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "2btag1toptag20chi2_reco", "1TopTagReconstruction", "Chi2"));
    hist_sel_2btag0toptag20chi2.reset(new AndHists(ctx, "2btag0toptag20chi2"));
    hist_sel_2btag0toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "2btag0toptag20chi2_reco", "0TopTagReconstruction", "Chi2"));
    hist_sel_0btag1toptag20chi2.reset(new AndHists(ctx, "0btag1toptag20chi2"));
    hist_sel_0btag1toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "0btag1toptag20chi2_reco", "1TopTagReconstruction", "Chi2"));
    hist_sel_0btag0toptag20chi2.reset(new AndHists(ctx, "0btag0toptag20chi2"));
    hist_sel_0btag0toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "0btag0toptag20chi2_reco", "0TopTagReconstruction", "Chi2"));

    // - Reweighted Hists for Background estimation -//
    // hist_sel_0btag0toptagreweighted->add_hist(new BstarToTWHypothesisHists(ctx, "0btag0toptagreweighted_reco", "0TopTagReconstruction", "Chi2"));
    
    string path_1b_background = ctx.get("BackgroundShapeNorm_1b0t");
    string path_1t_background = ctx.get("BackgroundShapeNorm_0b1t");
    string path_1b_signal = ctx.get("BackgroundShapeNorm_1b1t");
    string path_2b_signal = ctx.get("BackgroundShapeNorm_2b1t");
    
    hist_1btag0toptagreweighted.reset(new BstarToTWBackgroundHists(ctx, "1btag0toptagreweighted_reco", "0TopTagReconstruction", path_1b_background));
    hist_0btag1toptagreweighted.reset(new BstarToTWBackgroundHists(ctx, "0btag1toptagreweighted_reco", "0TopTagReconstruction", path_1t_background));
    hist_1btag1toptag20chi2reweighted.reset(new BstarToTWBackgroundHists(ctx, "1btag1toptagreweighted_reco", "0TopTagReconstruction", path_1b_signal));
    hist_2btag1toptag20chi2reweighted.reset(new BstarToTWBackgroundHists(ctx, "2btag1toptagreweighted_reco", "0TopTagReconstruction", path_2b_signal));
    
    
    // - Chi2 - //
    sel_1toptag20chi2.reset(new Chi2Selection(ctx, "1TopTagReconstruction", chi2_max));
    sel_0toptag20chi2.reset(new Chi2Selection(ctx, "0TopTagReconstruction", chi2_max));

    // - RecoMass - //
    sel_1toprecomass.reset(new RecoMassSelection(ctx, 500., "1TopTagReconstruction"));
    sel_0toprecomass.reset(new RecoMassSelection(ctx, 500., "0TopTagReconstruction"));

    // - Scale Factor Hists - //
    hist_BTagMCEfficiency.reset(new BTagMCEfficiencyHists(ctx, "BTagMCEfficiency", btag_wp));

  }

  bool BstarToTWAnalysisModule::process(Event & event) {
    // -- after PreSel -- //
    common->process(event);
    primary_lep->process(event);

    if (dataset_name.find("BstarToTW") == 0)
      {
	BstarToTWgenprod->process(event);
      }

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
	sf_ele_trigger->process(event);
      }

    // move this to presel!
    cl_topjet->process(event);
    if(!sel_1top->passes(event))
      return false;
    // L1_variation->process(event);
    if(do_scale_variation) scale_variation->process(event);
    if(do_top_pt_reweight) sf_top_pt_reweight->process(event);   
    hist_presel->fill(event);

    // -- Fill regions -- //
    hist_BTagMCEfficiency->fill(event);

    // - 2 btag regions - //
    if (sel_2btag->passes(event))
      {
	sf_btag->process(event);
	hist_sel_2btag->fill(event);
	if (sel_1toptag->passes(event))
	  {
	    reco_1toptag->process(event);
	    disc_1toptag->process(event);
	    if (!sel_1toprecomass->passes(event)) return false;
	    hist_sel_2btag1toptag->fill(event);
	    if (sel_1toptag20chi2->passes(event))
	      hist_sel_2btag1toptag20chi2->fill(event);
	  }
	else if (sel_1antitoptag->passes(event) && veto_toptag->passes(event))
	  {
	    reco_0toptag->process(event);
	    disc_0toptag->process(event);
	    if (!sel_0toprecomass->passes(event)) return false;
	    hist_sel_2btag0toptag->fill(event);
	    if (sel_0toptag20chi2->passes(event))
	      hist_sel_2btag0toptag20chi2->fill(event);
	  }
	return false;
      }
 
    // - 0 btag regions - //
    else if (veto_btag->passes(event))
      {
	hist_sel_0btag->fill(event);
	if (sel_1toptag->passes(event))
	  {
	    reco_1toptag->process(event);
	    disc_1toptag->process(event);
	    if (!sel_1toprecomass->passes(event)) return false;
	    hist_sel_0btag1toptag->fill(event);
	    if (sel_1toptag20chi2->passes(event))
	      hist_sel_0btag1toptag20chi2->fill(event);
	  }
	else if (sel_1antitoptag->passes(event) && veto_toptag->passes(event))
	  {
	    reco_0toptag->process(event);
	    disc_0toptag->process(event);
	    if (!sel_0toprecomass->passes(event)) return false;
	    hist_sel_0btag0toptag->fill(event);
	    //background_reweight->process(event);
	    //hist_sel_0btag0toptagreweighted->fill(event);
	    hist_1btag0toptagreweighted->fill(event);
	    hist_0btag1toptagreweighted->fill(event);
	    hist_1btag1toptag20chi2reweighted->fill(event);
	    hist_2btag1toptag20chi2reweighted->fill(event);
	    if (sel_0toptag20chi2->passes(event))
	      {	      
		hist_sel_0btag0toptag20chi2->fill(event);
	      }
	  }
	return false;
      }

    
    else if (sel_1btag->passes(event)) 
      {
	sf_btag->process(event);
	hist_sel_1btag->fill(event);
	// - Signal Region - //
	if(sel_1toptag->passes(event))
	  {	    
	    // -- Bstar Reconstructinon -- //	
	    reco_1toptag->process(event);
	    disc_1toptag->process(event);
	    if (!sel_1toprecomass->passes(event)) return false;
	    hist_sel_1btag1toptag->fill(event);
	    
	    if (dataset_name.find("BstarToTW") == 0)
	      {
		bstar_matchdisc->process(event);
		hist_bstar_matchreco->fill(event);
	      }

	    // -- Chi2 Cut -- //
	    if (!sel_1toptag20chi2->passes(event)) return false;
	    if (do_pdf_variations) pdf_variations->fill(event);
	    hist_sel_1btag1toptag20chi2->fill(event);
	  }
	// - Non-top control Region - //
	else if (sel_1antitoptag->passes(event) && veto_toptag->passes(event))
	  {
	    reco_0toptag->process(event);
	    disc_0toptag->process(event);
	    if (!sel_0toprecomass->passes(event)) return false;
	    hist_sel_1btag0toptag->fill(event);
	    if (sel_0toptag20chi2->passes(event))
	      {
		hist_sel_1btag0toptag20chi2->fill(event);
	      }
	  }
      }

    // -- Done -- //
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWAnalysisModule)

}
