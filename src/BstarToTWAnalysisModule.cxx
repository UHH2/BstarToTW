#include <iostream>
#include <memory>
#include <chrono>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/TopPtReweight.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/PrintingModules.h"

#include "UHH2/HOTVR/include/HOTVRHists.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"
#include "UHH2/HOTVR/include/HOTVRScaleFactor.h"
#include "UHH2/HOTVR/include/HadronicTop.h"

#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWModules.h"
#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/BstarToTW/include/BstarToTWSystematics.h"
#include "UHH2/BstarToTW/include/BstarToTWReconstruction.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisHists.h"
#include "UHH2/BstarToTW/include/BstarToTWPDFHists.h"
#include "UHH2/BstarToTW/include/BstarToTWHists.h"
#include "UHH2/BstarToTW/include/ElectroweakCorrections.h"
#include "UHH2/BstarToTW/include/BstarToTWSupplementPlotsModule.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWAnalysisModule: public AnalysisModule {
  public:
    
    explicit BstarToTWAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  
    bool is_ele, is_muo, do_for_cr, is_mc, do_scale_variation, do_pdf_variation, do_top_pt_reweight, do_supplement_plots;
    std::string dataset_name, prefiring_direction;
    // Event::Handle<float> h_tot_weight;

    // Modules for setting up event handles
    std::unique_ptr<AnalysisModule> primary_lep, hadronic_top, ht_calculator, bstar_gen;
    std::unique_ptr<AnalysisModule> output_module;
    std::unique_ptr<AnalysisModule> supplement_plots;
    std::unique_ptr<ObjectTagger> object_tagger;

    // Scale Factors & Systematics
    std::unique_ptr<AnalysisModule> sf_mc_lumiweight, sf_mc_pu_reweight;
    std::unique_ptr<AnalysisModule> sf_top_pt_reweight;
    std::unique_ptr<AnalysisModule> sf_lepton;
    std::unique_ptr<AnalysisModule> sf_btag;
    std::unique_ptr<AnalysisModule> sf_toptag;
    std::unique_ptr<AnalysisModule> sf_ewk;
    std::unique_ptr<AnalysisModule> sf_scale_variation;

    // Reconstruction Modules
    std::unique_ptr<AnalysisModule> reco_1toptag, disc_1toptag, disc_match;
    std::unique_ptr<Hists> hist_bstar_matchreco;

    // --- Selection + Hists
    // Pre Selection Histograms
    std::unique_ptr<Selection> sel_trigger, sel_badhcal;
    std::unique_ptr<AndHists> hist_trigger;
    std::unique_ptr<AndHists> hist_badhcal;
    std::unique_ptr<AndHists> hist_presel;
    // b tag MC efficiency histograms
    std::unique_ptr<Hists> hist_btag_mc_efficiency;
    // top tag
    std::unique_ptr<Selection> sel_1toptag;
    std::unique_ptr<AndHists> hist_1toptag;
    // reconstructed mass cut (to remove turn-on)
    std::unique_ptr<Selection> sel_recomass;

    // b tag categories
    std::unique_ptr<Selection> sel_0btag, sel_1btag, sel_2btag;
    std::unique_ptr<AndHists> hist_0btag1toptag, hist_1btag1toptag, hist_2btag1toptag;
    // 1+ btag region, for testing
    std::unique_ptr<Selection> sel_1plusbtag, sel_ptbal;
    std::unique_ptr<AndHists> hist_1plusbtag1toptag, hist_1plusbtag1toptag_tw, hist_1plusbtag1toptag_tt;
    std::unique_ptr<AndHists> hist_0btag1toptag_tw, hist_0btag1toptag_tt;
    // delta R (b, lepton)
    std::unique_ptr<Selection> sel_deltaR_bjetlep, sel_deltaR_bjetlep_cr; // delta R between primary lepton and b-jet/leading-jet
    std::unique_ptr<AndHists> hist_1btag1toptagdr;
    // chi2 
    std::unique_ptr<Selection> sel_chi2;
    std::unique_ptr<AndHists> hist_0btag1toptag20chi2, hist_1btag1toptag20chi2;
    std::unique_ptr<AndHists> hist_1btag1toptagdrcr, hist_1btag1toptagchi2cr;

  };

  BstarToTWAnalysisModule::BstarToTWAnalysisModule(Context & ctx) {
    // Set Flags
    // Year year;
    dataset_name = ctx.get("dataset_version");
    is_ele = ctx.get("analysis_channel") == "ELECTRON";
    is_muo = ctx.get("analysis_channel") == "MUON";
    do_for_cr = ctx.get("do_for_cr") == "true";
    is_mc = ctx.get("dataset_type") == "MC";
    prefiring_direction = ctx.get("Systematic_Prefiring", "nominal");
    do_supplement_plots = ctx.get("do_supplement_plots","false") == "true";

    // genprinter.reset(new GenParticlesPrinter(ctx));

    // --- Selection Parameters
    double top_fpt_max   = 0.8;    // maximum pt fraction of leading subjet
    double top_m_min     = 140.;   // minimum topjet mass
    double top_m_max     = 220.;   // maximum topjet mass
    double top_mpair_min = 50.;    // minimum pairwise mass of first three subjets
    double top_tau32_max = 0.56;   // maximum nsubjetiness tau_3/2
    double chi2_max      = 20.0;   // maximum chi2 of the reconstructed bstar hypothesis
    double deltaR_blep_min = 2.0;  // minimuim deltaR between lepton and (b)jet
    BTag::algo btag_algo = BTag::DEEPJET; // b tag algortihm
    BTag::wp btag_wp = BTag::WP_MEDIUM;   // b tag working point
    // b tag id
    JetId id_btag = BTag(btag_algo, btag_wp);
    // top tag ids
    TopJetId id_toptag = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), Tau32Groomed(top_tau32_max)); // standard top tag ID
    // TopJetId id_toptag_without_mass = AndId<TopJet>(HOTVRTopTag(top_fpt_max, 0.0, std::numeric_limits<double>::infinity(), top_mpair_min), Tau32Groomed(top_tau32_max)); // top tag ID without mass window cut

    // --- Reconstruction Modules
    string name_tw_reco = "tW_reco";
    // string name_discriminator = "Chi2"; 
    string name_discriminator = "closest_nu";
   reco_1toptag.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, name_tw_reco, id_toptag));
    disc_1toptag.reset(new BstarToTWChi2Discriminator(ctx, name_tw_reco));
    if (dataset_name.find("BstarToTW") == 0) // matching discriminator only works on signal
      {
    	bstar_gen.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    	disc_match.reset(new BstarToTWMatchDiscriminator(ctx, name_tw_reco));
      }
    hist_bstar_matchreco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWMatchedReco", name_tw_reco, "Match"));

    // --- Modules for setting up event handles
    primary_lep.reset(new PrimaryLepton(ctx));
    hadronic_top.reset(new HadronicTop(ctx));
    ht_calculator.reset(new HTCalculator(ctx));
    object_tagger.reset(new ObjectTagger());
    object_tagger->set_btag_loose_id(BTag(btag_algo, BTag::WP_LOOSE));
    object_tagger->set_btag_medium_id(BTag(btag_algo, BTag::WP_MEDIUM));
    object_tagger->set_btag_tight_id(BTag(btag_algo, BTag::WP_TIGHT));
    object_tagger->set_toptag_id(id_toptag);
    object_tagger->init(ctx);
    output_module.reset(new BstarToTWOutputModule(ctx, name_tw_reco));
    if (do_supplement_plots) supplement_plots.reset(new BstarToTWSupplementPlotsModule(ctx));

    // --- Scale Factors
    // parameters
    double par_top_pt_a = 0.0615; // top pt reweighting parameter a
    double par_top_pt_b = -0.000; // top pt reweighting parameter b
    do_scale_variation = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down" ||  ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");
    do_pdf_variation = (ctx.get("b_PDFVariation") == "true");
    // modules
    sf_top_pt_reweight.reset(new TopPtReweighting(ctx, par_top_pt_a, par_top_pt_b, ctx.get("Systematic_TopPt_a", "nominal"), ctx.get("Systematic_TopPt_b", "nominal"), "", "")); // https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting
    if (is_mc)
      {
	sf_mc_lumiweight.reset(new MCLumiWeight(ctx));
	sf_mc_pu_reweight.reset(new MCPileupReweight(ctx, ctx.get("Systematic_PU", "central")));
      }
    sf_lepton.reset(new LeptonScaleFactors(ctx));
    sf_btag.reset(new MCBTagScaleFactor(ctx, btag_algo, btag_wp, "jets", ctx.get("Systematic_BTag","central"), "comb", "incl", "MCBtagEfficiencies"));
    sf_toptag.reset(new HOTVRScaleFactor(ctx, id_toptag, ctx.get("Systematic_TopTag", "nominal"), "HadronicTop", "TopTagSF", "HOTVRTopTagSFs"));
    sf_ewk.reset(new ElectroweakCorrections(ctx));
    if (do_scale_variation) sf_scale_variation.reset(new MCScaleVariation(ctx));

    // ---- Selection + Hists
    // Trigger
    sel_trigger.reset(new BstarToTWTriggerSelection(ctx));
    hist_trigger.reset(new AndHists(ctx, "Trigger")); 
    // bad HCAL
    sel_badhcal.reset(new BadHCALSelection(ctx));
    hist_badhcal.reset(new AndHists(ctx, "BadHCAL"));
    hist_badhcal->add_hist(new HOTVRHists(ctx, "BadHCAL_HOTVR_tagged", id_toptag));
    // presel hists
    hist_presel.reset(new AndHists(ctx, "PreSel"));
    hist_presel->add_hist(new HOTVRPerformanceHists(ctx, "PreSel_Performance_HOTVR"));
    hist_btag_mc_efficiency.reset(new BTagMCEfficiencyHists(ctx,"BTagMCEfficiency", id_btag));    
    // toptag
    sel_1toptag.reset(new NTopJetSelection(1, 1, id_toptag));
    sel_recomass.reset(new RecoMassSelection(ctx, 500., name_tw_reco));    
    hist_1toptag.reset(new AndHists(ctx, "1toptag"));
    hist_1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "1toptag_reco", name_tw_reco, name_discriminator));
    hist_1toptag->add_hist(new HOTVRHists(ctx, "1toptag_HOTVR_tagged", id_toptag));
    // -- btag
    sel_0btag.reset(new NJetSelection(0, 0, id_btag));

    // - 0 btag
    hist_0btag1toptag.reset(new AndHists(ctx, "0btag1toptag"));
    hist_0btag1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "0btag1toptag_reco", name_tw_reco, name_discriminator));
    hist_0btag1toptag->add_hist(new HOTVRHists(ctx, "0btag1toptag_HOTVR_tagged", id_toptag));

    // - 1+ btag
    sel_1plusbtag.reset(new NJetSelection(1, -1, id_btag));
    hist_1plusbtag1toptag.reset(new AndHists(ctx, "1plusbtag1toptag"));
    hist_1plusbtag1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "1plusbtag1toptag_reco", name_tw_reco, name_discriminator));
    hist_1plusbtag1toptag->add_hist(new HOTVRMatchingHists(ctx, "1plusbtag1toptag_topmatch"));
    hist_1plusbtag1toptag->add_hist(new HOTVRHists(ctx, "1plusbtag1toptag_HOTVR_tagged", id_toptag));
    // PtBalance
    sel_ptbal.reset(new Chi2Selection(ctx, name_tw_reco, 0.05, "deltaPt_W"));
    hist_1plusbtag1toptag_tw.reset(new AndHists(ctx, "1plusbtag1toptag_tw"));
    hist_1plusbtag1toptag_tw->add_hist(new BstarToTWHypothesisHists(ctx, "1plusbtag1toptag_tw_reco", name_tw_reco, name_discriminator));
    hist_1plusbtag1toptag_tw->add_hist(new HOTVRMatchingHists(ctx, "1plusbtag1toptag_tw_topmatch"));
    hist_1plusbtag1toptag_tw->add_hist(new HOTVRHists(ctx, "1plusbtag1toptag_tw_HOTVR_tagged", id_toptag));
    hist_1plusbtag1toptag_tt.reset(new AndHists(ctx, "1plusbtag1toptag_tt"));
    hist_1plusbtag1toptag_tt->add_hist(new BstarToTWHypothesisHists(ctx, "1plusbtag1toptag_tt_reco", name_tw_reco, name_discriminator));
    hist_1plusbtag1toptag_tt->add_hist(new HOTVRMatchingHists(ctx, "1plusbtag1toptag_tt_topmatch"));
    hist_1plusbtag1toptag_tt->add_hist(new HOTVRHists(ctx, "1plusbtag1toptag_tt_HOTVR_tagged", id_toptag));
    hist_0btag1toptag_tw.reset(new AndHists(ctx, "0btag1toptag_tw"));
    hist_0btag1toptag_tw->add_hist(new BstarToTWHypothesisHists(ctx, "0btag1toptag_tw_reco", name_tw_reco, name_discriminator));
    hist_0btag1toptag_tw->add_hist(new HOTVRMatchingHists(ctx, "0btag1toptag_tw_topmatch"));
    hist_0btag1toptag_tw->add_hist(new HOTVRHists(ctx, "0btag1toptag_tw_HOTVR_tagged", id_toptag));
    hist_0btag1toptag_tt.reset(new AndHists(ctx, "0btag1toptag_tt"));
    hist_0btag1toptag_tt->add_hist(new BstarToTWHypothesisHists(ctx, "0btag1toptag_tt_reco", name_tw_reco, name_discriminator));
    hist_0btag1toptag_tt->add_hist(new HOTVRMatchingHists(ctx, "0btag1toptag_tt_topmatch"));
    hist_0btag1toptag_tt->add_hist(new HOTVRHists(ctx, "0btag1toptag_tt_HOTVR_tagged", id_toptag));

    // - 1 btag
    sel_1btag.reset(new NJetSelection(1, 1, id_btag));
    hist_1btag1toptag.reset(new AndHists(ctx, "1btag1toptag"));
    hist_1btag1toptag->add_hist(new HOTVRMatchingHists(ctx, "1btag1toptag_topmatch"));
    hist_1btag1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "1btag1toptag_reco", name_tw_reco, name_discriminator));
    hist_1btag1toptag->add_hist(new HOTVRHists(ctx, "1btag1toptag_HOTVR_tagged", id_toptag));
    // deltaR btag/leading jet
    sel_deltaR_bjetlep.reset(new LeadingJetDeltaRSelection(ctx, deltaR_blep_min, ctx.get_handle<vector<Jet> >("btag_medium")));
    sel_deltaR_bjetlep_cr.reset(new LeadingJetDeltaRSelection(ctx, 1.0, ctx.get_handle<vector<Jet> >("btag_medium")));
    hist_1btag1toptagdr.reset(new AndHists(ctx, "1btag1toptagdr"));
    hist_1btag1toptagdr->add_hist(new BstarToTWHypothesisHists(ctx, "1btag1toptagdr_reco", name_tw_reco, name_discriminator));
    hist_1btag1toptagdr->add_hist(new HOTVRHists(ctx, "1btag1toptagdr_HOTVR_tagged", id_toptag));
    hist_1btag1toptagdr->add_hist(new HOTVRMatchingHists(ctx, "1btag1toptagdr_topmatch"));
    // chi2
    sel_chi2.reset(new Chi2Selection(ctx, name_tw_reco, chi2_max, name_discriminator));
    hist_1btag1toptag20chi2.reset(new AndHists(ctx, "1btag1toptag20chi2"));
    hist_1btag1toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "1btag1toptag20chi2_reco", name_tw_reco, name_discriminator));
    hist_1btag1toptag20chi2->add_hist(new HOTVRHists(ctx, "1btag1toptag20chi2_HOTVR_tagged", id_toptag));
    hist_1btag1toptag20chi2->add_hist(new HOTVRMatchingHists(ctx, "1btag1toptag20chi2_topmatch"));
    if (do_pdf_variation) hist_1btag1toptag20chi2->add_hist(new BstarToTWPDFHists(ctx, "1btag1toptag20chi2_PDF_variations", true, do_pdf_variation));
    // validation regions
    hist_1btag1toptagchi2cr.reset(new AndHists(ctx, "1btag1toptagchi2cr"));
    hist_1btag1toptagchi2cr->add_hist(new BstarToTWHypothesisHists(ctx, "1btag1toptagchi2cr_reco", name_tw_reco, name_discriminator));
    hist_1btag1toptagchi2cr->add_hist(new HOTVRHists(ctx, "1btag1toptagchi2cr_HOTVR_tagged", id_toptag));
    hist_1btag1toptagchi2cr->add_hist(new HOTVRMatchingHists(ctx, "1btag1toptagchi2cr_topmatch"));
    hist_1btag1toptagdrcr.reset(new AndHists(ctx, "1btag1toptagdrcr"));
    hist_1btag1toptagdrcr->add_hist(new BstarToTWHypothesisHists(ctx, "1btag1toptagdrcr_reco", name_tw_reco, name_discriminator));
    hist_1btag1toptagdrcr->add_hist(new HOTVRHists(ctx, "1btag1toptagdrcr_HOTVR_tagged", id_toptag));
    hist_1btag1toptagdrcr->add_hist(new HOTVRMatchingHists(ctx, "1btag1toptagdrcr_topmatch"));

    hist_0btag1toptag20chi2.reset(new AndHists(ctx, "0btag1toptag20chi2"));
    hist_0btag1toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "0btag1toptag20chi2_reco", name_tw_reco, name_discriminator));
    hist_0btag1toptag20chi2->add_hist(new HOTVRHists(ctx, "0btag1toptag20chi2_HOTVR_tagged", id_toptag));
    hist_0btag1toptag20chi2->add_hist(new HOTVRMatchingHists(ctx, "0btag1toptag20chi2_topmatch"));

    // - 2 btag
    sel_2btag.reset(new NJetSelection(2, -1, id_btag));
    hist_2btag1toptag.reset(new AndHists(ctx, "2btag1toptag"));
    hist_2btag1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "2btag1toptag_reco", name_tw_reco, name_discriminator));
    hist_2btag1toptag->add_hist(new HOTVRMatchingHists(ctx, "2btag1toptag_topmatch"));
    hist_2btag1toptag->add_hist(new HOTVRHists(ctx, "2btag1toptag_HOTVR_tagged", id_toptag));
    if (do_pdf_variation) hist_2btag1toptag->add_hist(new BstarToTWPDFHists(ctx, "2btag1toptag_PDF_variations", true, do_pdf_variation));
  }

  bool BstarToTWAnalysisModule::process(Event & event) {
    // setup all the event handles
    primary_lep->process(event);
    hadronic_top->process(event);
    ht_calculator->process(event);
    object_tagger->process(event);
    if (dataset_name.find("BstarToTW") == 0)
      {
	bstar_gen->process(event);
      }
    
    // Scale Factors
    if(is_mc)
      {
	sf_mc_lumiweight->process(event);
	sf_mc_pu_reweight->process(event);
      }
    sf_lepton->process(event);
    sf_top_pt_reweight->process(event);
    sf_ewk->process(event);
    // prefiring weights
    if ( (!event.isRealData) ) 
      {
	if (prefiring_direction == "nominal") event.weight *= event.prefiringWeight;
	else if (prefiring_direction == "up") event.weight *= event.prefiringWeightUp;
	else if (prefiring_direction == "down") event.weight *= event.prefiringWeightDown;
      }
    if(do_scale_variation)
      {
	sf_scale_variation->process(event);
      }

   // Fill PreSel Histograms with all scale factors
    if (!sel_trigger->passes(event)) return false;
    hist_trigger->fill(event);

    // bad HCAL selection
    if (!sel_badhcal->passes(event)) {
      if (event.isRealData) return false;
      else event.weight *= (1 - 0.64844705699);
    }
    hist_badhcal->fill(event);

    hist_btag_mc_efficiency->fill(event);
    sf_btag->process(event);
    sf_toptag->process(event);

    hist_presel->fill(event);
    if (do_supplement_plots) supplement_plots->process(event);
    
    // - 1 toptag - //
    if (sel_1toptag->passes(event))
      {
	reco_1toptag->process(event);
	disc_1toptag->process(event);
	if (dataset_name.find("BstarToTW") == 0)
	  {
	    disc_match->process(event);
	    hist_bstar_matchreco->fill(event);
	  }	    
	if (!sel_recomass->passes(event)) return false;
	
	// - 1+ btag region for testign
	if (sel_1plusbtag->passes(event) && sel_deltaR_bjetlep->passes(event))
	  {
	    hist_1plusbtag1toptag->fill(event);
	    // split into tt and tw regions
	    if (sel_ptbal->passes(event))
	      hist_1plusbtag1toptag_tw->fill(event);
	    else
	      hist_1plusbtag1toptag_tt->fill(event);
	  }
	// - 0 btag regions - //
	if (sel_0btag->passes(event))
	  {
	    hist_0btag1toptag->fill(event);
	    // split into tt and tw regions
	    if (sel_ptbal->passes(event))
	      hist_0btag1toptag_tw->fill(event);
	    else
	      hist_0btag1toptag_tt->fill(event);
	    // - Chi2
	    if (sel_chi2->passes(event))
	      {
	    	hist_0btag1toptag20chi2->fill(event);
	      }
	    if (do_for_cr)
	      output_module->process(event);
	    return do_for_cr;
	  }
	else if (sel_1btag->passes(event)) 
	  {
	    hist_1btag1toptag->fill(event);
	    // - deltaR (b,lep)
	    if (sel_deltaR_bjetlep->passes(event)) 
	      {
		hist_1btag1toptagdr->fill(event);
		// - Chi2
		if (sel_chi2->passes(event))
		  {
		  hist_1btag1toptag20chi2->fill(event);
		if (!do_for_cr)
		  output_module->process(event);	
		return !do_for_cr;
		  }
		else
		  hist_1btag1toptagchi2cr->fill(event);
	      }
	    else if (sel_deltaR_bjetlep_cr->passes(event) && sel_chi2->passes(event))
	      {
		hist_1btag1toptagdrcr->fill(event);
	      }
	  }
	else if (sel_2btag->passes(event)) 
	  {
	    hist_2btag1toptag->fill(event);
	    if (!do_for_cr)
	      output_module->process(event);	
	    return !do_for_cr;
	  }
      }
    
    return false;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWAnalysisModule)

}
