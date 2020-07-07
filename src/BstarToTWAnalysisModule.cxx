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

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWAnalysisModule: public AnalysisModule {
  public:
    
    explicit BstarToTWAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  
    bool is_ele, is_muo, do_for_cr, is_mc, do_scale_variation, do_pdf_variation, do_top_pt_reweight;
    std::string dataset_name, prefiring_direction;
    // Event::Handle<float> h_tot_weight;

    // Modules for setting up event handles
    std::unique_ptr<AnalysisModule> primary_lep, hadronic_top, ht_calculator, bstar_gen;
    std::unique_ptr<AnalysisModule> output_module, genprinter;
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
    std::unique_ptr<AnalysisModule> reco_0toptag, disc_0toptag;
    std::unique_ptr<AnalysisModule> reco_toplep, disc_toplep;

    // Selection + Hists
    std::unique_ptr<AndHists> hist_trigger;
    std::unique_ptr<Selection> sel_trigger;
    std::unique_ptr<AndHists> hist_presel;
    std::unique_ptr<Selection> sel_badhcal;
    std::unique_ptr<AndHists> hist_badhcal;
    std::unique_ptr<Hists> hist_btag_mc_efficiency;

    std::unique_ptr<Selection> sel_deltaR_bjetlep; // delta R between primary lepton and b-jet/leading-jet
    std::unique_ptr<AndHists> hist_sel_1btagdr;
    
    std::unique_ptr<Selection> sel_1btag, sel_2btag; // btag
    std::unique_ptr<AndHists> hist_sel_1btag, hist_sel_2btag; 
    std::unique_ptr<Selection> veto_btag; // btag veto
    std::unique_ptr<AndHists> hist_sel_0btag;

    std::unique_ptr<Selection> sel_1toptag; // toptag
    std::unique_ptr<AndHists> hist_sel_1btag1toptag, hist_sel_0btag1toptag, hist_sel_2btag1toptag;
    std::unique_ptr<Selection> veto_toptag; //toptag veto
    std::unique_ptr<AndHists> hist_sel_1btag0toptag, hist_sel_0btag0toptag, hist_sel_2btag0toptag;

    std::unique_ptr<Selection> sel_subjet_btag;
    std::unique_ptr<AndHists> hist_sel_subjet_btag;
    
    std::unique_ptr<Selection> sel_1toptag20chi2, sel_0toptag20chi2; // chi2 cut (for 0 and 1 toptag reco)
    std::unique_ptr<Selection> sel_1toprecomass, sel_0toprecomass; // reco mass cut (for 0 and 1 toptag reco)
    std::unique_ptr<AndHists> hist_sel_1btag1toptag20chi2, hist_sel_1btag0toptag20chi2, hist_sel_0btag1toptag20chi2, hist_sel_0btag0toptag20chi2;

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

    // genprinter.reset(new GenParticlesPrinter(ctx));

    // Selection Parameters
    double top_fpt_max   = 0.8;    // maximum pt fraction of leading subjet
    double top_m_min     = 140.;   // minimum topjet mass
    double top_m_max     = 220.;   // maximum topjet mass
    double top_mpair_min = 50.;    // minimum pairwise mass of first three subjets
    double top_tau32_max = 0.56;   // maximum nsubjetiness tau_3/2
    double chi2_max      = 20.0;   // maximum chi2 of the reconstructed bstar hypothesis
    double deltaR_blep_min = 2.0;  // minimuim deltaR between lepton and (b)jet
    BTag::algo btag_algo = BTag::DEEPJET;
    BTag::wp btag_wp = BTag::WP_MEDIUM;

    // Toptags:
    TopJetId id_toptag = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), Tau32Groomed(top_tau32_max));
    // if (!is_hotvr)
    //   id_toptag = AndId<TopJet>(Type2TopTag(105., 220., Type2TopTag::MassType::groomed), Tau32(0.50));


    TopJetId id_toptag_veto = VetoId<TopJet>(id_toptag);
    // btag:
    JetId id_btag = BTag(btag_algo, btag_wp);

    // Modules for setting up event handles
    primary_lep.reset(new PrimaryLepton(ctx));
    hadronic_top.reset(new HadronicTop(ctx));
    ht_calculator.reset(new HTCalculator(ctx));
    object_tagger.reset(new ObjectTagger());
    object_tagger->set_btag_loose_id(BTag(btag_algo, BTag::WP_LOOSE));
    object_tagger->set_btag_medium_id(BTag(btag_algo, BTag::WP_MEDIUM));
    object_tagger->set_btag_tight_id(BTag(btag_algo, BTag::WP_TIGHT));
    object_tagger->set_toptag_id(id_toptag);
    object_tagger->init(ctx);

    output_module.reset(new BstarToTWOutputModule(ctx,"1TopTagReconstruction"));

    // - Trigger -
    sel_trigger.reset(new BstarToTWTriggerSelection(ctx));
    hist_trigger.reset(new AndHists(ctx, "Trigger")); 


    // - Scale Factors -

    do_scale_variation = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down" ||  ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");
    do_pdf_variation = (ctx.get("b_PDFVariation") == "true");

    double par_top_pt_a = 0.0615;
    double par_top_pt_b = -0.0005;
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

    // Reconstruction Modules
    reco_1toptag.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "1TopTagReconstruction", id_toptag));
    disc_1toptag.reset(new BstarToTWChi2Discriminator(ctx, "1TopTagReconstruction"));

    reco_toplep.reset(new LeptonicTopReconstruction(ctx, "1TopTagReconstruction", "btag_medium", "TopLepReco"));
    disc_toplep.reset(new LeptonicTopChi2Discriminator(ctx, "TopLepReco"));

    reco_0toptag.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "0TopTagReconstruction", id_toptag_veto));
    disc_0toptag.reset(new BstarToTWChi2Discriminator(ctx, "0TopTagReconstruction"));
    if (dataset_name.find("BstarToTW") == 0) // matching discriminator only works on signal
      {
    	bstar_gen.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    	disc_match.reset(new BstarToTWMatchDiscriminator(ctx, "1TopTagReconstruction"));
      }
    hist_bstar_matchreco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWMatchedReco", "1TopTagReconstruction", "Match"));

    // Selection + Hists
    hist_presel.reset(new AndHists(ctx, "PreSel"));
    hist_presel->add_hist(new HOTVRPerformanceHists(ctx, "PreSel_Performance_HOTVR"));
    hist_btag_mc_efficiency.reset(new BTagMCEfficiencyHists(ctx,"BTagMCEfficiency", id_btag));    
    // - bad HCAL
    sel_badhcal.reset(new BadHCALSelection(ctx));
    hist_badhcal.reset(new AndHists(ctx, "BadHCAL"));
    hist_badhcal->add_hist(new HOTVRHists(ctx, "BadHCAL_HOTVR_tagged", id_toptag));
    
    // btag
    sel_1btag.reset(new NJetSelection(1, 1, id_btag));
    hist_sel_1btag.reset(new AndHists(ctx, "1btag"));
    hist_sel_1btag->add_hist(new HOTVRMatchingHists(ctx, "1btag_topmatch"));

    sel_2btag.reset(new NJetSelection(2, -1, id_btag));
    hist_sel_2btag.reset(new AndHists(ctx, "2btag"));
    hist_sel_2btag->add_hist(new HOTVRMatchingHists(ctx, "2btag_topmatch"));

    veto_btag.reset(new NJetSelection(0, 0, id_btag));
    hist_sel_0btag.reset(new AndHists(ctx, "0btag"));

    // deltaR btag/leading jet
    sel_deltaR_bjetlep.reset(new LeadingJetDeltaRSelection(ctx, deltaR_blep_min, ctx.get_handle<vector<Jet> >("btag_medium")));
    hist_sel_1btagdr.reset(new AndHists(ctx, "1btagdr"));
    // toptag
    sel_1toptag.reset(new NTopJetSelection(1, 1, id_toptag));
    hist_sel_1btag1toptag.reset(new AndHists(ctx, "1btag1toptag"));
    hist_sel_1btag1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "1btag1toptag_reco", "1TopTagReconstruction", "Chi2"));
    hist_sel_1btag1toptag->add_hist(new LeptonicTopHypothesisHists(ctx, "1btag1toptag_toplep", "TopLepReco", "Chi2"));
    hist_sel_1btag1toptag->add_hist(new HOTVRMatchingHists(ctx, "1btag1toptag_topmatch"));
    hist_sel_1btag1toptag->add_hist(new HOTVRHists(ctx, "1btag1toptag_HOTVR_tagged", id_toptag));

    hist_sel_2btag1toptag.reset(new AndHists(ctx, "2btag1toptag"));
    hist_sel_2btag1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "2btag1toptag_reco", "1TopTagReconstruction", "Chi2"));
    hist_sel_2btag1toptag->add_hist(new LeptonicTopHypothesisHists(ctx, "2btag1toptag_toplep", "TopLepReco", "Chi2"));
    hist_sel_2btag1toptag->add_hist(new HOTVRMatchingHists(ctx, "2btag1toptag_topmatch"));
    hist_sel_2btag1toptag->add_hist(new HOTVRHists(ctx, "2btag1toptag_HOTVR_tagged", id_toptag));
    if (do_pdf_variation) hist_sel_2btag1toptag->add_hist(new BstarToTWPDFHists(ctx, "2btag1toptag_PDF_variations", true, do_pdf_variation));

    hist_sel_0btag1toptag.reset(new AndHists(ctx, "0btag1toptag"));
    hist_sel_0btag1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "0btag1toptag_reco", "1TopTagReconstruction", "Chi2"));
    hist_sel_0btag1toptag->add_hist(new HOTVRHists(ctx, "0btag1toptag_HOTVR_tagged", id_toptag));
    veto_toptag.reset(new NTopJetSelection(0, 0, id_toptag));
    hist_sel_1btag0toptag.reset(new AndHists(ctx, "1btag0toptag"));
    hist_sel_1btag0toptag->add_hist(new BstarToTWHypothesisHists(ctx, "1btag0toptag_reco", "0TopTagReconstruction", "Chi2"));
    hist_sel_2btag0toptag.reset(new AndHists(ctx, "2btag0toptag"));
    hist_sel_2btag0toptag->add_hist(new BstarToTWHypothesisHists(ctx, "2btag0toptag_reco", "0TopTagReconstruction", "Chi2"));
    hist_sel_0btag0toptag.reset(new AndHists(ctx, "0btag0toptag"));
    hist_sel_0btag0toptag->add_hist(new BstarToTWHypothesisHists(ctx, "0btag0toptag_reco", "0TopTagReconstruction", "Chi2"));
    // chi2
    sel_1toprecomass.reset(new RecoMassSelection(ctx, 500., "1TopTagReconstruction"));
    sel_1toptag20chi2.reset(new Chi2Selection(ctx, "1TopTagReconstruction", chi2_max, "Chi2"));
    hist_sel_1btag1toptag20chi2.reset(new AndHists(ctx, "1btag1toptag20chi2"));
    hist_sel_1btag1toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "1btag1toptag20chi2_reco", "1TopTagReconstruction", "Chi2"));
    hist_sel_1btag1toptag20chi2->add_hist(new HOTVRHists(ctx, "1btag1toptag20chi2_HOTVR_tagged", id_toptag));
    hist_sel_1btag1toptag20chi2->add_hist(new HOTVRMatchingHists(ctx, "1btag1toptag20chi2_topmatch"));
    if (do_pdf_variation) hist_sel_1btag1toptag20chi2->add_hist(new BstarToTWPDFHists(ctx, "1btag1toptag20chi2_PDF_variations", true, do_pdf_variation));

    hist_sel_0btag1toptag20chi2.reset(new AndHists(ctx, "0btag1toptag20chi2"));
    hist_sel_0btag1toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "0btag1toptag20chi2_reco", "1TopTagReconstruction", "Chi2"));
    sel_0toprecomass.reset(new RecoMassSelection(ctx, 500., "0TopTagReconstruction"));
    sel_0toptag20chi2.reset(new Chi2Selection(ctx, "0TopTagReconstruction", chi2_max, "Chi2"));
    hist_sel_1btag0toptag20chi2.reset(new AndHists(ctx, "1btag0toptag20chi2"));
    hist_sel_1btag0toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "1btag0toptag20chi2_reco", "0TopTagReconstruction", "Chi2"));
    hist_sel_0btag0toptag20chi2.reset(new AndHists(ctx, "0btag0toptag20chi2"));
    hist_sel_0btag0toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "0btag0toptag20chi2_reco", "0TopTagReconstruction", "Chi2"));

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
    
    // genprinter->process(event);

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

    hist_presel->fill(event);
    if (!sel_trigger->passes(event)) return false;
    hist_trigger->fill(event);

   // Fill PreSel Histograms with all scale factors
    


    // bad HCAL selection
    if (!sel_badhcal->passes(event)) return false;
    hist_badhcal->fill(event);

    hist_btag_mc_efficiency->fill(event);
    sf_btag->process(event);
    sf_toptag->process(event);

    // - 0 btag regions - //
    if (veto_btag->passes(event))
      {
	hist_sel_0btag->fill(event);
	// - 1 toptag - //
	if (sel_1toptag->passes(event))
	  {
	    reco_1toptag->process(event);
	    disc_1toptag->process(event);
	    if (!sel_1toprecomass->passes(event)) return false;
	    hist_sel_0btag1toptag->fill(event);
	    if (sel_1toptag20chi2->passes(event))
	      {
		hist_sel_0btag1toptag20chi2->fill(event);
	      }
	    if (do_for_cr)
	      output_module->process(event);
	    return do_for_cr;
	  }
	else if (veto_toptag->passes(event))
	  {
	    reco_0toptag->process(event);
	    disc_0toptag->process(event);
	    if (!sel_0toprecomass->passes(event)) return false;
	    hist_sel_0btag0toptag->fill(event);
	    if (sel_0toptag20chi2->passes(event))
	      {	      
		hist_sel_0btag0toptag20chi2->fill(event);
	      }
	  }
	return false;
      }

    else if (sel_1btag->passes(event)) 
      {
	hist_sel_1btag->fill(event);
	if (!sel_deltaR_bjetlep->passes(event)) return false;
	hist_sel_1btagdr->fill(event);
	// - Signal Region - //
	if(sel_1toptag->passes(event))
	  {	  	    
	    // -- Bstar Reconstructinon -- //	
	    reco_1toptag->process(event);
	    disc_1toptag->process(event);
	    reco_toplep->process(event);
	    disc_toplep->process(event);
	    if (!sel_1toprecomass->passes(event)) return false;
	    hist_sel_1btag1toptag->fill(event);
	    if (dataset_name.find("BstarToTW") == 0)
	      {
		disc_match->process(event);
		hist_bstar_matchreco->fill(event);
	      }	    
	    // -- Chi2 Cut -- //
	    if (!sel_1toptag20chi2->passes(event)) return false;
	    hist_sel_1btag1toptag20chi2->fill(event);
	    if (!do_for_cr)
	      output_module->process(event);	
	    return !do_for_cr;
	  }
	// - Non-top control Region - //
	else if (veto_toptag->passes(event))
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
	return false;
      }
    else if (sel_2btag->passes(event)) 
      {
	hist_sel_2btag->fill(event);
	// - TTbar Control Region - //
	if(sel_1toptag->passes(event))
	  {	    
	    // -- Bstar Reconstructinon -- //	
	    reco_1toptag->process(event);
	    disc_1toptag->process(event);
	    reco_toplep->process(event);
	    disc_toplep->process(event);
	    if (!sel_1toprecomass->passes(event)) return false;
	    hist_sel_2btag1toptag->fill(event);
	    if (dataset_name.find("BstarToTW") == 0)
	      {
		disc_match->process(event);
		hist_bstar_matchreco->fill(event);
	      }
	  }
	else if (veto_toptag->passes(event))
	  {
	    reco_0toptag->process(event);
	    disc_0toptag->process(event);
	    if (!sel_0toprecomass->passes(event)) return false;
	    hist_sel_2btag0toptag->fill(event);
	  }
	return false;
      }
    
    return false;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWAnalysisModule)

}
