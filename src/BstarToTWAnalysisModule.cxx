#include <iostream>
#include <memory>

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

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWAnalysisModule: public AnalysisModule {
  public:
    
    explicit BstarToTWAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  
    bool is_ele, is_muo, do_for_cr, is_mc;
    std::string dataset_name;
    // Event::Handle<float> h_tot_weight;

    // Modules for setting up event handles
    std::unique_ptr<AnalysisModule> primary_lep, hadronic_top, ht_calculator, bstar_gen;
    std::unique_ptr<AnalysisModule> output_module;
    std::unique_ptr<ObjectTagger> object_tagger;

    // Scale Factors
    std::unique_ptr<AnalysisModule> sf_mc_lumiweight, sf_mc_pu_reweight, sf_lepton;
    // TO DO!

    // Reconstruction Modules
    std::unique_ptr<AnalysisModule> reco_1toptag, disc_1toptag, disc_match;
    std::unique_ptr<Hists> hist_bstar_matchreco;
    std::unique_ptr<AnalysisModule> reco_0toptag, disc_0toptag;

    // Selection + Hists
    std::unique_ptr<AndHists> hist_presel;
    std::unique_ptr<Hists> hist_btag_mc_efficiency; // hists for calculating btag scale factors

    std::unique_ptr<Selection> sel_deltaPhi_metlep; // delta phi between primary lepton and MET
    std::unique_ptr<AndHists> hist_sel_deltaPhi_metlep;
    std::unique_ptr<Selection> sel_deltaR_jetlep; // delta R between primary lepton and ak4-jets
    std::unique_ptr<AndHists> hist_sel_deltaR_jetlep;
    std::unique_ptr<Selection> sel_deltaR_leadingjetlep, sel_deltaR_leadingbjetlep; // delta R between primary lepton and b-jet/leading-jet
    std::unique_ptr<AndHists> hist_sel_1btagdr, hist_sel_2btagdr, hist_sel_0btagdr;
    
    std::unique_ptr<Selection> sel_1btag; // btag
    std::unique_ptr<AndHists> hist_sel_1btag; 
    std::unique_ptr<Selection> veto_btag; // btag veto
    std::unique_ptr<AndHists> hist_sel_0btag;

    std::unique_ptr<Selection> sel_1toptag; // toptag
    std::unique_ptr<AndHists> hist_sel_1btag1toptag, hist_sel_0btag1toptag;
    std::unique_ptr<Selection> veto_toptag; //toptag veto
    std::unique_ptr<AndHists> hist_sel_1btag0toptag, hist_sel_0btag0toptag;

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

    // h_tot_weight = ctx.declare_event_output<float>("weight_total");

    // Modules for setting up event handles
    primary_lep.reset(new PrimaryLepton(ctx));
    hadronic_top.reset(new HadronicTop(ctx));
    ht_calculator.reset(new HTCalculator(ctx));
    object_tagger.reset(new ObjectTagger());
    object_tagger->init(ctx);
    // if (do_for_cr)
    //   output_module.reset(new BstarToTWOutputModule(ctx,"0TopTagReconstruction"));
    // else
    //   output_module.reset(new BstarToTWOutputModule(ctx,"1TopTagReconstruction"));

    // Scale Factor
    if (is_mc)
      {
	sf_mc_lumiweight.reset(new MCLumiWeight(ctx));
	std::string syst_pu = ctx.get("Systematic_PU");
	sf_mc_pu_reweight.reset(new MCPileupReweight(ctx, syst_pu));
      }
    // TO DO!
    //include years!
    if (is_muo)
      {
	sf_lepton.reset(new MuonScaleFactors2018(ctx));
      }
    else if (is_ele)
      {
	// sf_lepton.reset(MuonScaleFactors2018(ctx);
      }
    // Selection Parameters
    // double deltaPhi_min  = M_PI/2; // minimum delta phi between muon and top
    double top_fpt_max   = 0.8;    // maximum pt fraction of leading subjet
    double top_m_min     = 140.;   // minimum topjet mass
    double top_m_max     = 220.;   // maximum topjet mass
    double top_mpair_min = 50.;    // minimum pairwise mass of first three subjets
    double top_tau32_max = 0.56;   // maximum nsubjetiness tau_3/2
    double chi2_max      = 30.0;   // maximum chi2 of the reconstructed bstar hypothesis
    double deltaR_blep_min = 2.0;  // minimuim deltaR between lepton and (b)jet
    double deltaPhi_metlep_max = M_PI/2; // maximum deltaPhi between missign ET and lepton
    DeepCSVBTag::wp btag_wp  = DeepCSVBTag::WP_LOOSE;  // b-tag workingpoint
    // Toptags:
    TopJetId id_toptag = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), Tau32Groomed(top_tau32_max));
    TopJetId id_toptag_veto = VetoId<TopJet>(id_toptag);
    // btag:
    JetId id_btag  = DeepCSVBTag(btag_wp);
    
    // Reconstruction Modules
    reco_1toptag.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "1TopTagReconstruction", id_toptag));
    disc_1toptag.reset(new BstarToTWChi2Discriminator(ctx, "1TopTagReconstruction"));
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
    std::unique_ptr<BTagMCEfficiencyHists<DeepCSVBTag>> hists_btag_mc_efficiency;
    hists_btag_mc_efficiency.reset(new BTagMCEfficiencyHists<DeepCSVBTag>(ctx,"BTagLoose", btag_wp));
    // hist_BTagMCEfficiency.reset(new BTagMCEfficiencyHists(ctx, "BTagMCEfficiency", btag_wp));
    // deltaPhi met-lep
    sel_deltaPhi_metlep.reset(new METDeltaPhiSelection(ctx, deltaPhi_metlep_max));
    hist_sel_deltaPhi_metlep.reset(new AndHists(ctx, "deltaPhiMETLep"));
    // deltaR jet-lep
    sel_deltaR_jetlep.reset(new JetDeltaRSelection(ctx, 0.4));
    hist_sel_deltaR_jetlep.reset(new AndHists(ctx, "deltaR_jetlep"));
    // btag
    sel_1btag.reset(new NJetSelection(1, 1, id_btag));
    hist_sel_1btag.reset(new AndHists(ctx, "1btag"));
    veto_btag.reset(new NJetSelection(0, 0, id_btag));
    hist_sel_0btag.reset(new AndHists(ctx, "0btag"));
    // deltaR btag/leading jet
    sel_deltaR_leadingbjetlep.reset(new LeadingJetDeltaRSelection(ctx, deltaR_blep_min, ctx.get_handle<vector<Jet> >("btag_loose")));
    hist_sel_1btagdr.reset(new AndHists(ctx, "1btagdr"));
    sel_deltaR_leadingjetlep.reset(new LeadingJetDeltaRSelection(ctx, deltaR_blep_min));
    hist_sel_0btagdr.reset(new AndHists(ctx, "0btagdr"));
    // toptag
    sel_1toptag.reset(new NTopJetSelection(1, 1, id_toptag));
    hist_sel_1btag1toptag.reset(new AndHists(ctx, "1btag1toptag"));
    hist_sel_1btag1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "1btag1toptag_reco", "1TopTagReconstruction", "Chi2"));
    hist_sel_0btag1toptag.reset(new AndHists(ctx, "0btag1toptag"));
    hist_sel_0btag1toptag->add_hist(new BstarToTWHypothesisHists(ctx, "0btag1toptag_reco", "1TopTagReconstruction", "Chi2"));
    veto_toptag.reset(new NTopJetSelection(0, 0, id_toptag));
    hist_sel_1btag0toptag.reset(new AndHists(ctx, "1btag0toptag"));
    hist_sel_1btag0toptag->add_hist(new BstarToTWHypothesisHists(ctx, "1btag0toptag_reco", "0TopTagReconstruction", "Chi2"));
    hist_sel_0btag0toptag.reset(new AndHists(ctx, "0btag0toptag"));
    hist_sel_0btag0toptag->add_hist(new BstarToTWHypothesisHists(ctx, "0btag0toptag_reco", "0TopTagReconstruction", "Chi2"));
    // chi2
    sel_1toprecomass.reset(new RecoMassSelection(ctx, 500., "1TopTagReconstruction"));
    sel_1toptag20chi2.reset(new Chi2Selection(ctx, "1TopTagReconstruction", chi2_max, "Chi2"));
    hist_sel_1btag1toptag20chi2.reset(new AndHists(ctx, "1btag1toptag20chi2"));
    hist_sel_1btag1toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "1btag1toptag20chi2_reco", "1TopTagReconstruction", "Chi2"));
    hist_sel_0btag1toptag20chi2.reset(new AndHists(ctx, "0btag1toptag20chi2"));
    hist_sel_0btag1toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "0btag1toptag20chi2_reco", "1TopTagReconstruction", "Chi2"));
    sel_0toprecomass.reset(new RecoMassSelection(ctx, 500., "0TopTagReconstruction"));
    sel_0toptag20chi2.reset(new Chi2Selection(ctx, "0TopTagReconstruction", chi2_max, "Chi2"));
    hist_sel_1btag0toptag20chi2.reset(new AndHists(ctx, "1btag0toptag20chi2"));
    hist_sel_1btag0toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "1btag0toptag20chi2_reco", "0TopTagReconstruction", "Chi2"));
    hist_sel_0btag0toptag20chi2.reset(new AndHists(ctx, "0btag0toptag20chi2"));
    hist_sel_0btag0toptag20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "0btag0toptag20chi2_reco", "0TopTagReconstruction", "Chi2"));
    // subjet btag (only in signal region)
    sel_subjet_btag.reset(new TopJetDeltaRSelection(ctx, -1, ctx.get_handle<vector<TopJet> >("toptag")));
    hist_sel_subjet_btag.reset(new AndHists(ctx, "1btag1toptag20chi2subjet_btag"));
    hist_sel_subjet_btag->add_hist(new BstarToTWHypothesisHists(ctx, "1btag1toptag20chi2subjet_btag_reco", "1TopTagReconstruction", "Chi2"));

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
    // TO DO!

    // hist_BTagMCEfficiency->fill(event);
    hist_presel->fill(event);

    // deltaPhi met-lep
    if(!sel_deltaPhi_metlep->passes(event)) return false;
    hist_sel_deltaPhi_metlep->fill(event);
    // deltaR jet-lep
    if(!sel_deltaR_jetlep->passes(event)) return false;
    hist_sel_deltaR_jetlep->fill(event);

    // - 0 btag regions - //
    if (veto_btag->passes(event))
      {
	hist_sel_0btag->fill(event);
	// - deltaR jet, lepton - //
	if (!sel_deltaR_leadingjetlep->passes(event)) return false;
	hist_sel_0btagdr->fill(event);
	
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
	    // event.set(h_tot_weight, event.weight);
	    // if (do_for_cr)
	    //   output_module->process(event);
	    return do_for_cr;
	  }
	return false;
      }

    else if (sel_1btag->passes(event)) 
      {
	hist_sel_1btag->fill(event);
	if (!sel_deltaR_leadingbjetlep->passes(event)) return false;
	hist_sel_1btagdr->fill(event);
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
		disc_match->process(event);
		hist_bstar_matchreco->fill(event);
	      }
	    // -- Chi2 Cut -- //
	    if (!sel_1toptag20chi2->passes(event)) return false;
	    hist_sel_1btag1toptag20chi2->fill(event);
	    if (!sel_subjet_btag->passes(event)) return false;
	    hist_sel_subjet_btag->fill(event);
	    // event.set(h_tot_weight, event.weight);
	    // if (!do_for_cr)
	    //   output_module->process(event);
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
    return false;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWAnalysisModule)

}
