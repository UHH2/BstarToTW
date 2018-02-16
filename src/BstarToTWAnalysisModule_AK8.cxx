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

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWModules.h"
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

    // --- Additional Modules --- //
    // Common Modules
    std::unique_ptr<CommonModules> common;
    // Bstar Reconsturction
    std::unique_ptr<AnalysisModule> BstarToTWgenprod;
    std::unique_ptr<AnalysisModule> primary_lep, bstar_reco, bstar_disc, bstar_matchdisc;

    // Scale factors & Uncertainties
    std::unique_ptr<AnalysisModule> sf_muon_id, sf_muon_iso, sf_muon_trigger, sf_muon_trk;
    std::string sys_muon_id, sys_muon_iso, sys_muon_trigger, sys_muon_trk;
    std::unique_ptr<AnalysisModule> sf_btag_tight;
    std::string sys_btag_tight;
    std::unique_ptr<AnalysisModule> sf_top_pt_reweight;
    bool do_top_pt_reweight;
    std::unique_ptr<AnalysisModule> sf_toptag;
    std::string sys_toptag;
    std::unique_ptr<AnalysisModule> scale_variation;
    bool do_scale_variation, do_pdf_variations;
    string sys_pu;



    std::unique_ptr<Selection> sel_st;

    // --- Selections and Histogramms --- //    
    std::unique_ptr<AnalysisModule> cl_topjet;
    std::unique_ptr<Selection> sel_ntop;
    std::unique_ptr<AndHists> hist_presel; 

    // - toptag - //
    std::unique_ptr<Selection> sel_1toptag;
    std::unique_ptr<AndHists> hist_sel_1toptag;
    // ttbar CR //
    std::unique_ptr<Selection> sel_2btag_tight;
    std::unique_ptr<AndHists> hist_sel_2btag, hist_sel_not2btag, hist_toptag_efficiency, hist_toptag_sf_closure;

    // - chi2 - //
    std::unique_ptr<Selection> sel_20chi2;
    std::unique_ptr<AndHists> hist_sel_20chi2;

    // - additional hists - //
    std::unique_ptr<Hists> hist_bstar_reco, hist_bstar_matchreco;
    std::unique_ptr<Hists> hist_BTagMCEfficiency;

    // --- Scans --- //
    std::vector<std::unique_ptr<AnalysisModule>> rec_tau32, disc_tau32;
    std::vector<std::unique_ptr<Selection>> sel_tau32, sel_chi2;
    std::vector<std::unique_ptr<Hists>> hist_tau32, hist_chi2;
    int n_tau32, n_chi2;
    

  };

  BstarToTWAnalysisModule_AK8::BstarToTWAnalysisModule_AK8(Context & ctx) {


    dataset_name = ctx.get("dataset_version");

    
    do_scale_variation = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down") || (ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");
    do_pdf_variations  = ctx.get("b_PDFUncertainties") == "true";
    do_top_pt_reweight = ctx.get("b_TopPtReweight") == "true";

    sys_muon_id      = ctx.get("Systematic_MuonID");
    sys_muon_trigger = ctx.get("Systematic_MuonTrigger");
    sys_muon_iso     = ctx.get("Systematic_MuonIso");
    sys_muon_trk     = ctx.get("Systematic_MuonTrk");
    sys_btag_tight   = ctx.get("Systematic_BTag");
    sys_toptag       = ctx.get("Systematic_TopTag");
    sys_pu           = ctx.get("Systematic_PU");

    // --- Selection Variables --- //
    double deltaPhi_min = M_PI/2;  // minimum delta phi between muon and top

    double top_m_min     = 105.;   // minimum topjet mass
    double top_m_max     = 220.;   // maximum topjet mass
    double top_pt_min    = 400.0;
    double top_eta_max   = 2.5;

    double chi2_max      = 20;     // maximum chi2 of the reconstructed bstar hypothesis


    // - Object IDs - //
    TopJetId id_topjet =  PtEtaCut(top_pt_min, top_eta_max);
    cl_topjet.reset(new TopJetCleaner(ctx, id_topjet));
    sel_ntop.reset(new NTopJetSelection(1, -1));

    double tau32_wp = .57;
    TopJetId id_toptag =  AndId<TopJet>(Type2TopTag(top_m_min, top_m_max, Type2TopTag::MassType::groomed), Tau32(tau32_wp), DeltaPhiCut(ctx, deltaPhi_min));  

    CSVBTag::wp btag_wp_tight = CSVBTag::WP_TIGHT;
    JetId id_btag_tight  = CSVBTag(btag_wp_tight);


    // --- Setup Additional Modules --- //
    // - Common Modules - //
    common.reset(new CommonModules());
    common->disable_jec();
    common->disable_jersmear();
    common->init(ctx, sys_pu);

    sel_st.reset(new STSelection(400));

    // - Bstar reconstruction - //
    primary_lep.reset(new PrimaryLepton(ctx));
    bstar_reco.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "BstarToTWReconstruction", id_toptag));
    bstar_disc.reset(new BstarToTWChi2Discriminator(ctx, "BstarToTWReconstruction"));
    if (dataset_name.find("BstarToTW") == 0)
      {
    	// BstarToTWgenprod.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    	// bstar_matchdisc.reset(new BstarToTWMatchDiscriminator(ctx, "BstarToTWReconstruction"));
      }
    hist_bstar_reco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWReco", "BstarToTWReconstruction", "Chi2"));
    // hist_bstar_matchreco.reset(new BstarToTWHypothesisHists(ctx, "BstarToTWMatchedReco", "BstarToTWReconstruction", "Match"));

    // - Scale Factors & Uncertainties - //
    sf_muon_id.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1, "tightID", true, sys_muon_id));
    sf_muon_iso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1, "iso", true, sys_muon_iso));
    sf_muon_trigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger", true, sys_muon_trigger));
    sf_muon_trk.reset(new MCMuonTrkScaleFactor(ctx, "/nfs/dust/cms/user/froehlia/CMSSW_8_0_24_patch1/src/UHH2/common/data/Tracking_EfficienciesAndSF_BCDEFGH.root", 1, "track", sys_muon_trk, "muons"));

    sf_btag_tight.reset(new MCBTagScaleFactor(ctx, btag_wp_tight, "jets", sys_btag_tight, "mujets", "incl", "MCBtagEfficienciesTight"));

    sf_top_pt_reweight.reset(new TopPtReweight(ctx, 0.0615, -0.0005, "", "", do_top_pt_reweight, 1.0)); // https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting

    sf_toptag.reset(new CMSTTScaleFactor(ctx, "BstarToTW", sys_toptag));
   
    scale_variation.reset(new MCScaleVariation(ctx));

    // --- Selections and Histogramms --- //
    hist_presel.reset(new AndHists(ctx, "PreSel"));
    // hist_presel->add_hist(new HOTVRHists(ctx, "PreSel_HOTVRtagged", id_toptag));

    // - TopTag - //
    sel_1toptag.reset(new NTopJetSelection(1, 1, id_toptag));
    hist_sel_1toptag.reset(new AndHists(ctx, "1TopTag", id_toptag));
    // TopTag CR //

    sel_2btag_tight.reset(new NJetSelection(2, -1, id_btag_tight));
    hist_sel_2btag.reset(new AndHists(ctx, "2btag"));
    hist_sel_not2btag.reset(new AndHists(ctx, "less2btag"));
    // - Chi2 - //
    sel_20chi2.reset(new Chi2Selection(ctx, "BstarToTWReconstruction", chi2_max));
    hist_sel_20chi2.reset(new AndHists(ctx, "20Chi2", id_toptag));
    hist_sel_20chi2->add_hist(new BstarToTWHypothesisHists(ctx, "20Chi2_reco", "BstarToTWReconstruction", "Chi2"));

    // - Scale Factor Hists - //
    hist_BTagMCEfficiency.reset(new BTagMCEfficiencyHists(ctx, "BTagMCEfficiency", btag_wp_tight));    
    // hist_toptag_efficiency.reset(new AndHists(ctx, "TopTagEfficiency"));   
    // hist_toptag_efficiency->add_hist(new HOTVRHists(ctx, "TopTagEfficiency_HOTVRtagged", id_toptag));
    // hist_toptag_sf_closure.reset(new AndHists(ctx, "TopTagSFClosure"));
    // hist_toptag_sf_closure->add_hist(new HOTVRHists(ctx, "TopTagSFClosure_HOTVRtagged", id_toptag));


    // --- Scans --- //

    // - tau_3/2 Scan - //
    // n_tau32 = 10;
    // for (int i = 0; i < n_tau32; ++i)
    //   {
    // 	std::unique_ptr<Selection> sel, sel2;
    // 	TopJetId id = AndId<TopJet>(HOTVRTopTag(top_fpt_max, 105 + i*5, 220, top_mpair_min), DeltaPhiCut(ctx, deltaPhi_min), Tau32Groomed(0.56));
    // 	sel.reset(new NTopJetSelection(1, 1, id));
    // 	sel_tau32.push_back(std::move(sel));

    // 	std::string name = "Reconstruction" + to_string(i);

    // 	std::unique_ptr<AnalysisModule> rec;
    // 	rec.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, name, id));
    // 	rec_tau32.push_back(std::move(rec));

    // 	std::unique_ptr<AnalysisModule> disc;
    // 	disc.reset(new BstarToTWChi2Discriminator(ctx, name));
    // 	disc_tau32.push_back(std::move(disc));

    // 	std::unique_ptr<Hists> hist;
    // 	hist.reset(new BstarToTWHypothesisHists(ctx, "Tau32_Reco" + to_string(i), name, "Chi2"));
    // 	hist_tau32.push_back(std::move(hist));
    //   }

    // - Chi2 Scan - //
    // n_chi2 = 10;
    // double chi2[n_chi2] = {1, 5, 10, 15, 20, 30, 50, 75, 100, 200};
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

  bool BstarToTWAnalysisModule_AK8::process(Event & event) {
    // -- after PreSel -- //
    // apply scale factors for muons and trigger
    sf_muon_id->process(event);
    sf_muon_iso->process(event);
    sf_muon_trigger->process(event);
    sf_muon_trk->process(event);
    common->process(event);

    cl_topjet->process(event);
    if(!sel_ntop->passes(event)) return false;

    if(do_scale_variation) scale_variation->process(event);
    if(do_top_pt_reweight) sf_top_pt_reweight->process(event);   

    if(!sel_st->passes(event)) return false;
    hist_presel->fill(event);


    // -- Tau32 Scan -- //
    // for (int i = 0; i < n_tau32; ++i)
    //   {
    // 	if (sel_tau32.at(i)->passes(event))
    // 	  {
    // 	    rec_tau32.at(i)->process(event);
    // 	    disc_tau32.at(i)->process(event);
    // 	    hist_tau32.at(i)->fill(event);
    // 	  }
    //   }


    // -- TopTag Cut -- //    
    // Fill ttbar controll region
    hist_BTagMCEfficiency->fill(event);
    sf_btag_tight->process(event);
    if(sel_2btag_tight->passes(event)) 
      {
	hist_sel_2btag->fill(event);
	// hist_toptag_efficiency->fill(event);
       
	// sf_toptag->process(event);
	// hist_toptag_sf_closure->fill(event);
	return false;
      }

    hist_sel_not2btag->fill(event);
    // Cut on TopTag
    if(!sel_1toptag->passes(event)) return false;
    sf_toptag->process(event);
    hist_sel_1toptag->fill(event);

    // // -- Bstar Reconstructinon -- //	
    primary_lep->process(event);
    bstar_reco->process(event);
    bstar_disc->process(event);
    hist_bstar_reco->fill(event);
    if (dataset_name.find("BstarToTW") == 0)
      {
	// BstarToTWgenprod->process(event);
	// bstar_matchdisc->process(event);
	// hist_bstar_matchreco->fill(event);
      }

    // -- Chi2 Scan -- //
    // for (int i = 0; i < n_chi2; ++i)
    //   {
    // 	if (sel_chi2.at(i)->passes(event))
    // 	  {
    // 	    hist_chi2.at(i)->fill(event);
    // 	  }
    //   }


    // -- Chi2 Cut -- //
    if (!sel_20chi2->passes(event)) return false;
    hist_sel_20chi2->fill(event);
   
    // -- Done -- //
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWAnalysisModule_AK8)

}
