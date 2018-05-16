#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/common/include/CommonModules.h"
#include <UHH2/common/include/JetCorrections.h>
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/TTbarReconstruction.h"

#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/BstarToTW/include/BstarToTWReconstruction.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisHists.h"
#include "UHH2/BstarToTW/include/GenCleaningModules.h"

#include "UHH2/HOTVR/include/HOTVRHists.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"
#include "UHH2/HOTVR/include/HOTVRJetCorrector.h"

#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/TopTagPerformanceHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class TopTagPerformanceModule: public AnalysisModule {
  public:

    explicit TopTagPerformanceModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:

    std::unique_ptr<AnalysisModule> printer;

    Event::Handle<vector<TopJet>> h_ak8jets; // handle for ak8_SoftDrop collection
    std::unique_ptr<CommonModules> common;
    std::unique_ptr<AnalysisModule> genjet_cleaner, hotvr_jetlep_cleaner;
    std::unique_ptr<AnalysisModule> jec_subj_mc, jec_subj_BCD, jec_subj_EFearly, jec_subj_FlateG, jec_subj_H;
    std::unique_ptr<AnalysisModule> jec_topj_mc, jec_topj_BCD, jec_topj_EFearly, jec_topj_FlateG, jec_topj_H;

    std::unique_ptr<AnalysisModule> cl_topjet1, cl_topjet2;
    std::unique_ptr<Selection> sel_ntop1, sel_ntop2;
    std::unique_ptr<AnalysisModule> cl_muon;
    std::unique_ptr<AnalysisModule> BstarToTWgenprod;
    std::unique_ptr<AnalysisModule> primary_lep, w_reco, w_disc;
    std::unique_ptr<Hists> hist_w_reco;
    std::unique_ptr<Selection> sel_toptag1, sel_toptag2;
    std::unique_ptr<AndHists> hist_toptag1, hist_toptag2;


    std::unique_ptr<Selection> trig_IsoMu24, trig_IsoTkMu24;
    std::unique_ptr<Selection> sel_ngenjet,sel_lumi, sel_nmuo, sel_met, sel_1top;

    std::unique_ptr<Hists> hist_pileup;
    bool is_mc, is_qcd;

    const int runnr_BCD = 276811;
    const int runnr_EFearly = 278802;
    const int runnr_FlateG = 280385;

    // --- Scans --- //
    std::unique_ptr<Hists> hist_hotvr_pre, hist_cmstt_pre;
    std::vector<std::unique_ptr<Selection>> sel_hotvr, sel_cmstt;
    std::vector<std::unique_ptr<Hists>> hist_hotvr, hist_cmstt;
    int n_points = 100;

  };

  TopTagPerformanceModule::TopTagPerformanceModule(Context & ctx) {

    printer.reset(new GenParticlesPrinter(ctx)); 

    is_mc = ctx.get("dataset_type") == "MC";
    is_qcd = (ctx.get("dataset_version").find("QCD") == 0 || ctx.get("dataset_version").find("WJets") == 0);

    std::string ak8jets_name = "slimmedJetsAK8_SoftDrop";
    h_ak8jets = ctx.get_handle<vector<TopJet>>(ak8jets_name);

    double deltaPhi_min = M_PI/2;  // minimum delta phi between muon and top
    double top_pt_min = 200.0;
    double top_eta_max = 2.5;
    TopJetId id_topjet =  PtEtaCut(top_pt_min, top_eta_max);
    TopJetId id_topjet2 =  PtEtaCut(400., top_eta_max);
    TopJetId id_wreco  = AndId<TopJet>(id_topjet, DeltaPhiCut(ctx, deltaPhi_min));
    TopJetId id_hotvr_base = HOTVRTopTag(0.8, 140., 220., 50);
    TopJetId id_cmstt_base =  Type2TopTag(105., 220., Type2TopTag::MassType::groomed);
    TopJetId id_toptag1 = AndId<TopJet>(id_hotvr_base, Tau32Groomed(0.56));
    TopJetId id_toptag2 = AndId<TopJet>(id_cmstt_base, Tau32(0.57));

    double met_min     = 50.;
    double lep_eta_max = 2.4;
    double muo_pt_min  = 50.; 
    double muo_iso_max = 0.15;
    double ele_pt_min  = 30.0;
    double jet_pt_min  = 30.0;
    double jet_eta_max = 2.4;
    // MuonId id_muo_loose = AndId<Muon>(MuonIDLoose(), PtEtaCut(muo_pt_min, lep_eta_max), MuonIso(muo_iso_max));
    MuonId id_muo_loose = AndId<Muon>(MuonID(Muon::Selector::CutBasedIdLoose), MuonID(Muon::Selector::PFIsoTight));
    // MuonId id_muo_tight = AndId<Muon>(MuonIDTight(), PtEtaCut(muo_pt_min, lep_eta_max), MuonIso(muo_iso_max));
    MuonId id_muo_tight = AndId<Muon>(MuonID(Muon::Selector::CutBasedIdTight), MuonID(Muon::Selector::PFIsoTight));
    ElectronId id_ele = AndId<Electron>(ElectronID_Spring16_veto_noIso, PtEtaCut(ele_pt_min, lep_eta_max));
    JetId id_jet = AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_LEPVETO), PtEtaCut(jet_pt_min, jet_eta_max));

    common.reset(new CommonModules());
    common->switch_jetlepcleaner(true);
    common->set_muon_id(id_muo_tight);
    common->set_electron_id(id_ele);
    common->set_jet_id(id_jet);
    common->init(ctx);

    genjet_cleaner.reset(new GenJetCleaner(ctx, 100., 2.5));
    sel_ngenjet.reset(new NGenJetSelection(1, -1));

    cl_muon.reset(new MuonCleaner(id_muo_tight));

    trig_IsoMu24.reset(new TriggerSelection("HLT_IsoMu24_v*"));
    trig_IsoTkMu24.reset(new TriggerSelection("HLT_IsoTkMu24_v*"));
    sel_lumi.reset(new LumiSelection(ctx));
    sel_nmuo.reset(new NMuonSelection(1, 1));
    sel_met.reset(new METSelection(met_min));
    sel_1top.reset(new NTopJetSelection(1,1,id_wreco));
    // sel_toptag1.reset(new NTopJetSelection(1,1,id_toptag1));
    // sel_toptag2.reset(new NTopJetSelection(1,1,id_toptag2));
    // hist_toptag1.reset(new AndHists(ctx,,id_toptag1));
    // hist_toptag2.reset(new AndHists(1,1,id_toptag2));
		      

    primary_lep.reset(new PrimaryLepton(ctx));
    w_reco.reset(new BstarToTWReconstruction(ctx, NeutrinoReconstruction, "WReconstruction", id_wreco));
    w_disc.reset(new BstarToTWChi2Discriminator(ctx, "WReconstruction"));
    hist_w_reco.reset(new BstarToTWHypothesisHists(ctx, "WReco", "WReconstruction", "Chi2"));

    cl_topjet1.reset(new TopJetCleaner(ctx, id_topjet));
    cl_topjet2.reset(new TopJetCleaner(ctx, id_topjet, ak8jets_name));
    sel_ntop1.reset(new NTopJetSelection(1, -1));
    sel_ntop2.reset(new NTopJetSelection(1, -1, boost::none, h_ak8jets));

    if(is_mc)
      {
    	jec_subj_mc.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L23_AK4PFchs_MC));

    	jec_topj_mc.reset(new GenericTopJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L123_AK8PFchs_MC, ak8jets_name));
      }
    else
      { 
    	jec_subj_BCD.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_BCD_L23_AK4PFchs_DATA));
    	jec_subj_EFearly.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_EF_L23_AK4PFchs_DATA));
    	jec_subj_FlateG.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_G_L23_AK4PFchs_DATA));
    	jec_subj_H.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_H_L23_AK4PFchs_DATA));

    	jec_topj_BCD.reset(new GenericTopJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK8PFchs_DATA, ak8jets_name));
    	jec_topj_EFearly.reset(new GenericTopJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_EF_L123_AK8PFchs_DATA, ak8jets_name));
    	jec_topj_FlateG.reset(new GenericTopJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_G_L123_AK8PFchs_DATA, ak8jets_name));
    	jec_topj_H.reset(new GenericTopJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_H_L123_AK8PFchs_DATA, ak8jets_name));

      }

    hist_pileup.reset(new HOTVRPileUpHists(ctx, "HOTVR_PileUp"));
    // --- Scans --- //
    // TopJetId id_hotvr_base = AndId<TopJet>(HOTVRTopTag(0.8, 140., 220., 50), DeltaPhiCut(deltaPhi_min));
    // TopJetId id_cmstt_base =  AndId<TopJet>(Type2TopTag(105., 220., Type2TopTag::MassType::groomed), DeltaPhiCut(deltaPhi_min));   
    // TopJetId id_cmstt_base =  AndId<TopJet>(id_topjet2, Type2TopTag(105., 220., Type2TopTag::MassType::groomed), DeltaPhiCut(deltaPhi_min));


    for (int i = 0; i < n_points; ++i)
      {
    	std::unique_ptr<Selection> sel;
    	std::unique_ptr<Hists> hist;
	TopJetId id;
	double tau32 = 1 - (i * 1./n_points);
	// - HOTVR - //
    	id = AndId<TopJet>(id_hotvr_base, Tau32Groomed(tau32));
    	sel.reset(new NTopJetSelection(1, -1, id));
    	sel_hotvr.push_back(std::move(sel));

    	hist.reset(new TopTagPerformanceHists(ctx, "HOTVR_Performance_" + to_string(i), is_qcd,  id));
    	hist_hotvr.push_back(std::move(hist));

	// - CMSTT - //
    	id = AndId<TopJet>(id_cmstt_base, Tau32(tau32));
    	sel.reset(new NTopJetSelection(1, -1, id, h_ak8jets));
    	sel_cmstt.push_back(std::move(sel));

    	hist.reset(new TopTagPerformanceHists(ctx, "CMSTT_Performance_" + to_string(i), is_qcd, id, h_ak8jets));
    	hist_cmstt.push_back(std::move(hist));
      }


    hist_hotvr_pre.reset(new TopTagPerformanceHists(ctx, "HOTVR_Pre", is_qcd, id_topjet));
    hist_cmstt_pre.reset(new TopTagPerformanceHists(ctx, "CMSTT_Pre", is_qcd, id_topjet, h_ak8jets));

  }

  bool TopTagPerformanceModule::process(Event & event) {

    // printer->process(event);
   if(!is_mc)
      {
	if(!sel_lumi->passes(event)) return false;
      }
    if(!common->process(event)) return false;
    if(is_qcd) genjet_cleaner->process(event);
    // genjet_cleaner->process(event);
    // if(is_qcd && !sel_ngenjet->passes(event)) return false;

    // cl_muon->process(event);
    // if(!sel_nmuo->passes(event)) return false;

    if (is_mc)
      {
	jec_subj_mc->process(event);
	jec_topj_mc->process(event);
      }
    else
      {
	if(event.run <= runnr_BCD)         
	  {
	    jec_subj_BCD->process(event);
	    jec_topj_BCD->process(event);
	  }
	else if(event.run < runnr_EFearly) //< is correct, not <= 
	  {
	    jec_subj_EFearly->process(event);
	    jec_topj_EFearly->process(event);
	  }
	else if(event.run <= runnr_FlateG) 
	  {
	    jec_subj_FlateG->process(event);
	    jec_topj_FlateG->process(event);
	  }
	else if(event.run > runnr_FlateG)  
	  {
	    jec_subj_H->process(event);
	    jec_topj_H->process(event);
	  }
      }


    cl_topjet1->process(event);
    cl_topjet2->process(event);
    if(is_mc)
      {
	bool sel_pt_hotvr = sel_ntop1->passes(event);
	// if(sel_pt_hotvr) 
	hist_hotvr_pre->fill(event);
	bool sel_pt_cmstt = sel_ntop2->passes(event);
	// if(sel_pt_cmstt)
	hist_cmstt_pre->fill(event);

	if (!sel_pt_hotvr && !sel_pt_cmstt) return false;
	for (int i = 0; i < n_points; ++i)
	  {
	    if (sel_pt_hotvr)
	      {
		if (sel_hotvr.at(i)->passes(event)) hist_hotvr.at(i)->fill(event);
	      }
	    if (sel_pt_cmstt)
	      {
		if (sel_cmstt.at(i)->passes(event)) hist_cmstt.at(i)->fill(event);
	      }
	  }
      }


    if(!(trig_IsoMu24->passes(event) || trig_IsoTkMu24->passes(event))) return false;
    primary_lep->process(event);
    if(!sel_nmuo->passes(event)) return false;
    if(!sel_met->passes(event)) return false;

    if(sel_1top->passes(event))
      {
	w_reco->process(event);
	w_disc->process(event);
	hist_w_reco->fill(event);
      }

    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(TopTagPerformanceModule)
}
