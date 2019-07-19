#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/JetCorrectionSets.h>
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

    Event::Handle<vector<TopJet>> h_ak8jets; // handle for ak8_SoftDrop collection
    std::unique_ptr<CommonModules> common;
    std::unique_ptr<AnalysisModule> genjet_cleaner, hotvr_jetlep_cleaner;
    std::unique_ptr<AnalysisModule> jec_subj_mc;
    std::unique_ptr<AnalysisModule> jec_topj_mc;

    std::unique_ptr<AnalysisModule> cl_topjet_hotvr, cl_topjet_ak8;
    std::unique_ptr<Selection> sel_ntop_hotvr, sel_ntop_ak8;

    std::unique_ptr<Selection> sel_toptag_hotvr, sel_toptag_softdrop;
    std::unique_ptr<AndHists> hist_toptag_hotvr, hist_toptag_softdrop;

    std::unique_ptr<Selection> sel_trigger, sel_muo, sel_ngenjet;

    std::unique_ptr<Hists> hist_pileup;
    bool is_mc, is_qcd;

    // --- Scans --- //
    std::unique_ptr<Hists> hist_hotvr_pre, hist_softdrop_pre;
    std::vector<std::unique_ptr<Selection>> sel_hotvr, sel_softdrop;
    std::vector<std::unique_ptr<Hists>> hist_hotvr, hist_softdrop;
    int n_points = 100;

  };

  TopTagPerformanceModule::TopTagPerformanceModule(Context & ctx) {

    is_mc = ctx.get("dataset_type") == "MC";
    is_qcd = (ctx.get("dataset_version").find("QCD") == 0 || ctx.get("dataset_version").find("WJets") == 0);

    std::string ak8jets_name = "jetsAk8PuppiSubstructure_SoftDropPuppi";
    h_ak8jets = ctx.get_handle<vector<TopJet>>(ak8jets_name);

    double top_pt_min = 200.0;
    double top_eta_max = 2.5;
    TopJetId id_topjet =  PtEtaCut(top_pt_min, top_eta_max);
    TopJetId id_hotvr_base = HOTVRTopTag(0.8, 140., 220., 50);
    TopJetId id_softdrop_base =  Type2TopTag(105., 220., Type2TopTag::MassType::groomed);
    TopJetId id_hotvr = AndId<TopJet>(id_hotvr_base, Tau32Groomed(0.56));
    TopJetId id_softdrop = AndId<TopJet>(id_softdrop_base, Tau32(0.57));

    double jet_pt_min  = 30.0;
    double jet_eta_max = 2.4;
    double lep_pt_min  = 50.0;
    double lep_eta_max = 2.4;
    JetId id_jet = AndId<Jet>(JetPFID(JetPFID::WP_TIGHT_CHS), PtEtaCut(jet_pt_min, jet_eta_max));
    MuonId id_muo = AndId<Muon>(MuonID(Muon::Selector::CutBasedIdTight), PtEtaCut(lep_pt_min, lep_eta_max));

    common.reset(new CommonModules());
    common->set_jet_id(id_jet);
    common->set_muon_id(id_muo);
    common->disable_jec();
    common->disable_jersmear();
    common->init(ctx);

    genjet_cleaner.reset(new GenJetCleaner(ctx, jet_pt_min, jet_eta_max));
    sel_ngenjet.reset(new NGenJetSelection(1, -1));

    sel_muo.reset(new NMuonSelection(1,-1));

    cl_topjet_hotvr.reset(new TopJetCleaner(ctx, id_topjet));
    cl_topjet_ak8.reset(new TopJetCleaner(ctx, id_topjet, ak8jets_name));
    sel_ntop_hotvr.reset(new NTopJetSelection(1, -1));
    sel_ntop_ak8.reset(new NTopJetSelection(1, -1, boost::none, h_ak8jets));

    jec_subj_mc.reset(new HOTVRJetCorrector(ctx, JERFiles::Autumn18_V8_L123_AK4PFPuppi_MC));
    jec_topj_mc.reset(new GenericTopJetCorrector(ctx, JERFiles::Autumn18_V8_L123_AK8PFPuppi_MC, ak8jets_name));

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

	// - Softdrop - //
    	id = AndId<TopJet>(id_softdrop_base, Tau32(tau32));
    	sel.reset(new NTopJetSelection(1, -1, id, h_ak8jets));
    	sel_softdrop.push_back(std::move(sel));

    	hist.reset(new TopTagPerformanceHists(ctx, "Softdrop_Performance_" + to_string(i), is_qcd, id, h_ak8jets));
    	hist_softdrop.push_back(std::move(hist));
      }


    hist_hotvr_pre.reset(new TopTagPerformanceHists(ctx, "HOTVR_Pre", is_qcd, id_topjet));
    hist_softdrop_pre.reset(new TopTagPerformanceHists(ctx, "Softdrop_Pre", is_qcd, id_topjet, h_ak8jets));

  }

  bool TopTagPerformanceModule::process(Event & event) {

    if(!common->process(event)) return false;
    if(is_qcd) genjet_cleaner->process(event);

    // if (!sel_trigger->passes(event) && !sel_muo->passes(event)) return false;
    if (!sel_muo->passes(event)) return false;
    jec_subj_mc->process(event);
    jec_topj_mc->process(event);
    
    cl_topjet_hotvr->process(event);
    cl_topjet_ak8->process(event);

    bool sel_pt_hotvr = sel_ntop_hotvr->passes(event);
    hist_hotvr_pre->fill(event);
    bool sel_pt_softdrop = sel_ntop_ak8->passes(event);
    hist_softdrop_pre->fill(event);

    if (!sel_pt_hotvr && !sel_pt_softdrop) return false;
    for (int i = 0; i < n_points; ++i)
      {
	if (sel_pt_hotvr)
	  {
	    if (sel_hotvr.at(i)->passes(event)) hist_hotvr.at(i)->fill(event);
	  }
	if (sel_pt_softdrop)
	  {
	    if (sel_softdrop.at(i)->passes(event)) hist_softdrop.at(i)->fill(event);
	  }
      }

    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(TopTagPerformanceModule)
}
