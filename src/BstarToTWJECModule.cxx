#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/HOTVRHists.h"
#include "UHH2/BstarToTW/include/HOTVRJetCorrectionHists.h"
#include "UHH2/BstarToTW/include/HOTVRJetCorrector.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWJECModule: public AnalysisModule {
  public: 

    explicit BstarToTWJECModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:

    std::unique_ptr<AnalysisModule> jec_subj_mc, jec_subj_BCD, jec_subj_EFearly, jec_subj_FlateG, jec_subj_H;
    std::unique_ptr<AnalysisModule> cl_topjet;
    std::unique_ptr<Selection> sel_ntop;
    std::unique_ptr<Hists> hist_topcleaner, hist_ntop;
    std::unique_ptr<Hists> hist_jec_corr,hist_pileup;

    bool is_mc;

    const int runnr_BCD = 276811;
    const int runnr_EFearly = 278802;
    const int runnr_FlateG = 280385;

  };

  BstarToTWJECModule::BstarToTWJECModule(Context & ctx) {

    is_mc = ctx.get("dataset_type") == "MC";

    double top_pt_min = 200.0;
    double top_eta_max = 2.5;

    TopJetId id_topjet =  PtEtaCut(top_pt_min, top_eta_max);

    if(is_mc)
      {
    	jec_subj_mc.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_L23_AK4PFchs_MC));
      }
    else
      { 
    	jec_subj_BCD.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_BCD_L23_AK4PFchs_DATA));
    	jec_subj_EFearly.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_EF_L23_AK4PFchs_DATA));
    	jec_subj_FlateG.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_G_L23_AK4PFchs_DATA));
    	jec_subj_H.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_23Sep2016_V4_H_L23_AK4PFchs_DATA));

      }
    hist_jec_corr.reset(new HOTVRJetCorrectionHists(ctx, "JEC_corrections"));

    cl_topjet.reset(new TopJetCleaner(ctx, id_topjet));
    sel_ntop.reset(new NTopJetSelection(1, -1));
    hist_topcleaner.reset(new AndHists(ctx, "Topjet_Cleaning"));
    hist_ntop.reset(new AndHists(ctx, "1TopJetCut"));

    hist_pileup.reset(new HOTVRPileUpHists(ctx, "HOTVR_PileUp"));
  }

  bool BstarToTWJECModule::process(Event & event) {
    // TopJEC
    hist_pileup->fill(event);
    if (is_mc) jec_subj_mc->process(event);
    else
      {
	if(event.run <= runnr_BCD)         jec_subj_BCD->process(event);
	else if(event.run < runnr_EFearly) jec_subj_EFearly->process(event); //< is correct, not <= 
	else if(event.run <= runnr_FlateG) jec_subj_FlateG->process(event);
	else if(event.run > runnr_FlateG)  jec_subj_H->process(event);
      }
    if (is_mc)
      {
	hist_jec_corr->fill(event);
      }
    cl_topjet->process(event);
    hist_topcleaner->fill(event);

    // TopJet Selection
    if(!sel_ntop->passes(event)) return false;
    hist_ntop->fill(event);

    // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWJECModule)

}
