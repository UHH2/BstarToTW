#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"

#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/JetHists.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWGenSelections.h"
#include "UHH2/BstarToTW/include/HOTVRIds.h"
#include "UHH2/BstarToTW/include/HOTVRHists.h"
#include "UHH2/BstarToTW/include/HOTVRPerformanceHists.h"
#include "UHH2/BstarToTW/include/EfficiencyHists.h"

#include <vector>

using namespace std;
using namespace uhh2;

namespace uhh2 {

  /** \brief Module for comparing HOTVR with other top taggers.
   *
   * The efficiency of different top taggers are tested by applying
   * the corresponding TopJetIds and requiring ==1 tagged topjet. The
   * number of events passing this selection will then be filled into
   * a histogram, to calculate efficiencies.
   *
   */
  class HOTVRPerformanceModule: public AnalysisModule {
  public:
    
    explicit HOTVRPerformanceModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    std::unique_ptr<AnalysisModule> BstarToTWgenprod; // gen particle interpreter
    std::unique_ptr<Selection> sel_muonchannel; // select only semi-leptonic events on generator level.
    Event::Handle<vector<TopJet>> h_ak8jets; // handle for ak8_SoftDrop collection

    std::unique_ptr<Hists> all;

    // Cleaner
    TopJetId id_pt_eta;
    TopJetId id_jetlep;
    MuonId  id_muon;
    
    std::unique_ptr<AnalysisModule> cl_hotvr_pt_eta, cl_ak8_pt_eta;
    std::unique_ptr<AnalysisModule> cl_hotvr_jetlep, cl_ak8_jetlep;
    std::unique_ptr<AnalysisModule> cl_muo;

    // HOTVR Top Tag
    TopJetId id_hotvrtoptag;
    std::unique_ptr<Selection> sel_hotvrtoptag;
    std::unique_ptr<AnalysisModule> cl_hotvrtoptag;
    std::unique_ptr<Hists> hist_hotvrtoptag;
    std::unique_ptr<Hists> hist_hotvrtoptag_performance;
    std::unique_ptr<Hists> hist_hotvr;


    //CMS Top Tagger
    std::unique_ptr<Selection> sel_toptag;
    // // Top Tag 3%
    TopJetId id_toptag3;
    std::unique_ptr<AnalysisModule> cl_toptag3;
    std::unique_ptr<Hists> hist_toptag3;
    std::unique_ptr<Hists> hist_toptag3_performance;

    // // Top Tag 1%
    TopJetId id_toptag1;
    std::unique_ptr<AnalysisModule> cl_toptag1;
    std::unique_ptr<Hists> hist_toptag1;
    std::unique_ptr<Hists> hist_toptag1_performance;
 
    // // Top Tag 0.3%
    TopJetId id_toptag03;
    std::unique_ptr<AnalysisModule> cl_toptag03;
    std::unique_ptr<Hists> hist_toptag03;
    std::unique_ptr<Hists> hist_toptag03_performance;

    // // Top Tag 0.1%
    TopJetId id_toptag01;
    std::unique_ptr<AnalysisModule> cl_toptag01;
    std::unique_ptr<Hists> hist_toptag01;
    std::unique_ptr<Hists> hist_toptag01_performance;
 };

  HOTVRPerformanceModule::HOTVRPerformanceModule(Context & ctx) {
 
    BstarToTWgenprod.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    sel_muonchannel.reset(new MuonChannelSelection(ctx));
    
    std::string ak8jets_name = "slimmedJetsAK8_SoftDrop";
    h_ak8jets = ctx.get_handle<vector<TopJet>>(ak8jets_name);

    all.reset(new EfficiencyHists(ctx, "all"));

    double top_pt_min = 200.0;
    double top_eta_max = 2.5;

    id_pt_eta = PtEtaCut(top_pt_min, top_eta_max);
    // id_jetlep = GenDeltaPhiCut(ctx, M_PI/2);

    cl_hotvr_pt_eta.reset(new TopJetCleaner(ctx, id_pt_eta));
    cl_hotvr_jetlep.reset(new TopJetCleaner(ctx, GenDeltaPhiCut(ctx, M_PI/2)));
    cl_ak8_pt_eta.reset(new TopJetCleaner(ctx, id_pt_eta, ak8jets_name));
    cl_ak8_jetlep.reset(new TopJetCleaner(ctx, GenDeltaPhiCut(ctx, M_PI/2), ak8jets_name));


  // --- HOTVR
    id_hotvrtoptag = HOTVRTopTag(); // default settings for HOTVRTopTag, see HOTVRId.h

    sel_hotvrtoptag.reset(new NTopJetSelection(1, 1));
    cl_hotvrtoptag.reset(new TopJetCleaner(ctx, id_hotvrtoptag));
    hist_hotvrtoptag.reset(new EfficiencyHists(ctx, "HOTVRTopTag_Efficiencies"));
    hist_hotvrtoptag_performance.reset(new HOTVRPerformanceHists(ctx, "HOTVRTopTag_Performance"));
    hist_hotvr.reset(new HOTVRHists(ctx, "HOTVR"));

    // --- CMS Top Tagger
    sel_toptag.reset(new NTopJetSelection(1,1, boost::none, h_ak8jets));
   // TopTag 3%
    id_toptag3 = AndId<TopJet>(Type2TopTag(105,220, Type2TopTag::MassType::groomed), Tau32(0.81));

    cl_toptag3.reset(new TopJetCleaner(ctx, id_toptag3, ak8jets_name));
    hist_toptag3.reset(new EfficiencyHists(ctx, "TopTag3_Efficiencies", h_ak8jets));
    hist_toptag3_performance.reset(new HOTVRPerformanceHists(ctx, "TopTag3_Performance", h_ak8jets));

  //   // TopTag 1%
    id_toptag1 = AndId<TopJet>(Type2TopTag(105,220, Type2TopTag::MassType::groomed), Tau32(0.67));

    cl_toptag1.reset(new TopJetCleaner(ctx, id_toptag1, ak8jets_name));
    hist_toptag1.reset(new EfficiencyHists(ctx, "TopTag1_Efficiencies", h_ak8jets));
    hist_toptag1_performance.reset(new HOTVRPerformanceHists(ctx, "TopTag1_Performance", h_ak8jets));

  //   // TopTag 0.3%
    id_toptag03 = AndId<TopJet>(Type2TopTag(105,220, Type2TopTag::MassType::groomed), Tau32(0.54));

    cl_toptag03.reset(new TopJetCleaner(ctx, id_toptag03, ak8jets_name));
    hist_toptag03.reset(new EfficiencyHists(ctx, "TopTag0.3_Efficiencies", h_ak8jets));
    hist_toptag03_performance.reset(new HOTVRPerformanceHists(ctx, "TopTag0.3_Performance", h_ak8jets));
  //   // TopTag 0.1%
    id_toptag01 = AndId<TopJet>(Type2TopTag(105,220, Type2TopTag::MassType::groomed), Tau32(0.46));

    cl_toptag01.reset(new TopJetCleaner(ctx, id_toptag01, ak8jets_name));
    hist_toptag01.reset(new EfficiencyHists(ctx, "TopTag0.1_Efficiencies", h_ak8jets));
    hist_toptag01_performance.reset(new HOTVRPerformanceHists(ctx, "TopTag0.1_Performance", h_ak8jets));
  }

  bool HOTVRPerformanceModule::process(Event & event) {
    BstarToTWgenprod->process(event);

    // // Only semi-leptonic channel
    if(!sel_muonchannel->passes(event)) return false;
    all->fill(event);

    cl_hotvr_pt_eta->process(event);
    cl_hotvr_jetlep->process(event);
    cl_ak8_pt_eta->process(event);
    cl_ak8_jetlep->process(event);
    
    cl_hotvrtoptag->process(event);
    if(sel_hotvrtoptag->passes(event))
      {
	hist_hotvrtoptag->fill(event);
	hist_hotvrtoptag_performance->fill(event);
	hist_hotvr->fill(event);
      }
      
    cl_toptag3->process(event);
    if(sel_toptag->passes(event))
      {
	hist_toptag3->fill(event);
	hist_toptag3_performance->fill(event);
      }
    cl_toptag1->process(event);
    if(sel_toptag->passes(event)) 
      {
	hist_toptag1->fill(event);
	hist_toptag1_performance->fill(event);
      }
    cl_toptag03->process(event);
    if(sel_toptag->passes(event))
      {
	hist_toptag03->fill(event);
	hist_toptag03_performance->fill(event);
      }
    cl_toptag01->process(event);
    if(sel_toptag->passes(event))
      {
	hist_toptag01->fill(event);
	hist_toptag01_performance->fill(event);
      }

    // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(HOTVRPerformanceModule)

}
