#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"

#include "UHH2/BstarToTW/include/HOTVRIds.h"
#include "UHH2/BstarToTW/include/HOTVRPerformanceHists.h"
#include "UHH2/BstarToTW/include/HOTVRJetCorrector.h"
#include "UHH2/BstarToTW/include/AndHists.h"

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

    Event::Handle<vector<TopJet>> h_ak8jets; // handle for ak8_SoftDrop collection
    std::unique_ptr<AnalysisModule> jec_hotvr;

    // Cleaner
    std::unique_ptr<AnalysisModule> cl_muo, cl_hotvr, cl_ak8;

    // Selections
    std::unique_ptr<Selection> sel_muo;
    std::unique_ptr<Selection> tt_hotvr, tt_cms3, tt_cms1, tt_cms03, tt_cms01;

    // Histograms
    std::unique_ptr<Hists> hist_sel_muo;
    std::unique_ptr<Hists> hist_tt_hotvr, hist_tt_cms3, hist_tt_cms1, hist_tt_cms03, hist_tt_cms01;
    
  };

  HOTVRPerformanceModule::HOTVRPerformanceModule(Context & ctx) {

    std::string ak8jets_name = "slimmedJetsAK8_SoftDrop";
    h_ak8jets = ctx.get_handle<vector<TopJet>>(ak8jets_name);

    jec_hotvr.reset(new HOTVRJetCorrector(ctx));

    // Kinematic variables
    double muo_pt_min = 50.0;
    double muo_eta_max = 2.4;
    double muo_iso_max = 0.15;
   
    double top_pt_min = 200.0;
    double top_eta_max = 2.5;

    // Cleaner
    MuonId id_muo = AndId<Muon>(MuonIDTight(), PtEtaCut(muo_pt_min, muo_eta_max), MuonIso(muo_iso_max));
    cl_muo.reset(new MuonCleaner(id_muo));
    
    
    // TopJetId id_topjet = AndId<TopJet>(PtEtaCut(top_pt_min, top_eta_max), DeltaPhiCut(M_PI/2));
    TopJetId id_topjet = PtEtaCut(top_pt_min, top_eta_max);
    cl_hotvr.reset(new TopJetCleaner(ctx, id_topjet));
    cl_ak8.reset(new TopJetCleaner(ctx, id_topjet, ak8jets_name));

    // Selections
    sel_muo.reset(new NMuonSelection(1, -1));
    hist_sel_muo.reset(new AndHists(ctx, "MuonCut"));

    // HOTVR Top Tagger
    TopJetId id_hotvr = AndId<TopJet>(HOTVRTopTag(), Tau32Groomed(0.81));
    tt_hotvr.reset(new NTopJetSelection(1, -1, id_hotvr));
    hist_tt_hotvr.reset(new AndHists(ctx, "HOTVR"));

    // CMS Top Tagger
    // 3% misstag
    TopJetId id_cms3 = AndId<TopJet>(Type2TopTag(105,220, Type2TopTag::MassType::groomed), Tau32(0.81));
    tt_cms3.reset(new NTopJetSelection(1, -1, id_cms3, h_ak8jets));
    hist_tt_cms3.reset(new AndHists(ctx, "CMS3"));

    // 1% misstag
    TopJetId id_cms1 = AndId<TopJet>(Type2TopTag(105,220, Type2TopTag::MassType::groomed), Tau32(0.67));
    tt_cms1.reset(new NTopJetSelection(1, -1, id_cms1, h_ak8jets));
    hist_tt_cms1.reset(new AndHists(ctx, "CMS1"));

    // 0.3% misstag
    TopJetId id_cms03 = AndId<TopJet>(Type2TopTag(105,220, Type2TopTag::MassType::groomed), Tau32(0.54));
    tt_cms03.reset(new NTopJetSelection(1, -1, id_cms03, h_ak8jets));
    hist_tt_cms03.reset(new AndHists(ctx, "CMS03"));

    // 0.1% misstag
    TopJetId id_cms01 = AndId<TopJet>(Type2TopTag(105,220, Type2TopTag::MassType::groomed), Tau32(0.46));
    tt_cms01.reset(new NTopJetSelection(1, -1, id_cms01, h_ak8jets));
    hist_tt_cms01.reset(new AndHists(ctx, "CMS01"));

  }

  bool HOTVRPerformanceModule::process(Event & event) {
    
    // cl_muo->process(event);
    // if(!sel_muo->passes(event)) return false;
    hist_sel_muo->fill(event);

    jec_hotvr->process(event);

    cl_hotvr->process(event);
    cl_ak8->process(event);

    if(tt_hotvr->passes(event))
      {
	hist_tt_hotvr->fill(event);
      }

    if(tt_cms3->passes(event))
      {
	hist_tt_cms3->fill(event);
      }

    if(tt_cms1->passes(event))
      {
	hist_tt_cms1->fill(event);
      }

    if(tt_cms03->passes(event))
      {
	hist_tt_cms03->fill(event);
      }

    if(tt_cms01->passes(event))
      {
	hist_tt_cms01->fill(event);
      }
 
    return true;
  }

 UHH2_REGISTER_ANALYSIS_MODULE(HOTVRPerformanceModule)

}
