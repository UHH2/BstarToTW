#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/MuonHists.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/TopTagIndexer.h"
#include "UHH2/BstarToTW/include/IndexCleaner.h"
#include "UHH2/BstarToTW/include/HOTVRHists.h"
#include "UHH2/BstarToTW/include/HOTVRPerformanceHists.h"
#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/EfficiencyHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class HOTVRPerformanceModule: public AnalysisModule {
  public:
    
    explicit HOTVRPerformanceModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  
    // additional modules
    std::unique_ptr<AnalysisModule> BstarToTWGenProd, TopTagIndProd;

    // index cleaner
    std::unique_ptr<AnalysisModule> cl_PtTop, cl_MTop, cl_EtaTop, cl_NSubTop, cl_FptTop, cl_MpairTop, cl_Tau32Top;

    // selection
    MuonId id_Muo;

    std::unique_ptr<Selection> sel_NMuo;
    std::unique_ptr<Selection> sel_NTop;

    // histograms
    std::unique_ptr<AndHists> hist_nocuts, hist_nmuo, hist_pt, hist_eta, hist_nsub, hist_fpt, hist_mpair, hist_tau32, hist_ntop_pt, hist_ntop;
    std::unique_ptr<Hists> hist_efficiency;
  };

  HOTVRPerformanceModule::HOTVRPerformanceModule(Context & ctx) {

    // additional modules
    BstarToTWGenProd.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    TopTagIndProd.reset(new TopTagIndexerProducer(ctx, "TopTagIndexer"));

    // index cleaner
    double top_pt_min         = 200.;
    double top_eta_max        = 2.4;
    unsigned int top_nsub_min = 3;
    double top_fpt_max        = 0.8;
    double top_mpair_min      = 50.;
    double top_tau32_max      = 0.69;

    cl_PtTop.reset(new PtTopIndexCleaner(ctx, top_pt_min));
    cl_EtaTop.reset(new EtaTopIndexCleaner(ctx, top_eta_max));
    cl_NSubTop.reset(new NSubTopIndexCleaner(ctx, top_nsub_min));
    cl_FptTop.reset(new FptTopIndexCleaner(ctx, top_fpt_max));
    cl_MpairTop.reset(new MpairTopIndexCleaner(ctx, top_mpair_min));
    cl_Tau32Top.reset(new Tau32TopIndexCleaner(ctx, top_tau32_max));

    // selections
    id_Muo = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.4),MuonIso(0.15));

    sel_NMuo.reset(new NMuonSelection(1, -1, id_Muo));
    sel_NTop.reset(new NTopIndexSelection(ctx, 1, 1));

    // histogramms
    hist_nocuts.reset(new AndHists());
    hist_nocuts->add_hist(new HOTVRHists(ctx, "NoCuts_hotvr"));
    hist_nocuts->add_hist(new HOTVRPerformanceHists(ctx, "NoCuts_performance"));
    hist_nocuts->add_hist(new MuonHists(ctx, "NoCuts_muons"));

    hist_nmuo.reset(new AndHists());
    hist_nmuo->add_hist(new HOTVRHists(ctx, "NMuo_hotvr"));
    hist_nmuo->add_hist(new HOTVRPerformanceHists(ctx, "NMuo_performance"));

    hist_pt.reset(new AndHists());
    hist_pt->add_hist(new HOTVRHists(ctx, "PtTop_hotvr"));
    hist_pt->add_hist(new HOTVRPerformanceHists(ctx, "PtTop_performance"));

    hist_eta.reset(new AndHists());
    hist_eta->add_hist(new HOTVRHists(ctx, "EtaTop_hotvr"));
    hist_eta->add_hist(new HOTVRPerformanceHists(ctx, "EtaTop_performance"));

    hist_nsub.reset(new AndHists());
    hist_nsub->add_hist(new HOTVRHists(ctx, "NSubTop_hotvr"));
    hist_nsub->add_hist(new HOTVRPerformanceHists(ctx, "NSubTop_performance"));

    hist_fpt.reset(new AndHists());
    hist_fpt->add_hist(new HOTVRHists(ctx, "FptTop_hotvr"));
    hist_fpt->add_hist(new HOTVRPerformanceHists(ctx, "FptTop_performance"));

    hist_mpair.reset(new AndHists());
    hist_mpair->add_hist(new HOTVRHists(ctx, "MpairTop_hotvr"));
    hist_mpair->add_hist(new HOTVRPerformanceHists(ctx, "MpairTop_performance"));

    hist_tau32.reset(new AndHists());
    hist_tau32->add_hist(new HOTVRHists(ctx, "Tau32Top_hotvr"));
    hist_tau32->add_hist(new HOTVRPerformanceHists(ctx, "Tau32Top_performance"));

    hist_ntop_pt.reset(new AndHists());
    hist_ntop_pt->add_hist(new HOTVRHists(ctx, "NTop_pt_hotvr"));
    hist_ntop_pt->add_hist(new HOTVRPerformanceHists(ctx, "NTop_pt_performance"));

    hist_ntop.reset(new AndHists());
    hist_ntop->add_hist(new HOTVRHists(ctx, "NTop_hotvr"));
    hist_ntop->add_hist(new HOTVRPerformanceHists(ctx, "NTop_performance"));
    hist_ntop->add_hist(new MuonHists(ctx, "NTop_muons"));

    hist_efficiency.reset(new EfficiencyHists(ctx, "Efficiency"));

  }

  bool HOTVRPerformanceModule::process(Event & event) {
    // initialize additional modules
    BstarToTWGenProd->process(event);
    TopTagIndProd->process(event);
    hist_efficiency->fill(event);
    hist_nocuts->fill(event);

    // index cleaner
    if(!sel_NMuo->passes(event)) return false;
    hist_nmuo->fill(event);

    if(!cl_PtTop->process(event)) return false;
    hist_pt->fill(event);

    if(!cl_EtaTop->process(event)) return false;
    hist_eta->fill(event);

    if(!cl_NSubTop->process(event)) return false;
    hist_nsub->fill(event);

    if(!cl_FptTop->process(event)) return false;
    hist_fpt->fill(event);

    if(!cl_MpairTop->process(event)) return false;
    hist_mpair->fill(event);

    if (sel_NTop->passes(event)) hist_ntop_pt->fill(event);

    if(!cl_Tau32Top->process(event)) return false;
    hist_tau32->fill(event);

    if(!sel_NTop->passes(event)) return false;
    hist_ntop->fill(event);

    // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(HOTVRPerformanceModule)

}
