#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/TopTagIndexer.h"
#include "UHH2/BstarToTW/include/IndexCleaner.h"
#include "UHH2/BstarToTW/include/HOTVRPerformanceHists.h"
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
    std::unique_ptr<Selection> sel_NTop;

    // histograms
    std::unique_ptr<Hists> h_nocuts, h_pt, h_eta, h_m, h_nsub, h_fpt, h_mpair, h_tau32, h_n;
    std::unique_ptr<Hists> h_efficiency;
  };

  HOTVRPerformanceModule::HOTVRPerformanceModule(Context & ctx) {

    // additional modules
    BstarToTWGenProd.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    TopTagIndProd.reset(new TopTagIndexerProducer(ctx, "TopTagIndexer"));

    // index cleaner
    double top_pt_min         = 200.;
    double top_m_min          = 140.;
    double top_m_max          = 220.;
    double top_eta_max        = 2.4;
    unsigned int top_nsub_min = 3;
    double top_fpt_max        = 0.8;
    double top_mpair_min      = 50.;
    double top_tau32_max      = 0.69;

    cl_PtTop.reset(new PtTopIndexCleaner(ctx, top_pt_min));
    cl_MTop.reset(new MTopIndexCleaner(ctx, top_m_min, top_m_max));
    cl_EtaTop.reset(new EtaTopIndexCleaner(ctx, top_eta_max));
    cl_NSubTop.reset(new NSubTopIndexCleaner(ctx, top_nsub_min));
    cl_FptTop.reset(new FptTopIndexCleaner(ctx, top_fpt_max));
    cl_MpairTop.reset(new MpairTopIndexCleaner(ctx, top_mpair_min));
    cl_Tau32Top.reset(new MpairTopIndexCleaner(ctx, top_tau32_max));

    // selections
    sel_NTop.reset(new NTopIndexSelection(ctx, 1));

    // histogramms
    h_nocuts.reset(new HOTVRPerformanceHists(ctx, "NoCuts"));
    h_pt.reset(new HOTVRPerformanceHists(ctx, "PtCut"));
    h_eta.reset(new HOTVRPerformanceHists(ctx, "EtaCut"));
    h_m.reset(new HOTVRPerformanceHists(ctx, "MassCut"));
    h_nsub.reset(new HOTVRPerformanceHists(ctx, "NSubjetCut"));
    h_fpt.reset(new HOTVRPerformanceHists(ctx, "FptCut"));
    h_mpair.reset(new HOTVRPerformanceHists(ctx, "MpairCut"));
    h_tau32.reset(new HOTVRPerformanceHists(ctx, "Tau32Cut"));
    h_n.reset(new HOTVRPerformanceHists(ctx, "NTopCut"));

    h_efficiency.reset(new EfficiencyHists(ctx, "Efficiency"));

  }

  bool HOTVRPerformanceModule::process(Event & event) {
    // initialize additional modules
    BstarToTWGenProd->process(event);
    TopTagIndProd->process(event);
    h_nocuts->fill(event);
    h_efficiency->fill(event);

    // index cleaner
    if(!cl_PtTop->process(event)) return false;
    h_pt->fill(event);
    if(!cl_EtaTop->process(event)) return false;
    h_eta->fill(event);
    if(!cl_MTop->process(event)) return false;
    h_m->fill(event);
    if(!cl_NSubTop->process(event)) return false;
    h_nsub->fill(event);
    if(!cl_FptTop->process(event)) return false;
    h_fpt->fill(event);
    if(!cl_MpairTop->process(event)) return false;
    h_mpair->fill(event);
    if(!cl_Tau32Top->process(event)) return false;
    h_tau32->fill(event);
    if(!sel_NTop->passes(event)) return false;
    h_n->fill(event);

    // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(HOTVRPerformanceModule)

}
