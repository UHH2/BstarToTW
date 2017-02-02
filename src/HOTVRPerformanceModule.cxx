#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/common/include/EventHists.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWHists.h"
#include "UHH2/BstarToTW/include/BstarToTWCleaningModules.h"
#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
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
    // std::unique_ptr<CommonModules> common; 

    // GenParticle Interpreter
    std::unique_ptr<AnalysisModule> BstarToTWgenprod;

    // Cleaners
    std::unique_ptr<AnalysisModule> cl_top, cl_gentop;

    // Selections
    std::unique_ptr<Selection> sel_n_top, sel_pt_top;
    std::unique_ptr<Selection> sel_n_gentop;
    std::unique_ptr<Selection> sel_deltaR;
 
    // Histograms
    std::unique_ptr<Hists> h_nocuts, h_cl_top, h_sel_n_top;
    std::unique_ptr<Hists> h_nocuts_event, h_sel_n_top_event;
    std::unique_ptr<Hists> h_performance;
    std::unique_ptr<Hists> h_efficiency;
  };


  HOTVRPerformanceModule::HOTVRPerformanceModule(Context & ctx) {

    // GenParticle Interpreter
    BstarToTWgenprod.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));

    // Cleaners
    cl_top.reset(new HOTVRTopCleaner()); 
    cl_gentop.reset(new HOTVRGenTopCleaner());

    // Selections
    sel_n_top.reset(new NHotvrTopSelection(1)); // == 1 Top
    sel_n_gentop.reset(new NHotvrGenTopSelection(1)); // == 1 genTop

    // Histograms
    h_nocuts.reset(new BstarToTWHists(ctx, "No_Cuts"));
    h_nocuts_event.reset(new EventHists(ctx, "No_Cuts_Event"));
    h_cl_top.reset(new BstarToTWHists(ctx, "Topjet_Cleaner"));
    h_sel_n_top.reset(new BstarToTWHists(ctx, "NTop_Cut"));
    h_sel_n_top_event.reset(new EventHists(ctx, "NTop_Cut_Event"));
    h_performance.reset(new HOTVRPerformanceHists(ctx, "HOTVR_Performance"));
    h_efficiency.reset(new EfficiencyHists(ctx, "Efficiency"));
  }


  bool HOTVRPerformanceModule::process(Event & event) {
    BstarToTWgenprod->process(event);
    h_nocuts->fill(event);
    h_nocuts_event->fill(event);

    // applying cleaners
    cl_top->process(event);
    h_cl_top->fill(event);
    cl_gentop->process(event);

    h_efficiency->fill(event);
 
    // event selection

    if(!sel_n_top->passes(event)) return false;
    h_sel_n_top->fill(event);
    h_sel_n_top_event->fill(event);

    if(!sel_n_gentop->passes(event)) return false;
    h_performance->fill(event);

    // Done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(HOTVRPerformanceModule)

}
