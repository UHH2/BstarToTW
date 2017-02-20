#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/NSelections.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/TopTagIndexer.h"
#include "UHH2/BstarToTW/include/IndexCleaner.h"
#include "UHH2/BstarToTW/include/HOTVRHists.h"
#include "UHH2/BstarToTW/include/HOTVRPerformanceHists.h"
#include "UHH2/BstarToTW/include/EfficiencyHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWPreSelectionModule: public AnalysisModule {
  public:
    
    explicit BstarToTWPreSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  
    // additional modules
    std::unique_ptr<AnalysisModule> TopTagIndProd;

    // index cleaner
    std::unique_ptr<AnalysisModule> cl_HotvrTop, cl_EtaTop, cl_PtTop, cl_Tau32Top;

    // selection
    MuonId id_Muo;

    std::unique_ptr<Selection> sel_NTop;
    std::unique_ptr<Selection> sel_NMuo;

    // histograms
    std::unique_ptr<Hists> h_ntop_hotvr, h_ntop_performance;

    std::string version;
  };

  BstarToTWPreSelectionModule::BstarToTWPreSelectionModule(Context & ctx) {

    // additional modules
    TopTagIndProd.reset(new TopTagIndexerProducer(ctx, "TopTagIndexer"));

    // index cleaner
    double top_pt_min         = 200.;
    // double top_m_min          = 140.;
    // double top_m_max          = 220.;
    double top_eta_max        = 2.4;
    unsigned int top_nsub_min = 3;
    double top_fpt_max        = 0.8;
    double top_mpair_min      = 50.;
    double top_t32_max        = 0.69;

    cl_HotvrTop.reset(new HotvrTopIndexCleaner(ctx, top_pt_min, top_eta_max, top_nsub_min, top_fpt_max , top_mpair_min, top_t32_max));

    // selections
    id_Muo = AndId<Muon>(MuonIDTight(), PtEtaCut(130.0, 2.4),MuonIso(0.15));

    sel_NTop.reset(new NTopIndexSelection(ctx, 1, 1));
    sel_NMuo.reset(new NMuonSelection(1, -1, id_Muo));

    // histogramms
    h_ntop_hotvr.reset(new HOTVRHists(ctx, "NTopCut_HOTVR"));
    h_ntop_performance.reset(new HOTVRPerformanceHists(ctx, "NTopCut_Performance"));

    version = ctx.get("dataset_version");
  }

  bool BstarToTWPreSelectionModule::process(Event & event) {
    // initialize additional modules
    TopTagIndProd->process(event);
 

    // Selection
    if(!sel_NMuo->passes(event)) return false;
    if(!cl_HotvrTop->process(event)) return false;
    if(!sel_NTop->passes(event)) return false;
    h_ntop_hotvr->fill(event);
    if (version.find("BstarToTW") != string::npos) 
      {
	h_ntop_performance->fill(event);
      }

    // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWPreSelectionModule)

}
