#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/NSelections.h"

#include "UHH2/BstarToTW/include/BstarToTWSelections.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWPreSelectionModule: public AnalysisModule {
  public:
    
    explicit BstarToTWPreSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  

    MuonId id_muo;

    // Cleaner
    std::unique_ptr<AnalysisModule> cl_muo;

    // Selections
    std::unique_ptr<Selection> sel_NMuo;
    std::unique_ptr<Selection> sel_NTop;
  };

  BstarToTWPreSelectionModule::BstarToTWPreSelectionModule(Context & ctx) {

    // variables for selections
    // TODO: Implement triggers IsoMu24_OR_IsoTkMu24 (HLT_IsoMu24_v* || HLT_IsoTkMu24_v*)
    // Trigger is now in Selection. Maybe move this and MET cut to Preselection?
    double muo_pt_min  = 50.; 
    double muo_eta_max = 2.4;
    double muo_iso_max = 0.15;
    double top_pt_min  = 200.;
    double top_eta_max = 2.4;

    id_muo = AndId<Muon>(MuonIDTight(), PtEtaCut(muo_pt_min, muo_eta_max),MuonIso(muo_iso_max));

    // cleaner
    cl_muo.reset(new MuonCleaner(id_muo));

    // selections

    sel_NMuo.reset(new NMuonSelection(1, -1, id_muo));
    sel_NTop.reset(new NHotvrSelection(1, -1, top_pt_min, top_eta_max));
  }

  bool BstarToTWPreSelectionModule::process(Event & event) {

    cl_muo->process(event);

    // Muon selection
    if(!sel_NMuo->passes(event)) return false;

    // TopJet Selection
    if(!sel_NTop->passes(event)) return false;

    // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWPreSelectionModule)

}
