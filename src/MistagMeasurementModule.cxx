#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/LumiSelection.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/TTbarReconstruction.h"

#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/BstarToTW/include/BstarToTWModules.h"
#include "UHH2/BstarToTW/include/BstarToTWReconstruction.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisHists.h"
#include "UHH2/BstarToTW/include/GenCleaningModules.h"

#include "UHH2/HOTVR/include/HOTVRHists.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/TopTagPerformanceHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class MistagMeasurementModule: public AnalysisModule {
  public:

    explicit MistagMeasurementModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:

    std::unique_ptr<ObjectCleaner> cl_objects;
    std::unique_ptr<AnalysisModule> cl_topjet;
    std::unique_ptr<Selection> sel_1top;
    std::unique_ptr<Hists> hist_mistag;

  };

  MistagMeasurementModule::MistagMeasurementModule(Context & ctx) {

    // double deltaPhi_min = M_PI/2;  // minimum delta phi between muon and top
    double top_pt_min = 200.0;
    double top_eta_max = 2.5;
    TopJetId id_topjet =  PtEtaCut(top_pt_min, top_eta_max);
    TopJetId id_toptag = AndId<TopJet>(HOTVRTopTag(0.8, 140., 220., 50.), Tau32Groomed(0.56));

    cl_objects.reset(new ObjectCleaner(ctx));
    // cl_objects->switch_jet_lepton_cleaning(false);
    cl_topjet.reset(new TopJetCleaner(ctx, id_topjet));
    sel_1top.reset(new NTopJetSelection(1, -1));

    hist_mistag.reset(new MistagHists(ctx, "Mistag", id_toptag));

  }

  bool MistagMeasurementModule::process(Event & event) {
    // event cleaning
    if(!(cl_objects->process(event))) return false;
    cl_topjet->process(event);
    // require >1 top
    if (!sel_1top->passes(event)) return false;
    // measure mistag
    hist_mistag->fill(event);
    
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(MistagMeasurementModule)
}
