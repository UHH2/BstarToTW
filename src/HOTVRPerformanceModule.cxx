#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/NSelections.h"

#include "UHH2/BstarToTW/include/HOTVRIds.h"
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

  };

  HOTVRPerformanceModule::HOTVRPerformanceModule(Context & ctx) {

    std::string ak8jets_name = "slimmedJetsAK8_SoftDrop";
    h_ak8jets = ctx.get_handle<vector<TopJet>>(ak8jets_name);

    // jec_hotvr.reset(new HOTVRJetCorrector(ctx));

    // Kinematic variables
    double muo_pt_min = 50.0;
    double muo_eta_max = 2.4;
    double muo_iso_max = 0.15;
   
    double top_pt_min = 200.0;
    double top_eta_max = 2.5;
 
  }

  bool HOTVRPerformanceModule::process(Event & event) {
    
  
    return true;
  }

 UHH2_REGISTER_ANALYSIS_MODULE(HOTVRPerformanceModule)

}
