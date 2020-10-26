#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/NSelections.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWGenSelections.h"
#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/BstarToTWGenHists.h"

#include "UHH2/HOTVR/include/HOTVRGenIds.h"
#include "UHH2/HOTVR/include/HOTVRHists.h"
#include "UHH2/HOTVR/include/HOTVRIds.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class BstarToTWGenModule: public AnalysisModule {
  public:
    
    explicit BstarToTWGenModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    // Additional Modules
    // std::unique_ptr<AnalysisModule> printer;
    // std::unique_ptr<AnalysisModule> BstarToTWgenprod;
  
    std::unique_ptr<AndHists> hist_all;


  };


  BstarToTWGenModule::BstarToTWGenModule(Context & ctx){

    // Additional Modules
    // printer.reset(new GenParticlesPrinter(ctx));  
    // BstarToTWgenprod.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));

  }


  bool BstarToTWGenModule::process(Event & event) {
    // BstarToTWgenprod->process(event);
    // hist_all->fill(event);

    for (const auto & gp : *event.genparticles)
      {
	if (gp.status() == 21 && abs(gp.pdgId()) != 21 && abs(gp.pdgId()) != 5)
	  {
	    cout << "PDG_ID: " << gp.pdgId() << endl;
	    cout << "Daughters: " << gp.daughter1() << ", "<< gp.daughter2() << endl;
	  }
      }

    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWGenModule)

}
