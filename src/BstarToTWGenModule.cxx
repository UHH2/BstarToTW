#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/BstarToTW/include/BstarToTWGenHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class BstarToTWGenModule: public AnalysisModule {
  public:
    
    explicit BstarToTWGenModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    std::unique_ptr<AnalysisModule> printer;
    std::unique_ptr<AnalysisModule> BstarToTWgenprod;
    std::unique_ptr<Hists> h_BstarToTWgenhists;
  };


  BstarToTWGenModule::BstarToTWGenModule(Context & ctx){


    printer.reset(new GenParticlesPrinter(ctx));  

    BstarToTWgenprod.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    h_BstarToTWgenhists.reset(new BstarToTWGenHists(ctx, "BstarToTWgenhists"));
  }


  bool BstarToTWGenModule::process(Event & event) {
    BstarToTWgenprod->process(event);

    h_BstarToTWgenhists->fill(event);

    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWGenModule)

}
