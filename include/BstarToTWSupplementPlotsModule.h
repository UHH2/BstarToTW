#include <iostream>
#include <memory>
#include <chrono>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/NSelections.h"

#include "UHH2/HOTVR/include/HOTVRIds.h"

#include "UHH2/BstarToTW/include/AndHists.h"


using namespace std;
using namespace uhh2;

namespace uhh2 {

  /*
    Run an additional set of plots for supplemental studies
    * b tag categories without top tag applied
    * plots for 0b 1+b categories
    * top tag without mass window
    * ...
  //*/
  class BstarToTWSupplementPlotsModule: public AnalysisModule {
  public:
    
    explicit BstarToTWSupplementPlotsModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  
    // b tag categories
    std::unique_ptr<Selection> sel_0btag, sel_1plusbtag, sel_1btag, sel_2btag;
    std::unique_ptr<AndHists> hist_0btag, hist_1plusbtag, hist_1btag, hist_2btag; // btag only hists
    // top tag without mass window
    std::unique_ptr<Selection> sel_1toptag;
    std::unique_ptr<AndHists> hist_1toptag;

  };
}
