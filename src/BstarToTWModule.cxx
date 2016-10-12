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

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/BstarToTW/include/BstarToTWHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  /** \brief 
   * 
   * This is the central class which calls other AnalysisModules, Hists or Selection classes.
   * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
   */
  class BstarToTWModule: public AnalysisModule {
  public:
    
    explicit BstarToTWModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    
    std::unique_ptr<CommonModules> common;
    
    std::unique_ptr<JetCleaner> jetcleaner;
   
    // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
    // to avoid memory leaks.
    std::unique_ptr<Selection> njet_sel, nele_sel;
    
    // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
    std::unique_ptr<Hists> h_nocuts, h_njet, h_ele, h_nele;

    std::unique_ptr<AnalysisModule> genparticleprinter;
    MuonId MuId;
    ElectronId EleId;
    JetId JetId1;
  };


  BstarToTWModule::BstarToTWModule(Context & ctx){
    
    // If needed, access the configuration of the module here, e.g.:
    string testvalue = ctx.get("TestKey", "<not set>");
    cout << "TestKey in the configuration was: " << testvalue << endl;
    
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }
    
    // IDs
    EleId = AndId<Electron>(ElectronID_Spring15_25ns_medium, PtEtaCut(30.0, 2.4));
    MuId = AndId<Muon>(MuonIDTight(), PtEtaCut(30.0, 2.4),MuonIso(0.15));
    JetId1 = PtEtaCut(30.0,2.4);

    // Common Modules
    common.reset(new CommonModules());
    // common->switch_jetlepcleaner(true);
    common->set_electron_id(EleId);
    common->set_muon_id(MuId);
    common->set_jet_id(JetId1);
    common->init(ctx);  
    
    // Printer
    genparticleprinter.reset(new GenParticlesPrinter(ctx));
    
    // Selections
    njet_sel.reset(new NJetSelection(2)); // see common/include/NSelections.h
    nele_sel.reset(new NElectronSelection(1)); // see common/include/NSelections.h

    // Histograms
    h_nocuts.reset(new BstarToTWHists(ctx, "NoCuts"));
    h_njet.reset(new BstarToTWHists(ctx, "Njet"));
    h_nele.reset(new BstarToTWHists(ctx, "Nele"));
    h_ele.reset(new ElectronHists(ctx, "ele_nocuts"));
  }


  bool BstarToTWModule::process(Event & event) {
    // This is the main procedure, called for each event. Typically,
    // do some pre-processing by calling the modules' process method
    // of the modules constructed in the constructor (1).
    // Then, test whether the event passes some selection and -- if yes --
    // use it to fill the histograms (2).
    // Finally, decide whether or not to keep the event in the output (3);
    // this is controlled by the return value of this method: If it
    // returns true, the event is kept; if it returns false, the event
    // is thrown away.
    
 
    // 1. run all modules other modules.
    if(!common->process(event)) return false;
    
    // 2. test selections and fill histograms
    h_ele->fill(event);
    h_nocuts->fill(event);
    
    if(!njet_sel->passes(event)){
      return false;
    } 
    h_njet->fill(event);    
    
    if(!nele_sel->passes(event)){
      return false;
    }
    h_nele->fill(event);
    return true;
  }

  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the BstarToTWModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWModule)

}
