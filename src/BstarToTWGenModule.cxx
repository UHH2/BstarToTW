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
#include "UHH2/BstarToTW/include/HOTVRGenIds.h"
#include "UHH2/BstarToTW/include/HOTVRHists.h"
#include "UHH2/BstarToTW/include/HOTVRPerformanceHists.h"
#include "UHH2/BstarToTW/include/HOTVRIds.h"
#include "UHH2/BstarToTW/include/GenCleaningModules.h"
#include "UHH2/BstarToTW/include/GenNSelections.h"
#include "UHH2/BstarToTW/include/EfficiencyHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class BstarToTWGenModule: public AnalysisModule {
  public:
    
    explicit BstarToTWGenModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    // Additional Modules
    std::unique_ptr<AnalysisModule> printer;
    std::unique_ptr<AnalysisModule> BstarToTWgenprod;

    // Cleaner
    std::unique_ptr<AnalysisModule> cl_genhotvrtoptag;
    std::unique_ptr<AnalysisModule> cl_hotvrtoptag;
    std::unique_ptr<AnalysisModule> cl_muo;

    // Selections
    std::unique_ptr<Selection> sel_muonchannel, sel_ngentopjet, sel_ntopjet;
    std::unique_ptr<Selection> sel_muo;

    GenTopJetId id_genhotvrtoptag;

    TopJetId id_hotvrtoptag;

    MuonId id_muon;

    // Hists
    std::unique_ptr<AndHists> hist_all, hist_muonchannel, hist_genhotvr, hist_hotvr, hist_muo;


  };


  BstarToTWGenModule::BstarToTWGenModule(Context & ctx){

    // Additional Modules
    printer.reset(new GenParticlesPrinter(ctx));  
    BstarToTWgenprod.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));

    // Variables
    double top_pt_min = 200.0;
    double top_eta_max = 2.5;

    id_genhotvrtoptag = AndId<GenTopJet>(PtEtaCut(top_pt_min, top_eta_max), GenHOTVRTopTag(), GenTopMuonCleaner(ctx));
    id_hotvrtoptag = AndId<TopJet>(PtEtaCut(top_pt_min, top_eta_max), HOTVRTopTag(), GenLepOverlap(ctx)); // default settings for HOTVRTopTag, see HOTVRId.h
    id_muon = AndId<Muon>(MuonIDTight(), PtEtaCut(50.0, 2.4), MuonIso(0.15)) ;

    // Cleaner
    cl_genhotvrtoptag.reset(new GenTopJetCleaner(ctx, id_genhotvrtoptag));
    cl_hotvrtoptag.reset(new TopJetCleaner(ctx, id_hotvrtoptag));
    cl_muo.reset(new MuonCleaner(id_muon));

    // Selections
    sel_muonchannel.reset(new MuonChannelSelection(ctx));
    sel_ngentopjet.reset(new NGenTopJetSelection(1,1));
    sel_ntopjet.reset(new NTopJetSelection(1, 1));
    sel_muo.reset(new NMuonSelection(1,-1));
    
    // Hists
    hist_all.reset(new AndHists(ctx, "All"));
    hist_all->add_hist(new BstarToTWGenHists(ctx, "All_Gen")); 
    hist_all->add_hist(new HOTVRGenHists(ctx, "All_HOTVRGen")); 
   
    hist_muonchannel.reset(new AndHists(ctx, "MuonChannel"));
    hist_muonchannel->add_hist(new BstarToTWGenHists(ctx, "MuonChannel_Gen"));
    hist_muonchannel->add_hist(new HOTVRGenHists(ctx, "MuonChannel_HOTVRGen")); 
    hist_muonchannel->add_hist(new EfficiencyHists(ctx, "MuonChannel_Efficiencies")); 

    hist_genhotvr.reset(new AndHists(ctx, "GenHOTVR"));
    hist_genhotvr->add_hist(new BstarToTWGenHists(ctx, "GenHOTVR_Gen"));
    hist_genhotvr->add_hist(new HOTVRGenHists(ctx, "GenHOTVR_HOTVRGen")); 
    hist_genhotvr->add_hist(new GenEfficiencyHists(ctx, "GenHOTVR_Efficiencies")); 
    hist_genhotvr->add_hist(new GenHOTVRPerformanceHists(ctx, "GenHOTVR_Performance")); 

    hist_hotvr.reset(new AndHists(ctx, "HOTVR"));
    hist_hotvr->add_hist(new BstarToTWGenHists(ctx, "HOTVR_Gen"));
    hist_hotvr->add_hist(new EfficiencyHists(ctx, "HOTVR_Efficiencies"));
    hist_hotvr->add_hist(new HOTVRPerformanceHists(ctx, "HOTVR_Performance"));

    hist_muo.reset(new AndHists(ctx, "Muon"));
    hist_muo->add_hist(new BstarToTWGenHists(ctx, "Muon_Gen"));
    hist_muo->add_hist(new HOTVRGenHists(ctx, "Muon_HOTVRGen")); 
    hist_muo->add_hist(new EfficiencyHists(ctx, "Muon_Efficiencies")); 
  }


  bool BstarToTWGenModule::process(Event & event) {
    BstarToTWgenprod->process(event);
    hist_all->fill(event);

    cl_muo->process(event);
    if(sel_muo->passes(event)) hist_muo->fill(event);

    if(!sel_muonchannel->passes(event)) return false;
    hist_muonchannel->fill(event);
    
    cl_genhotvrtoptag->process(event);
    if(sel_ngentopjet->passes(event)) hist_genhotvr->fill(event);

    cl_hotvrtoptag->process(event);
    if(sel_ntopjet->passes(event)) hist_hotvr->fill(event);

    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWGenModule)

}
