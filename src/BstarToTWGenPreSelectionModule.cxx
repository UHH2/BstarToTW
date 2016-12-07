
#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/common/include/PrintingModules.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWGenHists.h"
#include "UHH2/BstarToTW/include/BstarToTWGenSelection.h"
#include "UHH2/BstarToTW/include/GenJetCluster.h"
#include "UHH2/BstarToTW/include/BstarToTWGenJetHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {  

  class BstarToTWGenPreSelectionModule: public AnalysisModule 
  {
  public:
    
    explicit BstarToTWGenPreSelectionModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    std::unique_ptr<AnalysisModule> printer;
    std::unique_ptr<AnalysisModule> GenProd, GenjetProd;

    // Selections
    std::unique_ptr<Selection> semilep_sel, tophad_sel, njet_sel, cmstoptag_sel, hotvrtoptag_sel;

    // Histograms
    std::unique_ptr<Hists> h_nocuts, h_semilep, h_tophad;
    std::unique_ptr<Hists> h_njet, h_ncmstop, h_nhotvrtop; 
    std::unique_ptr<Hists> h_njet_th, h_ncmstop_th, h_nhotvrtop_th; 
    std::unique_ptr<Hists> h_nocuts_genjets, h_semilep_genjets, h_tophad_genjets;
    std::unique_ptr<Hists> h_njet_genjets, h_ncmstop_genjets, h_nhotvrtop_genjets;
    std::unique_ptr<Hists> h_njet_genjets_th, h_ncmstop_genjets_th, h_nhotvrtop_genjets_th;
  };

  /* ----------------------------------------------------------------
   * BstarToTWGenPreSelectionModule - Analysis module for running
   * preselection on generator level.
   */
  BstarToTWGenPreSelectionModule::BstarToTWGenPreSelectionModule(Context & ctx)
  {
    printer.reset(new GenParticlesPrinter(ctx));  

    GenProd.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    GenjetProd.reset(new GenJetProducer(ctx, "BstarToTWgenjet"));
    
    semilep_sel.reset(new IsSemiLeptonic(ctx));
    tophad_sel.reset(new IsTopHad(ctx));

    njet_sel.reset(new NJet(ctx, 3, 25));
    cmstoptag_sel.reset(new NCMSTopJet(ctx, 1, 1));
    hotvrtoptag_sel.reset(new NHOTVRTopJet(ctx, 1, 1));

    h_nocuts.reset(new BstarToTWGenHists(ctx, "NoCuts"));
    h_semilep.reset(new BstarToTWGenHists(ctx, "SemiLep"));
    h_tophad.reset(new BstarToTWGenHists(ctx, "TopHad"));
    h_njet.reset(new BstarToTWGenHists(ctx, "NJet"));
    h_ncmstop.reset(new BstarToTWGenHists(ctx, "NCMSTopTagged"));
    h_nhotvrtop.reset(new BstarToTWGenHists(ctx, "NHOTVRTopTagged"));
    h_njet_th.reset(new BstarToTWGenHists(ctx, "NJet_th"));
    h_ncmstop_th.reset(new BstarToTWGenHists(ctx, "NCMSTopTagged_th"));
    h_nhotvrtop_th.reset(new BstarToTWGenHists(ctx, "NHOTVRTopTagged_th"));

    h_nocuts_genjets.reset(new BstarToTWGenJetHists(ctx, "NoCuts_genjets"));
    h_semilep_genjets.reset(new BstarToTWGenJetHists(ctx, "SemiLep_genjets"));
    h_tophad_genjets.reset(new BstarToTWGenJetHists(ctx, "TopHad_genjets"));
    h_njet_genjets.reset(new BstarToTWGenJetHists(ctx, "NJet_genjets"));
    h_ncmstop_genjets.reset(new BstarToTWGenJetHists(ctx, "NCMSTopTagged_genjets"));
    h_nhotvrtop_genjets.reset(new BstarToTWGenJetHists(ctx, "NHOTVRTopTagged_genjets"));
    h_njet_genjets_th.reset(new BstarToTWGenJetHists(ctx, "NJet_genjets_th"));
    h_ncmstop_genjets_th.reset(new BstarToTWGenJetHists(ctx, "NCMSTopTagged_genjets_th"));
    h_nhotvrtop_genjets_th.reset(new BstarToTWGenJetHists(ctx, "NHOTVRTopTagged_genjets_th"));
 
  }

  bool BstarToTWGenPreSelectionModule::process(Event & event) 
  {
    //printer->process(event);
    // --------------------------------------------------------------
    // GenParticle and Genjet Producer
    GenProd->process(event);
    GenjetProd->process(event);

    h_nocuts->fill(event);
    h_nocuts_genjets->fill(event);

    // --------------------------------------------------------------
    // Pre Selection

    // Select
    if(semilep_sel->passes(event))
      {	
	h_semilep->fill(event);
	h_semilep_genjets->fill(event);
	if(tophad_sel->passes(event))
	  {
	    h_tophad->fill(event);
	    h_tophad_genjets->fill(event);
	    if(njet_sel->passes(event))
	      {
		h_njet_th->fill(event);
		h_njet_genjets_th->fill(event);
	      }
	    if(cmstoptag_sel->passes(event))
	      {
		h_ncmstop_th->fill(event);
		h_ncmstop_genjets_th->fill(event);
	      }
	    if(hotvrtoptag_sel->passes(event))
	      {
		h_nhotvrtop_th->fill(event);
		h_nhotvrtop_genjets_th->fill(event);
	      }
	  }
      }

    // full signal
    if(njet_sel->passes(event))
      {
	h_njet_genjets->fill(event);
	h_njet_genjets->fill(event);
      }
    if(cmstoptag_sel->passes(event))
      {
	h_ncmstop->fill(event);
	h_ncmstop_genjets->fill(event);
      }
    if(hotvrtoptag_sel->passes(event))
      {
	h_nhotvrtop->fill(event);
	h_nhotvrtop_genjets->fill(event);
      }

    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWGenPreSelectionModule)

}
