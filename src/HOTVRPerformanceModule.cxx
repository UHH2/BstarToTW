#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/core/include/TopJet.h"

#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/JetHists.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWGenSelections.h"
#include "UHH2/BstarToTW/include/HOTVRIds.h"
#include "UHH2/BstarToTW/include/HOTVRHists.h"
#include "UHH2/BstarToTW/include/EfficiencyHists.h"

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
    std::unique_ptr<AnalysisModule> BstarToTWgenprod; // gen particle interpreter
    std::unique_ptr<Selection> sel_semilep; // select only semi-leptonic events on generator level.
    Event::Handle<vector<TopJet>> h_ak8jets; // handle for ak8_SoftDrop collection
    std::unique_ptr<Hists> hist_hotvr, hist_ak8;

    // HOTVR
    TopJetId id_hotvrtoptag;
    std::unique_ptr<Selection> sel_hotvrtoptag;
    std::unique_ptr<Hists> hist_hotvrtoptag;
    std::unique_ptr<TopJetHists> hist_cmstoptag;

    // CMSTopTag ungroomed
    TopJetId id_cmstoptag_ungroomed;
    std::unique_ptr<Selection> sel_cmstoptag_ungroomed;
    std::unique_ptr<Hists> hist_cmstoptag_ungroomed;

    // CMSTopTag groomed
    TopJetId id_cmstoptag_groomed;
    std::unique_ptr<Selection> sel_cmstoptag_groomed;
    std::unique_ptr<Hists> hist_cmstoptag_groomed;

  };

  HOTVRPerformanceModule::HOTVRPerformanceModule(Context & ctx) {
 
    BstarToTWgenprod.reset(new BstarToTWGenProducer(ctx, "BstarToTWgen"));
    sel_semilep.reset(new SemiLepSelection(ctx));
    
    std::string ak8jets_name = "slimmedJetsAK8_SoftDrop";

    h_ak8jets = ctx.get_handle<vector<TopJet>>(ak8jets_name);

    hist_hotvr.reset(new HOTVRHists(ctx, "HOTVR_Jets"));
    hist_ak8.reset(new TopJetHists(ctx, "AK8_Jets", 4, ak8jets_name));
    hist_cmstoptag.reset(new TopJetHists(ctx, "CMSTopTag_Jets", 4, ak8jets_name));

    // HOTVR
    id_hotvrtoptag = HOTVRTopTag(); // default settings for HOTVRTopTag, see HOTVRId.h
    sel_hotvrtoptag.reset(new NTopJetSelection(1,1, id_hotvrtoptag));
    hist_hotvrtoptag.reset(new EfficiencyHists(ctx, "HOTVRTopTag_Efficiencies"));

    // CMSTopTag ungroomed
    id_cmstoptag_ungroomed = CMSTopTag(); // default settings for CMSTopTag, see TopJetIds.h
    sel_cmstoptag_ungroomed.reset(new NTopJetSelection(1, 1, id_cmstoptag_ungroomed, h_ak8jets));
    hist_cmstoptag_ungroomed.reset(new EfficiencyHists(ctx, "CMSTopTag_ungroomed_Efficiencies",  h_ak8jets));
    hist_cmstoptag->set_TopJetId(id_cmstoptag_ungroomed);

    // CMSTopTag groomed
    id_cmstoptag_groomed = CMSTopTag(CMSTopTag::MassType::groomed); // as above but uses groomed mass
    sel_cmstoptag_groomed.reset(new NTopJetSelection(1, 1, id_cmstoptag_groomed, h_ak8jets));
    hist_cmstoptag_groomed.reset(new EfficiencyHists(ctx, "CMSTopTag_groomed_Efficiencies",  h_ak8jets));
  }

  bool HOTVRPerformanceModule::process(Event & event) {
    BstarToTWgenprod->process(event);
    // vector<TopJet> ak8jets = event.get(h_ak8jets);
    
    // cout << "N Jets = " << ak8jets.size() << endl;
    // for (TopJet topjet : ak8jets)
    //   {
    // 	vector<Jet> subjets = topjet.subjets();
    // 	cout << "N subjets = " << subjets.size() << endl;
    // 	if(subjets.size() < 3) continue;
	
    // 	float mjet = topjet.v4().M();
    // 	cout << "mjet = " << mjet << endl;

    // 	double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
    // 	double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
    // 	double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
    // 	double mmin = min(min(m12, m13), m23);
    // 	cout << "mmin = " << mmin << endl;
    //   }

    // Only semi-leptonic channel
    if(!sel_semilep->passes(event)) return false;
    hist_hotvr->fill(event);
    hist_ak8->fill(event);

    if(sel_hotvrtoptag->passes(event)) hist_hotvrtoptag->fill(event);

    if(sel_cmstoptag_ungroomed->passes(event)) hist_cmstoptag_ungroomed->fill(event);

    if(sel_cmstoptag_groomed->passes(event)) hist_cmstoptag_groomed->fill(event);


    // done
    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(HOTVRPerformanceModule)

}
