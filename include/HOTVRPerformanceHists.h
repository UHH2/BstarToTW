#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/TopTagIndexer.h"

#include <vector>

namespace uhh2 {

  /**  \brief Example class for booking and filling histograms
   * 
   * NOTE: This class uses the 'hist' method to retrieve histograms.
   * This requires a string lookup and is therefore slow if you have
   * many histograms. Therefore, it is recommended to use histogram
   * pointers as member data instead, like in 'common/include/ElectronHists.h'.
   */
  class HOTVRPerformanceHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    HOTVRPerformanceHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & event) override;
    virtual ~HOTVRPerformanceHists();
  protected:
    TH1F *DeltaR_Top_HotvrTopjets;
    TH2F *jet_area_vs_jet_pt;
    TH2F *pt_reco_over_pt_top_vs_pt_reco, *pt_reco_over_pt_top_vs_pt_top, *pt_reco_over_pt_top_vs_npv;

    uhh2::Event::Handle<BstarToTWGen> h_BstarToTWGen;
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    uhh2::Event::Handle<std::vector<TopJet>> h_AK8Jets;
  };
   
}
