#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"

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

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~HOTVRPerformanceHists();
  protected:
    
    TH1F *h_pt_reco_over_pt_top, *h_pt_reco_over_pt_gen;
    TH1F *h_deltaR_reco_top, *h_deltaR_reco_gen, *h_deltaR_gen_top;

    TH2F *h_pt_reco_over_pt_top_vs_pt, *h_pt_reco_over_pt_gen_vs_pt, *h_pt_reco_over_pt_top_vs_npv, *h_pt_reco_over_pt_gen_vs_npv;

    TH1F *h_matched_pt_reco, *h_matched_M_reco, *h_unmatched_pt_reco, *h_unmatched_M_reco;

    TH2F *h_jet_area_vs_jet_pt;
    
    TH1F *hotvr_counter, *cms_counter;

    uhh2::Event::Handle<BstarToTWGen> h_BstarToTWgen;
    uhh2::Event::Handle<std::vector<TopJet>> h_AK8Jets;
  };
   
}
