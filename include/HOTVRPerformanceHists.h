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
    HOTVRPerformanceHists(uhh2::Context & ctx, const std::string & dirname, const boost::optional<Event::Handle<std::vector<TopJet> > > &topjetcollection = boost::none);

    virtual void fill(const uhh2::Event & event) override;
    virtual ~HOTVRPerformanceHists();
  protected:
    TH1F *DeltaR_Top_HotvrTopjet, *EffPt_Top_HotvrTopjet;
    TH2F *EffPt_Top_HotvrTopjet_vs_pt_top, *EffPt_Top_HotvrTopjet_vs_npv;

    uhh2::Event::Handle<BstarToTWGen> h_BstarToTWGen;
    boost::optional<Event::Handle<std::vector<TopJet>>> h_topjetcollection;
  };


  class GenHOTVRPerformanceHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    GenHOTVRPerformanceHists(uhh2::Context & ctx, const std::string & dirname, const boost::optional<Event::Handle<std::vector<GenTopJet> > > &topjetcollection = boost::none);

    virtual void fill(const uhh2::Event & event) override;
    virtual ~GenHOTVRPerformanceHists();
  protected:
    TH1F *DeltaR_Top_HotvrTopjet, *EffPt_Top_HotvrTopjet;
    TH2F *EffPt_Top_HotvrTopjet_vs_pt_top, *EffPt_Top_HotvrTopjet_vs_npv;

    uhh2::Event::Handle<BstarToTWGen> h_BstarToTWGen;
    boost::optional<Event::Handle<std::vector<GenTopJet>>> h_topjetcollection;
  };
  
}
