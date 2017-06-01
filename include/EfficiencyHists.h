#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/TopJet.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"

#include <vector>

namespace uhh2 {

  class EfficiencyHists: public uhh2::Hists {
  public:
    EfficiencyHists(uhh2::Context & ctx, const std::string & dirname, const boost::optional<Event::Handle<std::vector<TopJet>>> &topjetcollection = boost::none);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~EfficiencyHists();
  protected:
    
    TH1F *All, *Matched;


    uhh2::Event::Handle<BstarToTWGen> h_bstargen;
    boost::optional<Event::Handle<std::vector<TopJet>>> h_topjetcollection;
  };


  class GenEfficiencyHists: public uhh2::Hists {
  public:
    GenEfficiencyHists(uhh2::Context & ctx, const std::string & dirname, const boost::optional<Event::Handle<std::vector<GenTopJet>>> &gentopjetcollection = boost::none);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~GenEfficiencyHists();
  protected:
    
    TH1F *All, *Matched;


    uhh2::Event::Handle<BstarToTWGen> h_bstargen;
    boost::optional<Event::Handle<std::vector<GenTopJet>>> h_gentopjetcollection;
  };
}
