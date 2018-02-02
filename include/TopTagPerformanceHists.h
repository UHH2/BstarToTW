#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"

#include "UHH2/common/include/TopJetIds.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include <vector>

namespace uhh2 {

  class TopTagPerformanceHists: public uhh2::Hists {
  public:

    TopTagPerformanceHists(uhh2::Context & ctx, const std::string & dirname, bool is_qcd, TopJetId id, const boost::optional<Event::Handle<std::vector<TopJet> > > &topjetcollection = boost::none);

    virtual void fill(const uhh2::Event & event) override;
    virtual ~TopTagPerformanceHists();

  protected:

    TH1F *N_evt, *N_mat, *pt_top, *pt_matched, *pt_mismatched, *pt_gen;
    TH2F *dR_max;
    bool b_is_qcd;
    TopJetId m_id;
    boost::optional<Event::Handle<std::vector<TopJet>>> h_topjetcollection;

  };

}
