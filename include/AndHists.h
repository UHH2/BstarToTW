#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/TopJetIds.h"

#include <vector>
#include <TH1F.h>

namespace uhh2 {
  class AndHists: public Hists {
  public:
    AndHists(uhh2::Context &ctx, const std::string & dirname, const boost::optional<TopJetId> &id_topjet = boost::none);

    virtual void fill(const uhh2::Event & event);
    virtual ~AndHists();
    void add_hist(uhh2::Hists *hist);

  protected:
    std::vector<uhh2::Hists*> hists_vector;

    TH1F *nevt;

  };

}
