#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include <vector>

namespace uhh2 {
  class AndHists {
  public:
    AndHists(uhh2::Context &ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & event);
    virtual ~AndHists();
    void add_hist(uhh2::Hists *hist);

  protected:
    std::vector<uhh2::Hists*> hists_vector;

  };

}
