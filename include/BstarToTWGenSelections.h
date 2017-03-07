#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"

namespace uhh2 {

  class SemiLepSelection: public uhh2::Selection {
  public:
    SemiLepSelection(uhh2::Context &ctx);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    uhh2::Event::Handle<BstarToTWGen> h_BstarToTWGen;
  };


}
