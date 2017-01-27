#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"

namespace uhh2 {
    
  /* NHotvrTopSelection
   *
   * Select events that have at least n_min TopJets following the
   * definition of HOTVR (1606.04961)
   */
  class NHotvrTopSelection: public uhh2::Selection {
  public:
    NHotvrTopSelection(unsigned int n);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    unsigned int n;

  };

  /* ptHotvrTopSelection
   *
   * Select events with pt of leading top > pt_min
   */
  class ptHotvrTopSelection: public uhh2::Selection {
  public:
    ptHotvrTopSelection(double pt_min_);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    double pt_min;

  };

  /* NHotvrGenTopSelection
   *
   * Select events that have at least n_min TopJets following the
   * definition of HOTVR (1606.04961)
   */
  class NHotvrGenTopSelection: public uhh2::Selection {
  public:
    NHotvrGenTopSelection(unsigned int n);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    unsigned int n;

  };

}
