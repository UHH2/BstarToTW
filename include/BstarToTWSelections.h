#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"

namespace uhh2 {

  /** 
   * This class selects event by requiring there to be at least n_min
   * and maximal n_max topjets, with pt > pt_min and eta < eta_max
   */
  class NHotvrSelection: public uhh2::Selection {
  public:
    NHotvrSelection(unsigned int n_min_, unsigned int n_max_, double pt_min_, double eta_max_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    unsigned int n_min, n_max;
    double pt_min, eta_max;
  };

  /**
   * This class selects events by requiring MET > met_min
   */
  class METSelection: public uhh2::Selection {
  public:
    METSelection(double met_min_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double met_min;
  };

}
