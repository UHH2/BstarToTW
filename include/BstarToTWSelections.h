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

  /* MassHotvrTopSelection
   *
   * Selection of Topjets within Masswindow m_min < m < m_max
   */
  class MHotvrTopSelection: public uhh2::Selection {
  public:
    MHotvrTopSelection(double m_min_ = 140, double m_max_ = 220);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    double m_min;
    double m_max;

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
