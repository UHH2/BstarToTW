#pragma once
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

namespace uhh2 {
    
  /* HOTVRTopCleaner
   *
   * Select events that have at least n_min TopJets following the
   * definition of HOTVR (1606.04961)
   */
  class HOTVRTopCleaner: public uhh2::AnalysisModule {
  public:
    explicit HOTVRTopCleaner(unsigned int nsub_min = 3, double fpt_max = 0.8, double m_min = 140, double m_max = 220, double mpairwise_min = 50);
    virtual bool process(uhh2::Event & event) override;
  private:
    unsigned int nsub_min;
    double fpt_max;
    double m_min, m_max;
    double mpairwise_min;

  };

  /* HOTVRGenTopCleaner
   *
   * Select events that have at least n_min TopJets following the
   * definition of HOTVR (1606.04961)
   */
  class HOTVRGenTopCleaner: public uhh2::AnalysisModule {
  public:
    explicit HOTVRGenTopCleaner(unsigned int nsub_min = 3, double fpt_max = 0.8, double m_min = 140, double m_max = 220, double mpairwise_min = 50);
    virtual bool process(uhh2::Event & event) override;
  private:
    unsigned int nsub_min;
    double fpt_max;
    double m_min, m_max;
    double mpairwise_min;

  };

}
