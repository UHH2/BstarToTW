#pragma once

#include "UHH2/core/include/Hists.h"

namespace uhh2 {

  /**  \brief Example class for booking and filling histograms
   * 
   * NOTE: This class uses the 'hist' method to retrieve histograms.
   * This requires a string lookup and is therefore slow if you have
   * many histograms. Therefore, it is recommended to use histogram
   * pointers as member data instead, like in 'common/include/ElectronHists.h'.
   */
  class CutflowHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    CutflowHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & event);
    void fill(const uhh2::Event & event, int i);
    virtual ~CutflowHists();
  protected:
    TH1F *cutflow;
      };
   
}
