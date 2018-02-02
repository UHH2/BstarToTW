#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/common/include/LuminosityHists.h"



  /**  \brief Example class for booking and filling histograms
   * 
   * NOTE: This class uses the 'hist' method to retrieve histograms.
   * This requires a string lookup and is therefore slow if you have
   * many histograms. Therefore, it is recommended to use histogram
   * pointers as member data instead, like in 'common/include/ElectronHists.h'.
   */

namespace uhh2 {
  class BstarToTWHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    BstarToTWHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~BstarToTWHists();
  protected:
    std::vector<run_lumi> upper_binborders;
    double lumi_per_bin;

    TH1F *MET, *HT_lep, *HT_jet, *ST, *rho;
    TH2F *LumiBlock_vs_NPV;

  };
   
}
