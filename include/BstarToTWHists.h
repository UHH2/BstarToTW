#pragma once

#include "UHH2/core/include/Hists.h"

namespace uhh2examples {

  /**  \brief Example class for booking and filling histograms
   * 
   * NOTE: This class uses the 'hist' method to retrieve histograms.
   * This requires a string lookup and is therefore slow if you have
   * many histograms. Therefore, it is recommended to use histogram
   * pointers as member data instead, like in 'common/include/ElectronHists.h'.
   */
  class BstarToTWHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    BstarToTWHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~BstarToTWHists();
  protected:
    // jet histograms
    TH1F *N_jets, *eta_jet1, *eta_jet2, *eta_jet3, *pt_jet1, *pt_jet2, *pt_jet3;

    // lepton histograms
    TH1F *N_mu, *eta_mu, *pt_mu, *reliso_mu;

    // misc
    TH1F *N_pv;
      };
   
}
