#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/Muon.h"
#include "UHH2/common/include/TopJetIds.h"

#include "UHH2/common/include/JetIds.h"

#include "UHH2/BstarToTW/include/HOTVRIds.h"

#include <vector>

namespace uhh2 {

  /**  \brief Example class for booking and filling histograms
   * 
   * NOTE: This class uses the 'hist' method to retrieve histograms.
   * This requires a string lookup and is therefore slow if you have
   * many histograms. Therefore, it is recommended to use histogram
   * pointers as member data instead, like in 'common/include/ElectronHists.h'.
   */
  class HOTVRHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    HOTVRHists(uhh2::Context & ctx, const std::string & dirname, const boost::optional<TopJetId> &topjetid = boost::none);

    virtual void fill(const uhh2::Event & event) override;
    virtual ~HOTVRHists();
  protected:
    // HOTVR
    TH1F *N_HotvrTopjets, *Pt_HotvrTopjets, *Eta_HotvrTopjets, *M_HotvrTopjets;
    TH1F *NSub_HotvrTopjets, *Fpt_HotvrTopjets, *Mpair_HotvrTopjets, *Tau32_HotvrTopjets;
    TH1F *DeltaR_L_HotvrTopjets, *DeltaPhi_L_HotvrTopjets;
    // b-jets
    TH1F *N_bjets_loose, *N_bjets_medium, *N_bjets_tight;

    boost::optional<TopJetId> m_topjetid;
  };
   
}
