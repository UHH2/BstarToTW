#pragma once

#include "UHH2/core/include/Hists.h"

#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/PrimaryLepton.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/ObjectIdUtils.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesis.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"

#include "UHH2/HOTVR/include/HadronicTop.h"

#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>


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
    uhh2::Event::Handle<FlavorParticle> h_primlep;
    JetId btag_loose = CSVBTag(CSVBTag::WP_LOOSE);

    TH1F *MET, *HT_lep, *HT_jet, *ST, *rho, *deltaPhi_blep, *deltaPhi_btop, *deltaPhi_hotvr_ak4, *deltaPhi_lep_ak4, *deltaPhi_hotvr_lak4, *deltaPhi_lep_lak4;
    TH2F *LumiBlock_vs_NPV;
  };

  class BstarToTWBackgroundHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    BstarToTWBackgroundHists(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const TString & path);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~BstarToTWBackgroundHists();
  protected:
    double *m_pa, *m_pb;
    /* double m_p0, m_p1; */
    uhh2::Event::Handle<std::vector<BstarToTWHypothesis>> h_hyps;
    TH1F *Bstar_reco_M, *Bstar_reco_M_rebin, *Bstar_reco_M_unbinned;
  };

  class BstarToTWAnalysisHists: public uhh2::Hists {
  public:
    BstarToTWAnalysisHists(uhh2::Context & ctx, const std::string & dirname);
    virtual void fill(const uhh2::Event &event) override;
    virtual ~BstarToTWAnalysisHists();
  protected:
    uhh2::Event::Handle<FlavorParticle> h_primlep;
    TH1F *deltaPhi_blep, *deltaPhi_btop, *deltaPhi_hotvr_ak4, *deltaPhi_lep_ak4, *deltaPhi_hotvr_lak4, *deltaPhi_lep_lak4;

  };

  class TopMatchHists: public uhh2::Hists {
  public:
    TopMatchHists(uhh2::Context &ctx, const std::string &dirname, const std::string & hyps_name);
    virtual void fill(const uhh2::Event &event) override;
    virtual ~TopMatchHists();
  protected:
    uhh2::Event::Handle<std::vector<GenTop> > h_tophad;
    uhh2::Event::Handle<std::vector<BstarToTWHypothesis> > h_hyps;
    TH1F *Bstar_reco_M, *HOTVR_pt, *HOTVR_eta;
  };
   
}
