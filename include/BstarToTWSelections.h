#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"

#include "UHH2/common/include/TopJetIds.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesis.h"

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

  /**
   * This class selects events by requiring ST > ST_min
   */
  class STSelection: public uhh2::Selection {
  public:
    STSelection(double st_min);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double m_st_min;
  };

  class Chi2Selection: public uhh2::Selection {
  public:
    Chi2Selection(uhh2::Context &ctx, std::string label, double chi2_max);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double m_chi2_max;
    uhh2::Event::Handle<std::vector<BstarToTWHypothesis>> h_hyp;    

  };
  
  class RecoMassSelection: public uhh2::Selection {
  public:
    RecoMassSelection(uhh2::Context &ctx, double m_min_, std::string label);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double m_min;
    uhh2::Event::Handle<std::vector<BstarToTWHypothesis>> h_hyp;
  };

  class JetDeltaPhiSelection: public uhh2::Selection {
  public:
    JetDeltaPhiSelection(uhh2::Context &ctx, double delta_phi_min, const boost::optional<uhh2::Event::Handle<std::vector<Jet> > > jet_collection = boost::none);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double m_delta_phi_min;
    boost::optional<uhh2::Event::Handle<std::vector<Jet> > > h_jets;
    uhh2::Event::Handle<FlavorParticle> h_primlep;
  };

  class MassCutSelection: public uhh2::Selection {
  public:
    MassCutSelection(double m_min_, double m_max_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double m_min, m_max;
  };


  class NGenJetSelection: public uhh2::Selection {
  public:
    NGenJetSelection(unsigned int n_min_, unsigned int n_max_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    unsigned int n_min, n_max;
  };
}
