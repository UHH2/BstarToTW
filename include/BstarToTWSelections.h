#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"

#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesis.h"

#include <TRandom.h>

namespace uhh2 {

  /** 
   * This class selects event by requiring there to be at least n_min
   * and maximal n_max topjets, with pt > pt_min and eta < eta_max
   */

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
    STSelection(uhh2::Context &ctx, double st_min);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    uhh2::Event::Handle<double> h_ht;
    uhh2::Event::Handle<FlavorParticle> h_primlep;
    double m_st_min;
  };

  class HTSelection: public uhh2::Selection {
  public:
    HTSelection(uhh2::Context &ctx, double ht_min);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    uhh2::Event::Handle<double> h_ht;
    double m_ht_min;
  };

  class Chi2Selection: public uhh2::Selection {
  public:
    Chi2Selection(uhh2::Context &ctx, std::string label, double chi2_max, const std::string disc_name = "Chi2");
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double m_chi2_max;
    uhh2::Event::Handle<std::vector<BstarToTWHypothesis>> h_hyp;    
    const std::string m_disc_name;
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
  class LeadingJetDeltaRSelection: public uhh2::Selection {
  public:
    LeadingJetDeltaRSelection(uhh2::Context &ctx, double delta_R_min, const boost::optional<uhh2::Event::Handle<std::vector<Jet> > > jet_collection = boost::none);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double m_delta_R_min;
    boost::optional<uhh2::Event::Handle<std::vector<Jet> > > h_jets;
    uhh2::Event::Handle<FlavorParticle> h_primlep;
  };

  class JetDeltaRSelection: public uhh2::Selection {
  public:
    JetDeltaRSelection(uhh2::Context &ctx, double delta_R_min, const boost::optional<uhh2::Event::Handle<std::vector<Jet> > > jet_collection = boost::none);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double m_delta_R_min;
    boost::optional<uhh2::Event::Handle<std::vector<Jet> > > h_jets;
    uhh2::Event::Handle<FlavorParticle> h_primlep;
  };

  class TopJetDeltaRSelection: public uhh2::Selection {
  public:
    TopJetDeltaRSelection(uhh2::Context &ctx, double delta_R_max, const boost::optional<uhh2::Event::Handle<std::vector<TopJet> > > topjet_collection = boost::none);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double m_delta_R_max;
    bool do_R_calc;
    boost::optional<uhh2::Event::Handle<std::vector<TopJet> > > h_topjets;
    uhh2::Event::Handle<std::vector<Jet> > h_bjets;
  };

  class METDeltaPhiSelection: public uhh2::Selection {
  public:
    METDeltaPhiSelection(uhh2::Context &ctx, double delta_phi_max);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double m_delta_phi_max;
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

  class BstarToTWTriggerSelection : public uhh2::Selection {
  public:
    BstarToTWTriggerSelection(uhh2::Context &ctx);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    Year year;
    bool is_muo, is_ele;
    std::unique_ptr<uhh2::Selection> trig_isomu24, trig_isotkmu24, trig_isomu27;
    std::unique_ptr<uhh2::Selection> trig_ele27, trig_ele32, trig_ele35;
    std::unique_ptr<uhh2::Selection> trig_photon175, trig_photon200;
  };

  class BadHCALSelection: public uhh2::Selection {
  public:
    BadHCALSelection(uhh2::Context &ctx, long int seed = 123456789);
    virtual bool passes(const uhh2::Event &event) override;

  private:
    TRandom *m_rng;
    long int m_seed;
    Year year;
    int m_runnumber = 319077;
    double m_lumi_ratio = 0.64844705699; // (Run 319077(17.370008 pb-1) + Run C + Run D) / all 2018

    double m_interval_eta = -1.3;
    double m_interval_phi_low = -1.57;
    double m_interval_phi_high = -0.87;

  };
}
