#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/GenJetCluster.h"

#include "fastjet/PseudoJet.hh"

#include <vector>

namespace uhh2 {
 
  class IsSemiLeptonic: public uhh2::Selection {
  public:
    IsSemiLeptonic(Context & ctx);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    Event::Handle<BstarToTWGen> h_BstarToTWGen;
  };
   
  class IsTopHad: public uhh2::Selection {
  public:
    IsTopHad(Context & ctx);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    Event::Handle<BstarToTWGen> h_BstarToTWGen;
  };

  class NJet: public uhh2::Selection {
  public:
    NJet(Context & ctx, unsigned int n_min, double pt_min);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    Event::Handle<GenJetCluster> h_genjetcluster;
    unsigned int _n_min;
    double _pt_min;
  };

  class NCMSTopJet: public uhh2::Selection {
  public:
    NCMSTopJet(Context & ctx, unsigned int n_min, unsigned int n_max);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    Event::Handle<GenJetCluster> h_genjetcluster;
    unsigned int _n_min;
    unsigned int _n_max;
  };

  class NHOTVRTopJet: public uhh2::Selection {
  public:
    NHOTVRTopJet(Context & ctx, unsigned int n_min, unsigned int n_max);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    Event::Handle<GenJetCluster> h_genjetcluster;
    unsigned int _n_min;
    unsigned int _n_max;
  };
}
