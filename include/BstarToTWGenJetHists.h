#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"
#include "UHH2/BstarToTW/include/GenJetCluster.h"
#include "fastjet/PseudoJet.hh"
#include "vector"

class BstarToTWGenJetHists: public uhh2::Hists {
public:
    BstarToTWGenJetHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & event) override;


protected:
    TH1F *N_ak4Jets, *N_ak8Jets, *N_topJets_CMS, *N_topJets_HOTVR;
    TH1F *M_top_CMS, *M_top_HOTVR;

    uhh2::Event::Handle<GenJetCluster> h_genjetcluster;

};
