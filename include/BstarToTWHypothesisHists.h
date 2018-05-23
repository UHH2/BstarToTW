#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/Utils.h"

#include "UHH2/BstarToTW/include/BstarToTWHypothesis.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"


/** \brief Common histograms for reconstruction hypotheses
 *
 * hyps_name is the name of the reconstruction hypothesis collection, for instance "HighMassReconstruction"
 * discriminator_name is the name of the discriminator used to choose the best reconstruction hypothesis, for instance "Chi2"
 */
class BstarToTWHypothesisHists: public uhh2::Hists {
public:
    BstarToTWHypothesisHists(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name);

    virtual void fill(const uhh2::Event & ev) override;

protected:
    TH1F *Discriminator, *Discriminator_2, *Discriminator_3;
    TH1F *Bstar_reco_M, *Bstar_reco_M_rebin, *Bstar_reco_M_unbinned, *Bstar_reco_Pt;
    TH1F *W_reco_M, *W_reco_Pt;
    TH1F *Top_reco_M, *Top_reco_Pt;
    TH1F *DeltaR_top_W, *DeltaPhi_top_W, *DeltaPt_top_W, *DeltaPt_top_W_over_pt, *PtTop_over_PtW;
    TH2F *Discriminator_vs_M_bstar;


    uhh2::Event::Handle<std::vector<BstarToTWHypothesis>> h_hyps;
    std::string m_name;
    std::string m_discriminator_name;
    
};
