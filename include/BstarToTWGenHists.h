#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/BstarToTW/include/BstarToTWGen.h"

/** \brief Histograms for BstarToTW quantities on generator (parton) level
 * 
 * BstarToTWgen container has to be filled before calling this histogram class
 */
class BstarToTWGenHists: public uhh2::Hists {
public:
    BstarToTWGenHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & event) override;



protected:
    TH1F *M_bstar, *pt_bstar;
    TH1F *pt_top, *eta_top, *phi_top, *M_top, *M_tophad, *M_toplep;
    TH1F *pt_W1, *eta_W1, *phi_W1, *M_W1;
    TH1F *pt_W2, *eta_W2, *phi_W2, *M_W2;
    TH1F *pt_b, *eta_b, *phi_b, *M_b;
    TH1F *pt_ele, *eta_ele, *phi_ele;
    TH1F *pt_muo, *eta_muo, *phi_muo;
    TH1F *cosThetastar_bstar, *cosThetastar_t, *deltaRMax;

    uhh2::Event::Handle<BstarToTWGen> h_BstarToTWgen;
};
