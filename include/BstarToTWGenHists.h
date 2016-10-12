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

    virtual void fill(const uhh2::Event & ev) override;

protected:
    TH1F *M_bstar, *pt_bstar, *pt_top, *eta_top, *phi_top, *M_top, *pt_W1, *eta_W1, *phi_W1, *M_W1, *pt_W2, *eta_W2, *phi_W2, *M_W2, *pt_b, *eta_b, *phi_b, *M_b, *pt_ele, *eta_ele, *phi_ele, *M_ele, *pt_muo, *eta_muo, *phi_muo, *M_muo, *pt_elenu, *eta_elenu, *phi_elenu, *pt_muonu, *eta_muonu, *phi_muonu, *pt_q1, *eta_q1, *phi_q1, *M_q1, *pt_q2, *eta_q2, *phi_q2, *M_q2;
    uhh2::Event::Handle<BstarToTWGen> h_BstarToTWgen;
};
