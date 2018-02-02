#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"
#include "UHH2/BstarToTW/include/BstarToTWHypothesis.h"
#include "UHH2/common/include/PDFWeights.h" 


namespace uhh2 {

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class BstarToTWPDFHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    BstarToTWPDFHists(uhh2::Context & ctx, const std::string & dirname, bool use_ntupleweights_, bool use_pdf_weights_ = false);

    virtual void fill(const uhh2::Event & ev) override;
    std::string histo_names[100];

  protected:
    uhh2::Event::Handle<std::vector<BstarToTWHypothesis>> h_hyps;
    std::string m_discriminator_name;
    bool use_ntupleweights, use_pdf_weights;
    bool is_mc, is_LO, take_ntupleweights;
    TString m_oname;
    std::unique_ptr<PDFWeights> m_pdfweights;


    virtual ~BstarToTWPDFHists();
};

}
