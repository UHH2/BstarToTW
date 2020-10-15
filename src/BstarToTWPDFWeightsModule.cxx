#include <iostream>
#include <memory>
#include <chrono>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/PDFWeights.h" 
#include "UHH2/common/include/Utils.h" 


using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BstarToTWPDFWeightsModule: public AnalysisModule {
  public:
    
    explicit BstarToTWPDFWeightsModule(Context & ctx);
    ~BstarToTWPDFWeightsModule();
    virtual bool process(Event & event) override;

  private:  
    std::vector<double> m_sumofweights;
    int m_N_tot;
    std::unique_ptr<PDFWeights> m_pdfweights;
    TString m_oname, m_pdfname, m_path_out;
  };

  BstarToTWPDFWeightsModule::BstarToTWPDFWeightsModule(Context & ctx) {
    m_oname = ctx.get("dataset_version");
    m_path_out = ctx.get("PDFWeightsPath");

    Year year = extract_year(ctx);

    if (year == Year::is2016v2 || year == Year::is2016v3)
      {
	if (m_oname.Contains("BstarToTW3000") || m_oname.Contains("BstarToTW2") || (m_oname.Contains("BstarToTW1") && !m_oname.Contains("BstarToTW11") && !m_oname.Contains("BstarToTW10")) ) // change pdfname for b*<=3TeV samples
	  m_pdfname = "MMHT2014lo68cl";
	else if (m_oname.Contains("BstarToTW")) // change pdfname for b* > 3TeV
	  m_pdfname = "PDF4LHC15_nnlo_30_pdfas";
      }
    else if(year == Year::is2017v1 || year == Year::is2017v2)
      {
	if (m_oname.Contains("BstarToTW")) // change pdfname for b*
	  m_pdfname = "PDF4LHC15_nnlo_30_pdfas";
      }
    else if(year == Year::is2018)
      {
	if (m_oname.Contains("BstarToTW")) // change pdfname for b*
	  m_pdfname = "PDF4LHC15_nnlo_30_pdfas";
      }

    cout << "Writing weights for :" << m_pdfname << endl;
    m_pdfweights.reset(new PDFWeights(m_pdfname)); 

    vector<double> vect(m_pdfweights->GetNWeights(),0);
    m_sumofweights = vect;
    m_N_tot = 0;

  }

  BstarToTWPDFWeightsModule::~BstarToTWPDFWeightsModule() {
    cout << "N kept events:  " << m_N_tot << endl;

    TString path_out = m_path_out;
    path_out += "/" + m_oname + "_" + m_pdfname + "_weights.txt";
    cout << "Writing to file: " << path_out << endl;

    std::ofstream outfile;
    outfile.open(((std::string)path_out).c_str());
    outfile << m_N_tot <<std::endl;
    for(unsigned int i=0; i< m_pdfweights->GetNWeights(); ++i){
      outfile<< m_sumofweights[i] << " ";
    }
    outfile.close();
    m_pdfweights.reset(nullptr);
  }

  bool BstarToTWPDFWeightsModule::process(Event & event) {
    
    std::vector<double> weights = m_pdfweights->GetWeightList(event);
    if(!weights.size()) 
      return false;
    m_N_tot++;

    for(unsigned int i=0; i< m_pdfweights->GetNWeights(); ++i)
      m_sumofweights[i]+=weights[i];

    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BstarToTWPDFWeightsModule)

}
