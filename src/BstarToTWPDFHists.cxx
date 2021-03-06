#include "UHH2/BstarToTW/include/BstarToTWPDFHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include <math.h>
#include <sstream>

#include "TH1F.h"
#include "TH2D.h"
#include <iostream>

using namespace std;
using namespace uhh2;


BstarToTWPDFHists::BstarToTWPDFHists(Context & ctx, const string & dirname, bool use_ntupleweights_, bool use_pdf_weights_): Hists(ctx, dirname), use_ntupleweights(use_ntupleweights_), use_pdf_weights(use_pdf_weights_){  

  is_mc = ctx.get("dataset_type") == "MC";
  m_oname = ctx.get("dataset_version");
  TString m_pdfname;
  TString weightpath = ctx.get("PDFWeightsPath") + "/" + m_oname;
  is_LO = m_oname.Contains("Diboson") || m_oname.Contains("DYJets") || m_oname.Contains("QCD"); // only non-madgraph

  // if (is_LO)
  m_pdfname = "NNPDF30_lo_as_0130";

  
  Year year = extract_year(ctx);
  // take_ntupleweights =  use_ntupleweights && (!is_LO || m_oname.Contains("DYJets")) && !(m_oname.Contains("ST_tW") && m_oname.Contains("2016v3"));

  if (year == Year::is2016v2 || year == Year::is2016v3)
    {
      take_ntupleweights = !(m_oname.Contains("QCD") || m_oname.Contains("ST_tW") || m_oname.Contains("Diboson") ||m_oname.Contains("BstarToTW") || m_oname.Contains("DYJets"));
      if (m_oname.Contains("BstarToTW3000") || m_oname.Contains("BstarToTW2") || (m_oname.Contains("BstarToTW1") && !m_oname.Contains("BstarToTW11") && !m_oname.Contains("BstarToTW10")) ) // change pdfname for 1200 <= b* <=3TeV samples
	m_pdfname = "MMHT2014lo68cl";
      else if (m_oname.Contains("BstarToTW")) // change pdfname for b* > 3TeV
	m_pdfname = "PDF4LHC15_nnlo_30_pdfas";
      // else if (m_oname.Contains("ST_tW"))
      // 	m_pdfname = "NNPDF30_nlo_as_0118";
	      
    }
  else if(year == Year::is2017v1 || year == Year::is2017v2)
    {
      take_ntupleweights = !(m_oname.Contains("QCD")|| m_oname.Contains("BstarToTW") || m_oname.Contains("Diboson") || m_oname.Contains("DYJets"));
      if (m_oname.Contains("BstarToTW")) // change pdfname for b*
	m_pdfname = "PDF4LHC15_nnlo_30_pdfas";
    }
  else if(year == Year::is2018)
    {
    take_ntupleweights = !(m_oname.Contains("QCD") || m_oname.Contains("Diboson") || m_oname.Contains("BstarToTW") || m_oname.Contains("Diboson") || m_oname.Contains("DYJets"));
    if (m_oname.Contains("BstarToTW")) // change pdfname for b*
      m_pdfname = "PDF4LHC15_nnlo_30_pdfas";
    }

  //For Mbstar reconstruction
  h_hyps = ctx.get_handle<std::vector<BstarToTWHypothesis>>("tW_reco");
  m_discriminator_name ="closest_nu"; 

  //if(!is_LO && !m_oname.Contains("SingleTop")) m_pdfname = "PDF4LHC15_nlo_mc"; 

  //take ntupleweights if
  //1) use_ntupleweights = true and the sample has ntupleweights stored

  //take weights from txt files if
  //1) the sample is LO (and doesn't have ntupleweights for that reason) (this assumption is protected by a runtime_error later)
  //2) the sample is NLO and yet doesn't have ntupleweights 


  cout << "For this sample '" << m_oname << "' is_LO is set to " << is_LO << endl;
  cout << "Are ntupleweights taken for this sample?: " << take_ntupleweights << endl;

  if(is_mc && !take_ntupleweights){
    if (m_oname.Contains("BstarToTW"))
      m_pdfweights.reset(new PDFWeights(m_pdfname, weightpath));
    else
      m_pdfweights.reset(new PDFWeights(m_pdfname)); 
  }
  int n_hists = (m_pdfweights) ? m_pdfweights->GetNWeights() : 100;

  for(int i=0; i<n_hists; i++){
    stringstream ss_name;
    ss_name << "Bstar_reco_M_fine_PDF_"  << i+1 ;

    stringstream ss_title;
    ss_title << "M_{tW} [GeV] for PDF No. "  << i+1 << " out of 100" ;

    string s_name = ss_name.str();
    string s_title = ss_title.str();
    const char* char_name = s_name.c_str();
    const char* char_title = s_title.c_str();
    histo_names[i] = s_name;

    
    // double xbins[21] = {500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 4000};
    double xbins[22] = {500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3200, 5500}; // new extended binning


    book<TH1F>(char_name, char_title, 21, xbins);
  }

}

void BstarToTWPDFHists::fill(const Event & event){
  double weight = event.weight;

  if(is_mc)
    {

      if(event.genInfo->systweights().size() == 0 && take_ntupleweights) throw runtime_error("In BstarToTWPDFHists.cxx: Systweights in event.genInfo() is empty but ntupleweights shall be taken. Is this correct? In this case add exception to take_ntupleweights.");    
      if(event.genInfo->systweights().size() != 0 && (is_LO && !m_oname.Contains("DYJets"))) throw runtime_error("In BstarToTWPDFHists.cxx: Systweights in event.genInfo() is NOT empty but this IS a LO sample. Is this correct? In this case Thomas says the genInfo weight should be used. Add this sample to take_ntupleweights");


      std::vector<BstarToTWHypothesis> hyps = event.get(h_hyps);
      const BstarToTWHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );
      if (!hyp)
	{
	  cout << "WARNING: BstarToTWHypothesis: " + m_discriminator_name + " No hypothesis was valid!" << endl;
	  return;
	}

      double mbstar_reco = 0;
      if((hyp->get_topjet() + hyp->get_w()).isTimelike())
	{    
	  mbstar_reco = (hyp->get_topjet() + hyp->get_w()).M();
	}
      else
	{
	  mbstar_reco = sqrt(-(hyp->get_topjet()+hyp->get_w()).mass2());
	}

      //Fill Mbstar (2 cases)
      if(take_ntupleweights)
	{
	  for(int i=0; i<100; i++)
	    {
	      if(use_pdf_weights)
		{
 
		  double pdf_weight = event.genInfo->systweights().at(i+9);
		  double fillweight = weight * pdf_weight/event.genInfo->originalXWGTUP();
		  const char* name = histo_names[i].c_str();

		  if (mbstar_reco <= 4000.) hist(name)->Fill(mbstar_reco, fillweight);
		  else hist(name)->Fill(3999., fillweight);

		}
	    }
	}

      else
	{ 
	  std::vector<double> weights = m_pdfweights->GetWeightList(event);
	  for(unsigned int i=0; i<m_pdfweights->GetNWeights(); i++)
	    {
	      if(use_pdf_weights)
		{
		  double fillweight = weight*weights[i];
		  const char* name = histo_names[i].c_str();

		  if (mbstar_reco <= 4000.) hist(name)->Fill(mbstar_reco, fillweight);
		  else hist(name)->Fill(3999., fillweight);

		}
	    }
	}
    }//is_mc
  
} 


BstarToTWPDFHists::~BstarToTWPDFHists(){}














