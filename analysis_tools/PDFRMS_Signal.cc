#include "../include/Tools.h"
#include <TStyle.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TRint.h>
#include <TClass.h>
#include <TKey.h>

using namespace std;

void AnalysisTool::PDFRMS_Signal(bool isHT){

  TString base_path, out_path;
  base_path = AnalysisTool::base_path_SR;
  out_path = base_path;
  base_path += "NOMINAL/";


  TString hist_name;
  if(isHT) hist_name = "H_T_comb_NoMLQ_from350_all_filled_rebin";
  else     hist_name = "M_LQ_MuEle_comb_all_filled";
  TString hist_path = "FinalSelection/" + hist_name ;

  cout << "Going to open files. Example: '" << base_path + "uhh2.AnalysisModuleRunner.MC.[LQsample].root', opening histogram '" << hist_path << "'" <<  endl;

  vector<TString> masspoints;
  masspoints.push_back("200");
  masspoints.push_back("300");
  masspoints.push_back("400");
  masspoints.push_back("500");
  masspoints.push_back("600");
  masspoints.push_back("700");
  masspoints.push_back("800");
  masspoints.push_back("900");
  masspoints.push_back("1000");
  masspoints.push_back("1200");
  masspoints.push_back("1400");
  masspoints.push_back("1700");
  masspoints.push_back("2000");

  vector<TFile*> files;
  for(unsigned int i=0; i<masspoints.size(); i++){
    files.emplace_back(new TFile(base_path + "uhh2.AnalysisModuleRunner.MC.LQtoTMuM" + masspoints.at(i) + "_norm.root","READ"));
  }

  vector<vector<TH1D*>> hists;
  
  cout << "1" << endl;
  // histograms
  vector<TH1D*> h_nominal;
  for(unsigned int i=0; i<masspoints.size(); i++){
    h_nominal.emplace_back((TH1D*)files.at(i)->Get(hist_path));
  }

  //4%0 -->In binning for alpha-extrapolation
  //4%1 -->For extrapolating TTbar(+DY) from CR to SR (only used in CR)
  //4%2 -->MLQ only for SR
  //4%3 -->HT no ele only for SR

  cout << "2" << endl;
  for(unsigned int i=0; i<masspoints.size(); i++){
    int n=0;
    files.at(i)->cd("PDF_variations");
    TDirectory* dir = gDirectory;
    TIter iter(dir->GetListOfKeys());
    TKey *key;
    TH1D* h;
    vector<TH1D*> hists_tmp;
    while ((key = (TKey*)iter())) {
      if( (isHT && n%4!=3) || (!isHT && n%4!=2)) {n++;continue;}
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      h = (TH1D*)key->ReadObj();
      hists_tmp.push_back(h);
      n++;
      //cout << "naechstes histogramm " <<  n <<  endl;
    }
    hists.push_back(hists_tmp);
  }

  cout << "3" << endl;
 
  /*
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  int n_low_edges;
  if(isHT)  n_low_edges = 16;
  else n_low_edges = 7;

  //if(hists_ht_DY.at(0)->GetNbinsX()+1 != 21) throw runtime_error("number of bin-low edges is not 21, please adjust.");
  double low_edges[n_low_edges];
  for(int i=1; i<hists_ht_DY.at(0)->GetNbinsX()+2; i++){
    low_edges[i-1] = hists_ht_DY.at(0)->GetBinLowEdge(i);
  }
  //if(hists_ht_DY.at(0)->GetNbinsX() != 20) throw runtime_error("number of bins is not 20, please adjust.");
  */

  vector<double> low_edges;
  for(int i=1; i<h_nominal.at(0)->GetNbinsX()+2; i++){
    low_edges.push_back(h_nominal.at(0)->GetBinLowEdge(i));
  }
  int n_low_edges = low_edges.size();
  cout << "3.1" << endl;



  vector<TH1D*> h_up, h_down;
  for(int i=0; i<masspoints.size(); i++){
    h_up.emplace_back(  new TH1D(hist_name,"LQ"+masspoints.at(i)+", PDF up",n_low_edges-1,&low_edges[0]));	       
    h_down.emplace_back(new TH1D(hist_name,"LQ"+masspoints.at(i)+", PDF down",n_low_edges-1,&low_edges[0]));          
  } 
       
  cout << "3.2" << endl;                                                 

  //calculate RMS and fill up and down histograms
  vector<vector<double>> v_rms;
  for(unsigned int k=0; k<masspoints.size(); k++){
    vector<double> v_rms_tmp;
    for(int j=1; j<h_nominal.at(0)->GetNbinsX()+1; j++){
      double nominal = h_nominal.at(k)->GetBinContent(j);

      double rms=0;
      for(unsigned int i=0; i<hists.at(0).size(); i++){
	rms += pow((hists.at(k).at(i)->GetBinContent(j)-nominal),2);
      }
      rms /= (hists.at(k).size()-1);

      v_rms_tmp.push_back(sqrt(rms));

      h_up.at(k)->SetBinContent(j,nominal + v_rms_tmp.at(j-1));
      h_up.at(k)->SetBinError(j,h_nominal.at(k)->GetBinError(j));
      h_down.at(k)->SetBinContent(j,max(0.,nominal - v_rms_tmp.at(j-1)));
      h_down.at(k)->SetBinError(j,h_nominal.at(k)->GetBinError(j));

      cout << "MassPoint: " << masspoints.at(k) << ", Content in bin " << j << ": " <<  nominal << " +- " << v_rms_tmp.at(j-1) << endl;
    }
    v_rms.push_back(v_rms_tmp);
    cout << "3.3" << endl; 
  }
 
  cout << "4" << endl;
  
  // for limit histograms
  for(unsigned int i=0; i<masspoints.size(); i++){                                                       
    TFile* f_out_up = new TFile(out_path + "PDF_up/LQtoTMuM" + masspoints.at(i) + "_norm.root","UPDATE");
    f_out_up->mkdir("FinalSelection");
    f_out_up->cd("FinalSelection");
    h_up.at(i)->Write();
    f_out_up->Close();
    TFile* f_out_dn = new TFile(out_path + "PDF_down/LQtoTMuM" + masspoints.at(i) + "_norm.root","UPDATE");
    f_out_dn->mkdir("FinalSelection");
    f_out_dn->cd("FinalSelection");
    h_down.at(i)->Write();
    f_out_dn->Close();

    delete f_out_up;
    delete f_out_dn;
  }
  
  cout << "5" << endl;
  cout << "done." << endl;












  for(unsigned int i=0; i<masspoints.size(); i++) delete h_down[i];
  for(unsigned int i=0; i<masspoints.size(); i++) delete h_up[i];
  for(unsigned int i=0; i<masspoints.size(); i++) delete h_nominal[i];
  for(unsigned int i=0; i<masspoints.size(); i++) delete files[i];

}









