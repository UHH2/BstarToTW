#pragma once
#include <TString.h>
#include <TVirtualFitter.h>
#include <iostream>
#include <TStyle.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TText.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <TFitResult.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
#include <sstream>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <string>


using namespace std;


void FindBiggestDeviations(){

  //In each bin find the lowest and the highest entry of all 6 variations + nominal
  //lowest: down-variation for theta, highest: up-variation for theta, nominal: nominal for theta

  // TString path = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Muon/2016/";
  TString path = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Electron/2018/";


  // TString signal_postfix = "_2016v3";
  // TString signal_postfix = "_2017v2";
  TString signal_postfix = "_2018";

  // TString channel = "1btag1toptag20chi2";
  TString channel = "2btag1toptag";
 
  // TString hist_name = channel + "_reco/Bstar_reco_M_rebin";

  const int nproc = 89;
  TString proc[nproc] = {"TT",  "ST",
			 "BstarToTW0700RH"+signal_postfix, "BstarToTW0800RH"+signal_postfix, "BstarToTW0900RH"+signal_postfix, "BstarToTW1000RH"+signal_postfix, "BstarToTW1100RH"+signal_postfix,
			 "BstarToTW1200RH"+signal_postfix, "BstarToTW1400RH"+signal_postfix, "BstarToTW1600RH"+signal_postfix, "BstarToTW1800RH"+signal_postfix,
			 "BstarToTW2000RH"+signal_postfix, "BstarToTW2200RH"+signal_postfix, "BstarToTW2400RH"+signal_postfix, "BstarToTW2600RH"+signal_postfix, 
			 "BstarToTW2800RH"+signal_postfix, "BstarToTW3000RH"+signal_postfix, "BstarToTW3200RH"+signal_postfix, "BstarToTW3400RH"+signal_postfix,
			 "BstarToTW3600RH"+signal_postfix, "BstarToTW3800RH"+signal_postfix, "BstarToTW4000RH"+signal_postfix, "BstarToTW4200RH"+signal_postfix,
			 "BstarToTW0700LH"+signal_postfix, "BstarToTW0800LH"+signal_postfix, "BstarToTW0900LH"+signal_postfix, "BstarToTW1000LH"+signal_postfix, "BstarToTW1100LH"+signal_postfix,
			 "BstarToTW1200LH"+signal_postfix, "BstarToTW1400LH"+signal_postfix, "BstarToTW1600LH"+signal_postfix, "BstarToTW1800LH"+signal_postfix,
			 "BstarToTW2000LH"+signal_postfix, "BstarToTW2200LH"+signal_postfix, "BstarToTW2400LH"+signal_postfix, "BstarToTW2600LH"+signal_postfix, 
			 "BstarToTW2800LH"+signal_postfix, "BstarToTW3000LH"+signal_postfix, "BstarToTW3200LH"+signal_postfix, "BstarToTW3400LH"+signal_postfix,
			 "BstarToTW3600LH"+signal_postfix, "BstarToTW3800LH"+signal_postfix, "BstarToTW4000LH"+signal_postfix, "BstarToTW4200LH"+signal_postfix,
			 "BstarToTW0700VL"+signal_postfix, "BstarToTW0800VL"+signal_postfix, "BstarToTW0900VL"+signal_postfix, "BstarToTW1000VL"+signal_postfix, "BstarToTW1100VL"+signal_postfix,
			 "BstarToTW1200VL"+signal_postfix, "BstarToTW1400VL"+signal_postfix, "BstarToTW1600VL"+signal_postfix, "BstarToTW1800VL"+signal_postfix,
 			 "BstarToTW2000VL"+signal_postfix, "BstarToTW2200VL"+signal_postfix, "BstarToTW2400VL"+signal_postfix, "BstarToTW2600VL"+signal_postfix, 
			 "BstarToTW2800VL"+signal_postfix, "BstarToTW3000VL"+signal_postfix, "BstarToTW3200VL"+signal_postfix, "BstarToTW3400VL"+signal_postfix,
			 "BstarToTW3600VL"+signal_postfix, "BstarToTW3800VL"+signal_postfix, "BstarToTW4000VL"+signal_postfix, "BstarToTW4200VL"+signal_postfix,
			 "BprimeTToTW0700RH_2016v3", "BprimeTToTW0800RH_2016v3", "BprimeTToTW0900RH_2016v3", "BprimeTToTW1000RH_2016v3", "BprimeTToTW1100RH_2016v3", "BprimeTToTW1200RH_2016v3",
			 "BprimeTToTW1300RH_2016v3", "BprimeTToTW1400RH_2016v3", "BprimeTToTW1500RH_2016v3", "BprimeTToTW1600RH_2016v3", "BprimeTToTW1700RH_2016v3", "BprimeTToTW1800RH_2016v3",
			 "BprimeTToTW0700LH_2016v3", "BprimeTToTW0800LH_2016v3", "BprimeTToTW0900LH_2016v3", "BprimeTToTW1000LH_2016v3", "BprimeTToTW1100LH_2016v3", "BprimeTToTW1200LH_2016v3",
			 "BprimeTToTW1300LH_2016v3", "BprimeTToTW1400LH_2016v3", "BprimeTToTW1500LH_2016v3", "BprimeTToTW1600LH_2016v3", "BprimeTToTW1700LH_2016v3", "BprimeTToTW1800LH_2016v3"};


  for (int n = 0; n < nproc; ++n)
    {

      TString process = proc[n];

      //Files & histograms for other processes
      unique_ptr<TFile> f_uu, f_un, f_nu, f_nd, f_dn, f_dd, f_nom;
      unique_ptr<TH1F> h_uu, h_un, h_nu, h_nd, h_dn, h_dd, h_nom;
      cout << (path + "SCALE_upup/uhh2.AnalysisModuleRunner.MC." + process  + ".root") << endl;
 
      f_uu.reset (new TFile(path + "SCALE_upup/uhh2.AnalysisModuleRunner.MC." + process  + ".root","READ"));
      h_uu.reset((TH1F*)f_uu->Get(channel+"_reco/Bstar_reco_M_rebin")); 
      f_un.reset (new TFile(path + "SCALE_upnone/uhh2.AnalysisModuleRunner.MC." + process  + ".root","READ"));
      h_un.reset((TH1F*)f_un->Get(channel+"_reco/Bstar_reco_M_rebin")); 
      f_nu.reset (new TFile(path + "SCALE_noneup/uhh2.AnalysisModuleRunner.MC." + process  + ".root","READ"));
      h_nu.reset((TH1F*)f_nu->Get(channel+"_reco/Bstar_reco_M_rebin")); 
      f_nd.reset (new TFile(path + "SCALE_nonedown/uhh2.AnalysisModuleRunner.MC." + process  + ".root","READ"));
      h_nd.reset((TH1F*)f_nd->Get(channel+"_reco/Bstar_reco_M_rebin")); 
      f_dn.reset (new TFile(path + "SCALE_downnone/uhh2.AnalysisModuleRunner.MC." + process  + ".root","READ"));
      h_dn.reset((TH1F*)f_dn->Get(channel+"_reco/Bstar_reco_M_rebin")); 
      f_dd.reset (new TFile(path + "SCALE_downdown/uhh2.AnalysisModuleRunner.MC." + process  + ".root","READ"));
      h_dd.reset((TH1F*)f_dd->Get(channel+"_reco/Bstar_reco_M_rebin")); 
      f_nom.reset(new TFile(path + "NOMINAL/uhh2.AnalysisModuleRunner.MC." + process  + ".root","READ"));
      h_nom.reset((TH1F*)f_nom->Get(channel+"_reco/Bstar_reco_M_rebin")); 


      const int nbins = h_nom->GetNbinsX();
      cout << "Number of bins: " << nbins << endl;
      vector<double> min_bins, max_bins, max_err, min_err;

      for(int i=1; i<h_nom->GetNbinsX()+1; i++){
	double entries[7] = {h_uu->GetBinContent(i),h_un->GetBinContent(i),h_nu->GetBinContent(i),h_nd->GetBinContent(i),h_dn->GetBinContent(i),h_dd->GetBinContent(i),h_nom->GetBinContent(i)};
	double errors[7] = {h_uu->GetBinError(i),h_un->GetBinError(i),h_nu->GetBinError(i),h_nd->GetBinError(i),h_dn->GetBinError(i),h_dd->GetBinError(i),h_nom->GetBinError(i)};
	double min = DBL_MAX;
	double max_error = 0, min_error = 0;
	double max = 0;
	for(int j=0; j<7; j++){
	  if(entries[j] > max) {max = entries[j]; max_error = errors[j];}
	  if(entries[j] < min) {min = entries[j]; min_error = errors[j];}
	}

	min_bins.push_back(min);
	max_bins.push_back(max);
	min_err.push_back(min_error);
	max_err.push_back(max_error);
      }

      cout << "minimum entries: " << endl;
      for(int i=0; i<h_nom->GetNbinsX(); i++){
	cout << min_bins[i] << " ";
      }
      cout << endl;

      cout << "maximum entries: " << endl;
      for(int i=0; i<h_nom->GetNbinsX(); i++){
	cout << max_bins[i] << " ";
      }
      cout << endl;

      double xbins[21] = {500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 4000}; // new extended binning

      unique_ptr<TH1F> hist_out_up, hist_out_dn;
      hist_out_up.reset(new TH1F("Bstar_reco_M_rebin", "M_{tW} [GeV]", 20, xbins));
      for(int i=1; i<h_nom->GetNbinsX()+1; i++){
	hist_out_up->SetBinContent(i,max_bins[i-1]);
	hist_out_up->SetBinError(i, max_err[i-1]);
      }

      cout << "test: " << hist_out_up->GetBinContent(3) << endl;
      cout << "out_path: " << path << "SCALE_up/uhh2.AnalysisModuleRunner.MC." << process << ".root" << endl;

      unique_ptr<TFile> f_out_up, f_out_dn;
      f_out_up.reset(new TFile(path + "SCALE_up/uhh2.AnalysisModuleRunner.MC." + process  + ".root","UPDATE"));
      f_out_up->mkdir(channel+"_reco");
      f_out_up->cd(channel+"_reco");
      hist_out_up->Write();
      f_out_up->Close();
      
      hist_out_up.reset();
      
      hist_out_dn.reset(new TH1F("Bstar_reco_M_rebin", "M_{tW} [GeV]", 20, xbins));
      for(int i=1; i<h_nom->GetNbinsX()+1; i++){
	hist_out_dn->SetBinContent(i,min_bins[i-1]);
	hist_out_dn->SetBinError(i, min_err[i-1]);
      }

      cout << "test: " << hist_out_dn->GetBinContent(3) << endl;
      cout << "out_path: " << path << "SCALE_down/uhh2.AnalysisModuleRunner.MC." << process << ".root" << endl;

      f_out_dn.reset(new TFile(path + "SCALE_down/uhh2.AnalysisModuleRunner.MC." + process  + ".root","UPDATE"));
      f_out_dn->mkdir(channel+"_reco");
      f_out_dn->cd(channel+"_reco");
      hist_out_dn->Write();
      f_out_dn->Close();

      hist_out_dn.reset();

      cout << "wrote files." << endl;
    }
}
