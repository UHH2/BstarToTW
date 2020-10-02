#pragma once
#include <TStyle.h>
#include <TH1.h>
#include <TH1F.h>
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

void PDFRMS_Background(){

  TString base_path = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Muon/2016/";
  // TString base_path = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Electron/2016/";
  TString out_path = base_path;
  base_path += "NOMINAL/";

  TString signal_postfix = "_2016v3";
  // TString signal_postfix = "_2017v2";
  // TString signal_postfix = "_2018";
  // TString signal_postfix = "";

  TString hist_path = "1btag1toptag20chi2";
  // TString hist_path = "2btag1toptag";

  TString hist_name = hist_path + "_reco/Bstar_reco_M_rebin";

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
      cout << "Processing " << proc[n] << "..." << endl;
      unique_ptr<TFile> file;
      file.reset(new TFile(base_path + "uhh2.AnalysisModuleRunner.MC." + proc[n] + ".root","READ"));
  
  
      // histograms
      unique_ptr<TH1F> h_nominal;
      h_nominal.reset((TH1F*)file->Get(hist_name));

      vector<TH1F*> hists;
      file->cd(hist_path + "_PDF_variations");
      TDirectory* dir = gDirectory;
      TIter iter(dir->GetListOfKeys());
      TKey *key;
      TH1F* h;
      while ((key = (TKey*)iter())) 
	{
	  TClass *cl = gROOT->GetClass(key->GetClassName());
	  if (!cl->InheritsFrom("TH1")) continue;
	  h = (TH1F*)key->ReadObj();
	  hists.push_back(std::move(h));
	}

      //calculate RMS 
      vector<double> rms_up, rms_dn, err;
      for(int j=1; j < h_nominal->GetNbinsX()+1; j++)
	{
	  double nominal = h_nominal->GetBinContent(j);

	  double rms=0;
	  for(unsigned int i=0; i<hists.size(); i++)
	    {
	      rms += pow((hists.at(i)->GetBinContent(j)-nominal),2);
	    }
	  rms /= (hists.size());
	  rms_up.push_back(nominal + sqrt(rms));
	  rms_dn.push_back(nominal - sqrt(rms));
	  err.push_back(h_nominal->GetBinError(j));
	}
      //and fill up and down histograms
      double xbins[21] = {500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 4000}; // new extended binning

      unique_ptr<TH1F> h_up;
      h_up.reset(new TH1F("Bstar_reco_M_rebin", proc[n]+", PDF up", 20, xbins));
      for (int i = 0; i < rms_up.size(); ++i)
	{
	  h_up->SetBinContent(i+1, rms_up.at(i));
	  h_up->SetBinError(i+1, err.at(i));
	}                                             
      unique_ptr<TFile> f_out_up;
      f_out_up.reset(new TFile(out_path + "PDF_up/uhh2.AnalysisModuleRunner.MC." + proc[n] + ".root","UPDATE"));
      f_out_up->mkdir(hist_path+"_reco");
      f_out_up->cd(hist_path+"_reco");
      h_up->Write();
      f_out_up->Close();
      h_up.reset();

      unique_ptr<TH1F> h_dn;
      h_dn.reset(new TH1F("Bstar_reco_M_rebin", proc[n]+", PDF down", 20, xbins));  
      for (int i = 0; i < rms_dn.size(); ++i)
	{
	  h_dn->SetBinContent(i+1, rms_dn.at(i));
	  h_dn->SetBinError(i+1, err.at(i));
	}  
      unique_ptr<TFile> f_out_dn;
      f_out_dn.reset(new TFile(out_path + "PDF_down/uhh2.AnalysisModuleRunner.MC." + proc[n] + ".root","UPDATE"));
      f_out_dn->mkdir(hist_path+"_reco");
      f_out_dn->cd(hist_path+"_reco");
      h_dn->Write();
      f_out_dn->Close();
      h_dn.reset();

      cout << "done." << endl;

    }
}
