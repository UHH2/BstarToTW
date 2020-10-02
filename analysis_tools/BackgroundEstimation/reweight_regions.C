#include "../cosmetics.h"
#include "../plot_styles.C"
#include <vector>
#include <cmath>
#include <TRandom.h>
#include <TRandomGen.h>

void fill_histogram_from_file(TH1F* hist, const TString &fName, const TString &region, const TString &histName, const double &weightSign = 1);
TFitResultPtr fit_ratio(TH1F* hSignal, TH1F* hControl, TF1* ffit);
void fill_reweighted(TH1F* hReweight,TH1F* hReweight_up,TH1F* hReweight_down, TF1* fitfun1, TF1* fitfun2, TFitResultPtr r1, TFitResultPtr r2, const TString &fname, const double &weightSign = 1, double fraction = 1, long int seed = 123456);
void create_output(const TString fout_name, const TString subdir_name, TH1F* hist);

void reweight_regions(TString channel, TString year) {

  const TString controlRegion = "0btag1toptag";
  vector<TString> regions;
  regions.push_back("1btag1toptag20chi2");
  regions.push_back("2btag1toptag");
  vector<TString> inDirs;
  if (channel == "Combined")
    {
      inDirs.push_back("/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Electron/"+year+"/BACKGROUND/");
      inDirs.push_back("/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Muon/"+year+"/BACKGROUND/");
    }
  else
    {
      inDirs.push_back("/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/"+channel+"/"+year+"/BACKGROUND/");
    }
  const TString prefix = "uhh2.AnalysisModuleRunner.";
   
  const int nSamples = 3;
  const TString sampleNames[nSamples] = {"WJets", "DYJets", "Diboson"};
  // const TString sampleNames[nSamples] = {"WJets"};
  const int nTopSamples = 2;
  const TString topSampleNames[nTopSamples] = {"TT", "ST"};

  // read in histograms
  // const int nbins = 17;
  // const double xbins[nbins] = {0, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1700, 1900, 2100, 2400, 3500};
  const int nbins = 21;
  double xbins[nbins] = {500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 4000}; // new extended binning
  

  for (TString region : regions)
    {
      TH1F *hSignal_M = new TH1F("signal_m", channel+"_"+region, nbins-1, xbins);
      TH1F *hControl_M = new TH1F("control_m", channel+"_"+controlRegion, nbins-1, xbins);
      TH1F *hReweight_fit_M = new TH1F("reweight_fit", "_reweight_fit", nbins-1, xbins);
      TH1F *hReweight_fit_M_up = new TH1F("reweight_fit_up", "_reweight_fit", nbins-1, xbins);
      TH1F *hReweight_fit_M_down = new TH1F("reweight_fit_down", "_reweight_fit", nbins-1, xbins);
      TH1F *hReweightInject_fit_M = new TH1F("reweightinject_fit", "_reweight_fit", nbins-1, xbins);
      TH1F *hReweightInject_fit_M_up = new TH1F("reweightinject_fit_up", "_reweight_fit", nbins-1, xbins);
      TH1F *hReweightInject_fit_M_down = new TH1F("reweightinject_fit_down", "_reweight_fit", nbins-1, xbins);
      hReweight_fit_M->Sumw2();
      hReweight_fit_M_up->Sumw2();
      hReweight_fit_M_down->Sumw2();
      hReweightInject_fit_M->Sumw2();
      hReweightInject_fit_M_up->Sumw2();
      hReweightInject_fit_M_down->Sumw2();

      // fill histograms
      std::cout << "get histograms" << std::endl;
      vector<TString> fNames;
      for (TString inDir : inDirs)
	{
	  for (int i = 0; i < nSamples; ++i)
	    {
	      TString fName = inDir + prefix + "MC." + sampleNames[i] + ".root";
	      fNames.push_back(fName);
	      fill_histogram_from_file(hSignal_M, fName, region, "_reco/Bstar_reco_M_rebin");
	      fill_histogram_from_file(hControl_M, fName, controlRegion, "_reco/Bstar_reco_M_rebin");
	    }
	}

      // fit_ratio
      std::cout << "fit ratio" << std::endl;

      TF1 *ffit_gaus = new TF1("ffit_gaus", "gaus(0)+pol0(3)", 500, 5000);
      // if (region == "1btag1toptag20chi2") 
      // 	{
      ffit_gaus->SetParameters(0.01, 1450.0, 1500.0, 0.005);
      ffit_gaus->SetParLimits(0, 0.0, 0.1);
      ffit_gaus->SetParLimits(3, 0.0, 0.1);
      // 	}
      // else 
      // 	{
      // 	  ffit_gaus->SetParameters(0.002, 1450.0, 160.0, 0.005);
      // 	  ffit_gaus->SetParLimits(0, 0.0, 0.1);
      // 	  ffit_gaus->SetParLimits(3, 0.0, 0.1);
      // 	}
      TFitResultPtr r_gaus = fit_ratio(hSignal_M, hControl_M, ffit_gaus);

      TF1 *ffit_landau = new TF1("ffit_landau", "landau", 500, 5000);
      TFitResultPtr r_landau = fit_ratio(hSignal_M, hControl_M, ffit_landau);
      
      // reweight from control region
      double split = 0.66;
      double fraction = (region == "1btag1toptag20chi2") ? split : 1 - split;
      
      std::cout << "start reweighting" << std::endl;
      for (TString inDir : inDirs)
	{
	  std::cout << "reweight data..." << std::endl;
	  fill_reweighted(hReweight_fit_M, hReweight_fit_M_up, hReweight_fit_M_down, ffit_gaus, ffit_landau, r_gaus, r_landau, inDir + prefix+ "DATA.DATA.root", 1, fraction);
	  bool inject_signal = false;
	  double signal_strength = 0.05;
	  TString signal_sample = "MC.BstarToTW3000LH.root";
	  if (inject_signal) 
	    {
	      fill_reweighted(hReweight_fit_M, hReweight_fit_M_up, hReweight_fit_M_down, ffit_gaus, ffit_landau, r_gaus, r_landau, inDir + prefix+ signal_sample, signal_strength, 1.);
	    }
	  for (int i = 0; i < nTopSamples; ++i)
	    {
	      std::cout << "subtract reweighted " <<  topSampleNames[i] <<std::endl;
	      TString fName = inDir + prefix + "MC." + topSampleNames[i] + ".root";
	      fill_reweighted(hReweight_fit_M, hReweight_fit_M_up, hReweight_fit_M_down, ffit_gaus, ffit_landau, r_gaus, r_landau, fName, -1, 1.);
	    }
	}


      create_output(prefix+"MC.Other.root",region+"_reco", hReweight_fit_M);
      create_output(prefix+"MC.Other_UP.root",region+"_reco", hReweight_fit_M_up);
      create_output(prefix+"MC.Other_DOWN.root",region+"_reco", hReweight_fit_M_down);
      
    }
}

void create_output(const TString fout_name, const TString subdir_name, TH1F* hist) {
  TFile *fout = new TFile(fout_name, "update");
  std::cout << "write to file" << std::endl;
  TDirectory* subdir = fout->mkdir(subdir_name);
  subdir->cd();
  TH1F *hOut =(TH1F*) hist->Clone("Bstar_reco_M_rebin");
  hOut->Write();
  fout->Close();
  delete fout;
}

void fill_histogram_from_file(TH1F* hist, const TString &fName, const TString &region, const TString &histName, const double &weightSign) {

  TFile* f = new TFile(fName);
  hist->Add((TH1F*)  f->Get(region + histName), weightSign);
  f->Close();
  delete f;
}

TFitResultPtr fit_ratio(TH1F* hSignal, TH1F* hControl, TF1* ffit) {
  setCMSStyle();

  TGraphAsymmErrors* ratio = new TGraphAsymmErrors();
  ratio->Divide(hSignal, hControl, "pois");
  
  TFitResultPtr r = ratio->Fit(ffit, "S0", "", 500, 5000);
  
  // Draw and safe
  TCanvas* c = new TCanvas("cfit","cfit",600,600);
  TPad* pad = SetupPad();

  gStyle->SetOptFit(1);
  pad->Draw();
  pad->cd();
  ratio->Draw("AP");

  ffit->SetLineColor(kAzure-3);
  ffit->SetLineWidth(2);
  ffit->Draw("same");

  ratio->GetXaxis()->SetTitle("M_{tW} [GeV]");
  ratio->GetYaxis()->SetTitle("#alpha");
  ratio->GetYaxis()->CenterTitle();
  ratio->GetHistogram()->GetYaxis()->SetRangeUser(0.00, 0.4);
  TLegend* leg = new TLegend(0.25,0.75,0.5,0.9);
  leg->SetBorderSize(0);
  leg->AddEntry(ratio,"#alpha_{i}", "lpe");
  leg->AddEntry(ffit, "fit","l");
  
  leg->Draw();

  return r;
}

void fill_reweighted(TH1F* hReweight,TH1F* hReweight_up,TH1F* hReweight_down, TF1* fitfun1, TF1* fitfun2,  TFitResultPtr r1, TFitResultPtr r2, const TString &fname, const double &weightSign, double fraction, long int seed) {
  TFile *f = new TFile(fname);
  TTreeReader tReader("AnalysisTree", f);
  TTreeReaderValue<float> preWeight(tReader, "final_weight");
  TTreeReaderValue<float> recoMass(tReader, "reco_mass");
  TRandom* rng = new TRandomMixMax();
  rng->SetSeed(seed);
  float scale = 1. / fraction;
  while (tReader.Next()) 
    {
      if ( rng->Uniform() < fraction )
	{
	  float mass = (*recoMass);
	  mass = (mass >= 4000.) ? 3999.9 : mass;
	  double x[1] = { mass };
	  double err_1[1], err_2[1];

	  // evaluate fit at masspoint and get 1 sigma deviation from 68.3% confidence interval
	  double diff = fitfun1->Eval(*recoMass) - fitfun2->Eval(*recoMass);
	  float weight = fitfun1->Eval(*recoMass);
	  r1->GetConfidenceIntervals(1, 1, 1, x, err_1, 0.683, false);
	  r2->GetConfidenceIntervals(1, 1, 1, x, err_2, 0.683, false);
	  double err = sqrt( pow(diff/2, 2) + pow(err_1[0], 2) + pow(err_2[0], 2));
	  float weight_up = weight - diff/2 + err;
	  float weight_down = weight - diff/2 - err;
	  weight = (weight < 0) ? 0 : weight * (*preWeight);
	  weight_up = (weight_up < 0) ? 0 : weight_up * (*preWeight);
	  weight_down = (weight_down < 0) ? 0 : weight_down * (*preWeight);
	  hReweight->Fill(mass, weight * weightSign * scale );
	  hReweight_up->Fill(mass, weight_up * weightSign * scale );
	  hReweight_down->Fill(mass, weight_down * weightSign * scale );
	}
    }
  f->Close();
  delete f;
}


