#include "../cosmetics.h"
#include "../plot_styles.C"
#include <vector>
#include <cmath>

void fill_histogram_from_file(TH1F* hist, const TString &fName, const TString &region, const TString &histName, const int &weightSign = 1);
TGraphAsymmErrors* fit_ratio(TH1F* hSignal, TH1F* hControl);


void plot_transfer(TString channel, const TString year = "2018", const TString signalRegion = "1btag1toptag20chi2") {

  const TString controlRegion = "0btag1toptag_tw";

  vector<TString> inDirs;
  if (channel == "Combined")
    {
      inDirs.push_back("/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Electron/"+year+"/WORKING/");
      inDirs.push_back("/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Muon/"+year+"/WORKING/");
    }
  else
    {
      inDirs.push_back("/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/"+channel+"/"+year+"/WORKING/");
    }
  const TString prefix = "uhh2.AnalysisModuleRunner.";
   
  const int nSamples = 4;
  const TString sampleNames[nSamples] = {"WJets", "DYJets", "Diboson", "QCD"};
  const int nTopSamples = 6 - nSamples;
  const TString topSampleNames[nTopSamples] = {"TT", "ST"};
  // read in histograms
  const int nbins = 21;
  double xbins[nbins] = {500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 4000}; // new extended binning
  double xbins_tev[nbins] = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 4.0}; // new extended binning
  
  const int nbins_rebin = 18;  
  double xbins_rebin[nbins] = {500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 2000, 2400, 2800, 4000}; // new extended binning


  TH1F *hSignal_ST = new TH1F("signal_st", channel+"_"+signalRegion, nbins-1, xbins);
  TH1F *hControl_ST = new TH1F("control_st", channel+"_"+controlRegion, nbins-1, xbins);

  vector<TString> fNames;
  for (TString inDir : inDirs)
    {
      for (int i = 0; i < nSamples; ++i)
	{
	  TString fName = inDir + prefix + "MC." + sampleNames[i] + ".root";
	  fNames.push_back(fName);
	  if (sampleNames[i] == "QCD") continue;
	  fill_histogram_from_file(hSignal_ST, fName, signalRegion, "_reco/Bstar_reco_M_rebin");
	  fill_histogram_from_file(hControl_ST, fName, controlRegion, "_reco/Bstar_reco_M_rebin");
	}
    }

  hSignal_ST->SetBinErrorOption( TH1::kPoisson);
  hControl_ST->SetBinErrorOption( TH1::kPoisson);

  // hSignal_ST->SetBins(nbins-1, xbins_tev);
  // hControl_ST->SetBins(nbins-1, xbins_tev);

  // hSignal_ST = (TH1F*)hSignal_ST->Rebin(nbins_rebin-1, "signal_rebin", xbins_rebin);
  // hControl_ST = (TH1F*)hControl_ST->Rebin(nbins_rebin-1, "control_rebin", xbins_rebin);
  // fit_ratio
  TGraphAsymmErrors* ratio = fit_ratio(hSignal_ST, hControl_ST);
  

}

void fill_histogram_from_file(TH1F* hist, const TString &fName, const TString &region, const TString &histName, const int &weightSign) {

  TFile* f = new TFile(fName);
  hist->Add((TH1F*)  f->Get(region + histName), weightSign);
  f->Close();
  delete f;
}

TGraphAsymmErrors* fit_ratio(TH1F* hSignal, TH1F* hControl) {
  setCMSStyle();

  TGraphAsymmErrors* ratio = new TGraphAsymmErrors();
  TGraphAsymmErrors* fit_err = new TGraphAsymmErrors();
  ratio->Divide(hSignal, hControl, "pois");
  TGraphAsymmErrors* ratio2 =(TGraphAsymmErrors*) ratio->Clone("ratio2");
  
  TF1 *ffit_gaus = new TF1("ffit_gaus", "gaus(0)+pol0(3)", 500, 5000);
  ffit_gaus->SetParameters(0.01, 1450.0, 1500.0, 0.005);
  ffit_gaus->SetParLimits(0, 0.0, 0.1);
  ffit_gaus->SetParLimits(3, 0.0, 0.1);
  TFitResultPtr r_gaus =  ratio2->Fit(ffit_gaus, "S0+", "", 500, 5000);

  TF1 *ffit_landau = new TF1("ffit_landau", "landau", 500, 5000);
  TFitResultPtr r_landau = ratio->Fit(ffit_landau, "S0+", "", 500, 5000);

  // TF1 *ffit = new TF1("ffit", "pol2", 500, 5000);
  
 
  int n_points = 10000;
  double xmax,xmin;
  xmin = hSignal->GetXaxis()->GetXmin();
  xmax = hSignal->GetXaxis()->GetXmax();
  for (int i = 0; i < n_points; ++i)
    {
      double x, y, diff, err_up, err_down;
      double err_gaus[1];
      double err_landau[1];
      x = xmin + xmax - ( xmax * (float) i / n_points);
      double xi[1] = {x};
      diff = ffit_gaus->Eval(x) - ffit_landau->Eval(x);
      y = ffit_gaus->Eval(x) - (diff / 2.);
      r_gaus->GetConfidenceIntervals(1, 1, 1, xi, err_gaus, 0.683, false);
      r_landau->GetConfidenceIntervals(1, 1, 1, xi, err_landau, 0.683, false);
      err_up = sqrt( pow(diff/2, 2) + TMath::Power(err_gaus[0], 2) + TMath::Power(err_landau[0], 2));
      err_down = sqrt( pow(diff/2, 2) + TMath::Power(err_gaus[0], 2) + TMath::Power(err_landau[0], 2));
      fit_err->SetPoint(i,x, y);
      fit_err->SetPointError(i,0,0,err_up,err_down);
    }
  

 // Draw and safe
  TCanvas* c = new TCanvas("cfit","cfit",600,600);
  TPad* pad = SetupPad();

  gStyle->SetOptFit(1);
  pad->Draw();
  pad->cd();

  fit_err->SetFillColor(kGreen+1);
  fit_err->Draw("A3");
  fit_err->GetXaxis()->SetTitle("M_{tW} [GeV]");
  fit_err->GetYaxis()->SetTitle("#alpha");
  fit_err->GetYaxis()->CenterTitle();
  fit_err->GetHistogram()->GetYaxis()->SetRangeUser(0, 0.24);

  ffit_gaus->SetLineColor(kAzure-3);
  ffit_gaus->SetLineWidth(2);
  ffit_gaus->Draw("same");
  ffit_landau->SetLineColor(kViolet-2);
  ffit_landau->SetLineWidth(2);
  ffit_landau->Draw("same");
  ratio->Draw("P same");
  ratio2->Draw("P sames");
  c->Update();
  TPaveStats *p2 = (TPaveStats*)ratio2->FindObject("stats");
  p2->SetY1NDC(0.7);
  p2->SetY2NDC(0.85);
  c->Modified();
  TLegend* leg = new TLegend(0.25,0.66,0.5,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(ratio,"#alpha", "lpe");
  leg->AddEntry(ffit_gaus, "gaus fit","l");
  leg->AddEntry(ffit_landau, "pol0 fit","l");
  leg->AddEntry(fit_err, "uncertainty","f");
  
  leg->Draw();

  TString outputName; outputName.Form("ratiofit_%s.eps", hSignal->GetTitle());
  c->Print(outputName);

  return ratio;
}

