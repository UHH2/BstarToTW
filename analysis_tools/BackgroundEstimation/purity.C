#include "../cosmetics.h"
#include "../plot_styles.C"

void fill_histogram_from_file(TH1F* hist, const TString &fName, const TString &region, const TString &histName, const int &weightSign = 1);

void purity(TString channel, const TString year = "2018", const TString region = "1btag1toptag20chi2") {

  TString inDir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/"+channel+"/"+year+"/NOMINAL/";
  const TString prefix = "uhh2.AnalysisModuleRunner.";

  const int nSamples = 4;
  const TString sampleNames[nSamples] = {"WJets", "DYJets", "Diboson", "QCD"};
  const int nTopSamples = 6 - nSamples;
  const TString topSampleNames[nTopSamples] = {"TT", "ST"};

  const int nbins_rebin = 23;
  double xbins_rebin[nbins_rebin] = {0, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3000, 3500};

  TH1F *h_all = new TH1F("all", "all", nbins_rebin-1, xbins_rebin);
  TH1F *h_purity = new TH1F("purity", "purity", nbins_rebin-1, xbins_rebin);

  for (int i = 0; i < nSamples; ++i)
    {
      TString fName = inDir + prefix + "MC." + sampleNames[i] + ".root";
      fill_histogram_from_file(h_all, fName, region, "_reco/Bstar_reco_M_rebin");      
      fill_histogram_from_file(h_purity, fName, region, "_reco/Bstar_reco_M_rebin");      
   }

  for (int i = 0; i < nTopSamples; ++i)
    {
      TString fName = inDir + prefix + "MC." + topSampleNames[i] + ".root";
      fill_histogram_from_file(h_all, fName, region, "_reco/Bstar_reco_M_rebin");
    }

  setCMSStyle();
  TGraphAsymmErrors* ratio = new TGraphAsymmErrors();
  ratio->Divide(h_purity, h_all);
  ratio->GetHistogram()->GetXaxis()->SetTitle("M_{tW} [GeV]");
  ratio->GetHistogram()->GetYaxis()->SetTitle("purity [arb.u.]");
  TCanvas* c = new TCanvas("cfit","cfit",600,600);
  TPad* pad = SetupPad();
  
  pad->Draw();
  pad->cd();
  ratio->Draw("AP");
  c->Print("purity_"+region+"_"+year+".eps");
}

void fill_histogram_from_file(TH1F* hist, const TString &fName, const TString &region, const TString &histName, const int &weightSign) {

  TFile* f = new TFile(fName);
  hist->Add((TH1F*)  f->Get(region + histName), weightSign);
  f->Close();
  delete f;
}
