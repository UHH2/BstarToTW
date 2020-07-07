#include "cosmetics.h"
#include "plot_styles.C"

#include "TH2.h"
#include "TMath.h"


Double_t median1(TH1 *h1) {
  //compute the median for 1-d histogram h1
  Int_t nbins = h1->GetXaxis()->GetNbins();
  Double_t *x = new Double_t[nbins];
  Double_t *y = new Double_t[nbins];
  for (Int_t i=0;i<nbins;i++) {
    x[i] = h1->GetXaxis()->GetBinCenter(i+1);
    y[i] = h1->GetBinContent(i+1);
  }
  Double_t median = TMath::Median(nbins,x,y);
  delete [] x;
  delete [] y;
  return median;
}

void median2(TH2 *h2) {
  //compute and print the median for each slice along X of h2

  Int_t nbins = h2->GetXaxis()->GetNbins();
  for (Int_t i=1;i<=nbins;i++) {
    TH1 *h1 = h2->ProjectionY("",i,i);
    Double_t median = median1(h1);
    Double_t mean = h1->GetMean();
    printf("Median of Slice %d, Median=%g, Mean = %g\n",i,median,mean);
    delete h1;
  }
}

TH2F* get_hist_from_file(const TString file_name, const TString hist_name) {
  TFile *f = new TFile(file_name);
  TH2F *hist = (TH2F*)( f->Get(hist_name))->Clone();
  hist->SetDirectory(0);
  f->Close();
  delete f;
  return hist;
}

void hotvr_resolution() {
  
  TString in_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/HOTVR_JEC/2018/";
  const TString prefix = "uhh2.AnalysisModuleRunner.MC.";
  
  const TString hist_name = "delta_pt_gen_reco";
  // const TString hist_name = "delta_pt_gen_reco_subjet";

  TH2F *h_2016_pre = get_hist_from_file(in_dir+prefix+"TT_2016v3.root", "HOTVR_Pre_Performance/"+hist_name);
  h_2016_pre->SetName("h_2016_pre");
  // median2(h_2016);
  h_2016_pre->FitSlicesY();
  TH1F *h_2016_pre_1 = (TH1F*)gDirectory->Get("h_2016_pre_1");
  // TH1F *h_2016_err = (TH1F*)gDirectory->Get("h_2016_2");
  // for (int i = 1; i < h_2016_1->GetNbinsX(); ++i) h_2016_1->SetBinError(i, h_2016_err->GetBinContent(i));

  TH2F *h_2017_pre = get_hist_from_file(in_dir+prefix+"TT_2017v2.root", "HOTVR_Pre_Performance/"+hist_name);
  h_2017_pre->SetName("h_2017_pre");
  // median2(h_2017);
  h_2017_pre->FitSlicesY();
  TH1F *h_2017_pre_1 = (TH1F*)gDirectory->Get("h_2017_pre_1");
  // TH1F *h_2017_err = (TH1F*)gDirectory->Get("h_2017_2");
  // for (int i = 1; i < h_2017_1->GetNbinsX(); ++i) h_2017_1->SetBinError(i, h_2017_err->GetBinContent(i));

  TH2F *h_2018_pre = get_hist_from_file(in_dir+prefix+"TT_2018.root", "HOTVR_Pre_Performance/"+hist_name);
  h_2018_pre->SetName("h_2018_pre");
  // median2(h_2018);
  h_2018_pre->FitSlicesY();
  TH1F *h_2018_pre_1 = (TH1F*)gDirectory->Get("h_2018_pre_1");
  // TH1F *h_2018_err = (TH1F*)gDirectory->Get("h_2018_2");
  // for (int i = 1; i < h_2018_1->GetNbinsX(); ++i) h_2018_1->SetBinError(i, h_2018_err->GetBinContent(i));


  TH2F *h_2016 = get_hist_from_file(in_dir+prefix+"TT_2016v3.root", "HOTVR_Performance/"+hist_name);
  h_2016->SetName("h_2016");
  // median2(h_2016);
  h_2016->FitSlicesY();
  TH1F *h_2016_1 = (TH1F*)gDirectory->Get("h_2016_1");
  // TH1F *h_2016_err = (TH1F*)gDirectory->Get("h_2016_2");
  // for (int i = 1; i < h_2016_1->GetNbinsX(); ++i) h_2016_1->SetBinError(i, h_2016_err->GetBinContent(i));

  TH2F *h_2017 = get_hist_from_file(in_dir+prefix+"TT_2017v2.root", "HOTVR_Performance/"+hist_name);
  h_2017->SetName("h_2017");
  // median2(h_2017);
  h_2017->FitSlicesY();
  TH1F *h_2017_1 = (TH1F*)gDirectory->Get("h_2017_1");
  // TH1F *h_2017_err = (TH1F*)gDirectory->Get("h_2017_2");
  // for (int i = 1; i < h_2017_1->GetNbinsX(); ++i) h_2017_1->SetBinError(i, h_2017_err->GetBinContent(i));

  TH2F *h_2018 = get_hist_from_file(in_dir+prefix+"TT_2018.root", "HOTVR_Performance/"+hist_name);
  h_2018->SetName("h_2018");
  // median2(h_2018);
  h_2018->FitSlicesY();
  TH1F *h_2018_1 = (TH1F*)gDirectory->Get("h_2018_1");
  // TH1F *h_2018_err = (TH1F*)gDirectory->Get("h_2018_2");
  // for (int i = 1; i < h_2018_1->GetNbinsX(); ++i) h_2018_1->SetBinError(i, h_2018_err->GetBinContent(i));


  setCMSStyle();
  h_2016_1->GetXaxis()->SetTitle("p_{T}^{gen} [GeV]");
  h_2016_1->GetXaxis()->SetTitleOffset(1.2);
  // h_2016_1->GetYaxis()->SetTitle("N_{gensubjets}");
  h_2016_1->GetYaxis()->SetTitle("#LT p_{T}^{rec}/p_{T}^{gen} #GT");
  h_2016_1->GetYaxis()->CenterTitle();
  h_2016_1->GetYaxis()->SetRangeUser(0.8,1.2);  
  // h_2016_1->GetYaxis()->SetRangeUser(0,6);

  TCanvas* c = new TCanvas("c","c",600,600);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();
  h_2016_1->SetLineColor(kRed);
  h_2016_1->SetLineWidth(2); 
  h_2016_1->Draw();
  h_2017_1->SetLineColor(kBlue+1);
  h_2017_1->SetLineWidth(2);
  h_2017_1->Draw("SAME");
  h_2018_1->SetLineColor(kGreen+2);
  h_2018_1->SetLineWidth(2);
  h_2018_1->Draw("SAME");

  h_2016_pre_1->SetLineColor(kRed);
  h_2016_pre_1->SetLineWidth(2); 
  h_2016_pre_1->SetLineStyle(2);
  h_2016_pre_1->Draw("SAME");
  h_2017_pre_1->SetLineColor(kBlue+1);
  h_2017_pre_1->SetLineWidth(2);
  h_2017_pre_1->SetLineStyle(2);
  h_2017_pre_1->Draw("SAME");
  h_2018_pre_1->SetLineColor(kGreen+2);
  h_2018_pre_1->SetLineWidth(2);
  h_2018_pre_1->SetLineStyle(2);
  h_2018_pre_1->Draw("SAME");

  TLegend* leg = new TLegend(0.6,0.57,0.95,0.92);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_2016_1, "2016", "l");
  leg->AddEntry(h_2017_1, "2017", "l");
  leg->AddEntry(h_2018_1, "2018", "l");
  leg->AddEntry(h_2016_pre_1, "2016 raw", "l");
  leg->AddEntry(h_2017_pre_1, "2017 raw", "l");
  leg->AddEntry(h_2018_pre_1, "2018 raw", "l");
  leg->Draw();

  c->Print("hotvr_resolution.eps");
  c->Print("hotvr_resolution-eps-converted-to.pdf");

}
