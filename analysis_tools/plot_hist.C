#include "cosmetics.h"
#include "plot_styles.C"

TH1F* get_hist_from_file(const TString file_name, const TString hist_name);

void plot_hist(const TString channel, const TString year, const TString region, const TString process, const TString hist_name) {

  TString in_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/"+channel+"/"+year+"/";
  const TString prefix = "uhh2.AnalysisModuleRunner.MC.";
  
  TH1F *h = get_hist_from_file(in_dir+region+"/"+prefix+process+".root", hist_name);
  h->SetLineColor(kBlue);
  h->SetLineWidth(2);
  h->SetMarkerStyle(20);
  h->GetXaxis()->SetTitle(h->GetTitle());
  
  setCMSStyle();
  TCanvas* c = new TCanvas("cfit","cfit",600,600);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();
  h->Draw();

  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h, h->GetName(), "p");
  leg->Draw();

  c->Print("hists.eps");
}

TH1F* get_hist_from_file(const TString file_name, const TString hist_name) {
  TFile *f = new TFile(file_name);
  TH1F *hist = (TH1F*)( f->Get(hist_name))->Clone();
  hist->SetDirectory(0);
  f->Close();
  delete f;
  return hist;
}
