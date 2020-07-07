#include "cosmetics.h"
#include "plot_styles.C"

TH1F* get_hist_from_file(const TString file_name, const TString hist_name);

void compare_hist(const TString channel, const TString year, const TString process, const TString hist_name, const TString region1, const TString  region2) {

  TString in_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/"+channel+"/"+year+"/";
  const TString prefix = "uhh2.AnalysisModuleRunner.MC.";
  
  TH1F *h_region1 = get_hist_from_file(in_dir+region1+"/"+prefix+process+".root", hist_name);
  h_region1->SetName("region1");
  h_region1->SetLineColor(kBlue);
  h_region1->SetLineWidth(2);
  h_region1->SetMarkerStyle(0);
  h_region1->GetXaxis()->SetTitle(h_region1->GetTitle());
  TH1F *h_region2 = get_hist_from_file(in_dir+region2+"/"+prefix+process+".root", hist_name);
  h_region2->SetName("region2");
  h_region2->SetLineColor(kRed);
  h_region2->SetLineWidth(2);
  h_region2->SetMarkerStyle(0);
  
  setCMSStyle();
  TCanvas* c = new TCanvas("cfit","cfit",600,600);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();
  h_region1->Draw("HIST");
  h_region2->Draw("HIST SAME");

  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_region1, region1, "l");
  leg->AddEntry(h_region2, region2, "l");
  leg->Draw();

  c->Print("hists.eps");
  gPad->SetLogy(true);
  c->Print("hists_logy.eps");
}

TH1F* get_hist_from_file(const TString file_name, const TString hist_name) {
  TFile *f = new TFile(file_name);
  TH1F *hist = (TH1F*)( f->Get(hist_name))->Clone();
  hist->SetDirectory(0);
  f->Close();
  delete f;
  return hist;
}
