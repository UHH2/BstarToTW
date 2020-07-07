#include "cosmetics.h"
#include "plot_styles.C"

TH1F* get_hist_from_file(const TString file_name, const TString hist_name);

void trigger_efficiency_plots(const TString & channel, const TString & year){
  // TString year = "all";
  TString in_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/"+channel+"/"+year+"/NOMINAL/";
  const TString prefix = "uhh2.AnalysisModuleRunner.";
  
  const TString num_name = "Trigger_Analysis";
  const TString den_name = "PreSel_Analysis";

  const TString hist_name = "lepton_pt";

  // get histogram of numerator
  TH1F *h_num_data = get_hist_from_file(in_dir+prefix+"DATA.DATA.root", num_name+"/"+hist_name);
  TH1F *h_den_data = get_hist_from_file(in_dir+prefix+"DATA.DATA.root", den_name+"/"+hist_name);
  TGraphAsymmErrors* ratio_data = new TGraphAsymmErrors();
  ratio_data->BayesDivide(h_num_data, h_den_data);
  ratio_data->SetMarkerColor(kBlue);
  ratio_data->SetLineColor(kBlue);

  TH1F *h_num_mc = get_hist_from_file(in_dir+prefix+"MC.TT.root", num_name+"/"+hist_name);
  TH1F *h_den_mc = get_hist_from_file(in_dir+prefix+"MC.TT.root", den_name+"/"+hist_name);
  TGraphAsymmErrors* ratio_mc = new TGraphAsymmErrors();
  ratio_mc->BayesDivide(h_num_mc, h_den_mc);
  ratio_mc->SetMarkerColor(kRed);
  ratio_mc->SetLineColor(kRed);

  setCMSStyle();
  TCanvas* c = new TCanvas("c","c",600,600);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();
  
  ratio_data->Draw("AP");
  ratio_mc->Draw("SAME P");

  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(ratio_data, "DATA", "l");
  leg->AddEntry(ratio_mc, "MC", "l");
  leg->Draw();


  c->Print("efficiency.eps");
  c->Print("efficiency.png");

}


TH1F* get_hist_from_file(const TString file_name, const TString hist_name) {
  TFile *f = new TFile(file_name);
  TH1F *hist = (TH1F*)( f->Get(hist_name))->Clone();
  hist->SetDirectory(0);
  f->Close();
  delete f;
  return hist;
}
