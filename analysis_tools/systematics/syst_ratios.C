#include "../cosmetics.h"
#include "../plot_styles.C"

TH1F* get_hist_from_file(const TString file_name, const TString hist_name);

void syst_ratios(){
  TString year = "all";
  TString in_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Muon/"+year+"/";
  const TString prefix = "uhh2.AnalysisModuleRunner.MC.";
  
  const TString num_name = "SCALE";
  const TString den_name = "NOMINAL";

  const TString hist_name = "1btag1toptag20chi2_reco/Bstar_reco_M_rebin";

  // get histogram of numerator
  TH1F *h_num_up = get_hist_from_file(in_dir+num_name+"_up/"+prefix+"ST.root", hist_name);
  TH1F *h_num_down = get_hist_from_file(in_dir+num_name+"_down/"+prefix+"ST.root", hist_name);
  TH1F *h_den = get_hist_from_file(in_dir+den_name+"/"+prefix+"ST.root", hist_name);

  TGraphAsymmErrors* ratio_up = new TGraphAsymmErrors();
  ratio_up->Divide(h_num_up, h_den, "pois");
  ratio_up->SetMarkerColor(kRed);
  ratio_up->SetLineColor(kRed);
  ratio_up->GetHistogram()->GetYaxis()->SetRangeUser(0.0, 2.0);

  TGraphAsymmErrors* ratio_down = new TGraphAsymmErrors();
  ratio_down->Divide(h_num_down, h_den, "pois");
  ratio_down->SetMarkerColor(kBlue);
  ratio_down->SetLineColor(kBlue);

  setCMSStyle();
  TCanvas* c = new TCanvas("c","c",600,600);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();
  
  ratio_up->Draw("AP");
  ratio_down->Draw("SAME P");

  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(ratio_up, num_name+" up", "l");
  leg->AddEntry(ratio_down, num_name+" down", "l");
  leg->Draw();


  c->Print("ratio_"+num_name+"_"+den_name+"_"+year+".eps");
  c->Print("ratio_"+num_name+"_"+den_name+"_"+year+".png");

}


TH1F* get_hist_from_file(const TString file_name, const TString hist_name) {
  TFile *f = new TFile(file_name);
  TH1F *hist = (TH1F*)( f->Get(hist_name))->Clone();
  hist->SetDirectory(0);
  f->Close();
  delete f;
  return hist;
}
