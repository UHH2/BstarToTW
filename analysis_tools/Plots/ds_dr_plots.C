#include "../cosmetics.h"
#include "../plot_styles.C"

void reweight_hist(TH1F* h, TH1F* h_up, TH1F* h_down, TF1* ffit);
void create_output(const TString fout_name, const TString subdir_name, TH1F* hist);
void ds_dr_plots() {

  // const TString channel = "Electron";
  const TString channel = "Muon";
  
  TString in_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/"+channel+"/2016/Working/";
  TString in_dir_all = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/"+channel+"/all/NOMINAL/";
  const TString prefix = "uhh2.AnalysisModuleRunner.";

  const TString subdir_name= "1btag1toptag20chi2_reco";
  // const TString subdir_name= "2btag1toptag_reco";

  const TString hist_name = "Bstar_reco_M_rebin";
  TFile *f_ds = new TFile(in_dir+prefix+"MC.DS_ST.root");
  TFile *f_dr = new TFile(in_dir+prefix+"MC.DR_ST.root");

  TH1F* h_ds = (TH1F*) f_ds->Get(subdir_name+"/"+hist_name);
  TH1F* h_dr = (TH1F*) f_dr->Get(subdir_name+"/"+hist_name);

  // TH1F* h_ds = (TH1F*) f_ds->Get("2btag1toptag_reco/Bstar_reco_M_rebin");
  // TH1F* h_dr = (TH1F*) f_dr->Get("2btag1toptag_reco/Bstar_reco_M_rebin");

  h_ds->SetFillColor(796);
  h_ds->Scale(1,"width");
  h_dr->SetLineColor(810);
  h_dr->SetMarkerStyle(8);
  h_dr->SetMarkerColor(810);
  h_dr->Scale(1,"width");

  setCMSStyle();
  
  TGraphAsymmErrors *ratio_plot = new TGraphAsymmErrors(h_ds, h_dr, "pois");

  TCanvas* c = new TCanvas("c","c",600,600);
  // auto ratio_plot = new TRatioPlot(hs,h_data);
  TPad* pad_top = SetupRatioPadTop();
  TPad* pad_bot = SetupRatioPad();
  pad_top->Draw();
  pad_bot->Draw();
  pad_top->cd();
  gPad->SetLogy(true);
  h_ds->SetMaximum(max(h_ds->GetMaximum(), h_dr->GetMaximum()) * 3.0);
  h_ds->GetYaxis()->SetTitleSize(0.05 / 0.65);
  h_ds->GetYaxis()->SetLabelSize(0.04 / 0.65);
  h_ds->GetXaxis()->SetTitleSize(0.05 / 0.65);
  h_ds->GetXaxis()->SetLabelSize(0.04 / 0.65);
  h_ds->Draw("HIST");
  h_dr->Draw("SAME LPE");

  TLegend* leg = new TLegend(0.6,0.5,0.95,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_ds, "single top DS", "f");
  leg->AddEntry(h_dr, "single top DR", "ple");
  leg->Draw();

  pad_bot->cd();
  ratio_plot->GetXaxis()->SetLimits(h_ds->GetXaxis()->GetXmin(), h_ds->GetXaxis()->GetXmax());
  double xmin = ratio_plot->GetXaxis()->GetXmin();
  double xmax = ratio_plot->GetXaxis()->GetXmax();
  TLine *line = new TLine(xmin, 1, xmax, 1);
  line->SetLineColor(kBlack);
  ratio_plot->GetYaxis()->CenterTitle();
  ratio_plot->GetYaxis()->SetRangeUser(0.3, 1.7);
  ratio_plot->GetYaxis()->SetTitle("DS/DR");
  ratio_plot->GetYaxis()->SetTitleSize(0.05 / 0.34);
  ratio_plot->GetYaxis()->SetTitleOffset(1.1 * 0.34);
  ratio_plot->GetYaxis()->SetLabelSize(0.04 / 0.34);
  ratio_plot->GetXaxis()->SetTitleSize(0.05 / 0.34);
  ratio_plot->GetXaxis()->SetLabelSize(0.04 / 0.34);
  ratio_plot->GetXaxis()->SetTitle("M_{tW} [GeV]");
  ratio_plot->GetYaxis()->SetNdivisions(505);
  ratio_plot->Draw("AP0");
  line->Draw();

  c->Print(subdir_name+".eps");
  c->Print(subdir_name+".png");

  TCanvas* cfit = new TCanvas("cfit","cfit",600,600);
  TGraphAsymmErrors *ratio_plot_fit = (TGraphAsymmErrors*) ratio_plot->Clone();
  ratio_plot_fit->GetYaxis()->SetRangeUser(-0.1,2.4);
  ratio_plot_fit->GetYaxis()->SetTitleSize(0.05);
  ratio_plot_fit->GetYaxis()->SetTitleOffset(1.1);
  ratio_plot_fit->GetYaxis()->SetLabelSize(0.04);
  ratio_plot_fit->GetXaxis()->SetTitleSize(0.05);
  ratio_plot_fit->GetXaxis()->SetLabelSize(0.04);
  TPad* pad = SetupPad();
  ratio_plot_fit->Draw("AP0");
  TF1* fitfun = new TF1("fitfun","expo");
  ratio_plot_fit->Fit(fitfun);

  TFile* f_all = new TFile(in_dir_all+prefix+"MC.ST.root");
  TH1F* h_all= (TH1F*) f_all->Get(subdir_name+"/"+hist_name);
  TH1F* h_up = (TH1F*) h_all->Clone();
  h_up->Reset();
  TH1F* h_down = (TH1F*) h_all->Clone();
  h_down->Reset();
  reweight_hist(h_all, h_up, h_down, fitfun);
  create_output(prefix+"UP"+"MC.ST.root", subdir_name, h_up);
  create_output(prefix+"DOWN"+"MC.ST.root", subdir_name, h_down);
  
  cfit->Print(subdir_name+"_fit.eps");
  cfit->Print(subdir_name+"_fit.png");
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


void reweight_hist(TH1F* h, TH1F* h_up, TH1F* h_down, TF1* ffit) {
  std::cout << "reweight hist ..." << std::endl;
    
  for (int i = 0; i <= h->GetNbinsX(); ++i)
    {
      h_down->SetBinContent(i, h->GetBinContent(i) * ffit->Eval(h->GetBinCenter(i)));
      h_up->SetBinContent(i, h->GetBinContent(i) + abs(h_down->GetBinContent(i) - h->GetBinContent(i)));
    } 
}
