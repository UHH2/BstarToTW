#include "../cosmetics.h"
#include "../plot_styles.C"

TH1F* get_hist_from_file(const TString file_name, const TString hist_name);

void shape_ratios() {

  TString f_num = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/analysis_tools/CombineTools/control_region_only/fitDiagnostics.root";
  TString n_num = "shapes_fit_b/muon_region_2b1t/TT";

  TString f_den = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/analysis_tools/CombineTools/control_region_only/fitDiagnostics.root";
  TString n_den = "shapes_prefit/muon_region_2b1t/TT";

  // TString f_num = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Muon/all/TOPPT_beta_down/uhh2.AnalysisModuleRunner.MC.TT.root";
  // TString n_num = "2btag1toptag_reco/Bstar_reco_M_rebin";

  // TString f_den = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Muon/all/NOMINAL/uhh2.AnalysisModuleRunner.MC.TT.root";
  // TString n_den = "2btag1toptag_reco/Bstar_reco_M_rebin";
  
  TH1F* h_num = get_hist_from_file(f_num, n_num);
  TH1F* h_den = get_hist_from_file(f_den, n_den);

  h_den->SetFillColor(796);
  h_den->Scale(1,"width");
  h_num->SetLineColor(810);
  h_num->SetMarkerStyle(8);
  h_num->SetMarkerColor(810);
  h_num->Scale(1,"width");

  setCMSStyle();
  
  TGraphAsymmErrors *ratio_plot = new TGraphAsymmErrors(h_num, h_den, "pois");

  TCanvas* c = new TCanvas("c","c",600,600);
  // auto ratio_plot = new TRatioPlot(hs,h_data);
  TPad* pad_top = SetupRatioPadTop();
  TPad* pad_bot = SetupRatioPad();
  pad_top->Draw();
  pad_bot->Draw();
  pad_top->cd();
  gPad->SetLogy(true);
  h_den->SetMaximum(max(h_den->GetMaximum(), h_num->GetMaximum()) * 3.0);
  h_den->GetYaxis()->SetTitleSize(0.05 / 0.65);
  h_den->GetYaxis()->SetLabelSize(0.04 / 0.65);
  h_num->GetXaxis()->SetTitleSize(0.05 / 0.65);
  h_num->GetXaxis()->SetLabelSize(0.04 / 0.65);
  h_den->Draw("HIST");
  h_num->Draw("SAME LPE");

  TLegend* leg = new TLegend(0.6,0.5,0.95,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_num, n_num, "f");
  leg->AddEntry(h_den, n_den, "ple");
  leg->Draw();

  pad_bot->cd();
  ratio_plot->GetXaxis()->SetLimits(h_den->GetXaxis()->GetXmin(), h_den->GetXaxis()->GetXmax());
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

  c->Print("ratio.eps");
  c->Print("ratio.png");
}

TH1F* get_hist_from_file(const TString file_name, const TString hist_name) {
  TFile *f = new TFile(file_name);
  TH1F *hist = (TH1F*)( f->Get(hist_name))->Clone();
  hist->SetDirectory(0);
  f->Close();
  delete f;
  return hist;
}
