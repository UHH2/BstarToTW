#include "../cosmetics.h"
#include "../plot_styles.C"


void postFitPlot() {

  
  const TString file_name = "Postfit.root";
  // const TString channel = "muon2b1t_postfit";
  // const TString channel = "electron2b1t_postfit";

  const TString channel = "muon2b1t_prefit";
  // const TString channel = "electron2b1t_prefit";


  const int nbins = 20;
  double xbins[] = {500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 4000}; // new extended binning

  TFile *fin = new TFile(file_name);

  TH1F *h_data = (TH1F*) fin->Get(channel + "/data_obs");
  h_data->SetBins(nbins, xbins);
  TH1F *h_tt = (TH1F*) fin->Get(channel + "/TT");
  h_tt->SetBins(nbins, xbins);
  TH1F *h_st = (TH1F*) fin->Get(channel + "/ST");
  h_st->SetBins(nbins, xbins);
  TH1F *h_other = (TH1F*) fin->Get(channel + "/Other");
  h_other->SetBins(nbins, xbins);
  TH1F *h_err =  (TH1F*) fin->Get(channel + "/TotalBkg");
  h_err->SetBins(nbins, xbins);

  setCMSStyle();

  THStack *h_bkg = new THStack("hs","");

  h_data->SetLineColor(1);
  h_data->SetMarkerStyle(8);
  h_tt->SetFillColor(810);
  h_st->SetFillColor(800);
  h_other->SetFillColor(594);

  h_bkg->Add(h_other);
  h_bkg->Add(h_st);
  h_bkg->Add(h_tt);

  TGraphAsymmErrors *ratio_plot = new TGraphAsymmErrors(h_data, h_err, "pois");

  TCanvas* c = new TCanvas("cfit","cfit",600,600);
  TPad* pad_top = SetupRatioPadTop();
  TPad* pad_bot = SetupRatioPad();
  pad_top->Draw();
  pad_bot->Draw();
  pad_top->cd();
  pad_top->SetLogy(true);

  h_bkg->SetMaximum(max(h_bkg->GetMaximum(), h_data->GetMaximum())*1.1);
  h_err->SetFillColorAlpha(12, 0.3);
  h_err->SetMarkerSize(0);
  h_bkg->Draw("HIST");
  h_data->Draw("PESAME");
  h_err->Draw("PSAME");

  TLegend* leg = new TLegend(0.6,0.5,0.95,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(h_data, "Data", "p");
  leg->AddEntry(h_tt, "TT", "f");
  leg->AddEntry(h_st, "ST", "f");
  leg->AddEntry(h_other, "Other", "f");
  leg->Draw();

  pad_bot->cd();
  pad_bot->SetLogy(false);

  ratio_plot->GetXaxis()->SetLimits(h_data->GetXaxis()->GetXmin(), h_data->GetXaxis()->GetXmax());
  double xmin = ratio_plot->GetXaxis()->GetXmin();
  double xmax = ratio_plot->GetXaxis()->GetXmax();
  TLine *line = new TLine(xmin, 1, xmax, 1);
  line->SetLineColor(kBlack);
  ratio_plot->GetYaxis()->CenterTitle();
  ratio_plot->GetYaxis()->SetRangeUser(0.3, 1.7);
  ratio_plot->GetYaxis()->SetTitle("Data/MC");
  ratio_plot->GetXaxis()->SetTitle("M_{tW} [GeV]");
  ratio_plot->GetYaxis()->SetNdivisions(505);
  ratio_plot->Draw("AP0");
  line->Draw();

  TString hist_name = "pre_fit_"+channel;
  c->Print(hist_name+".eps");
  c->Print(hist_name+".png");
      
}
