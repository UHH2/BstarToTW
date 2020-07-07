#include "../cosmetics.h"
#include "../plot_styles.C"


void top_mass_plots() {

  const TString channel = "Electron";
  // const TString channel = "Muon";
  
  TString in_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/"+channel+"/all/Working/";
  const TString prefix = "uhh2.AnalysisModuleRunner.";


  

  TFile *f_tt = new TFile(in_dir+prefix+"MC.TT.root");
  TFile *f_st = new TFile(in_dir+prefix+"MC.ST.root");
  TFile *f_wjets = new TFile(in_dir+prefix+"MC.WJets.root");
  TFile *f_data = new TFile(in_dir+prefix+"DATA.DATA.root");

  TH2F *h2_tt_merged = (TH2F*) f_tt->Get("2btag_topmatch/M_Matched");
  TH2F *h2_tt_semi = (TH2F*) f_tt->Get("2btag_topmatch/M_Semi");
  TH2F *h2_tt_non = (TH2F*) f_tt->Get("2btag_topmatch/M_Non");

  TH2F *h2_st_merged = (TH2F*) f_st->Get("2btag_topmatch/M_Matched");
  TH2F *h2_st_semi = (TH2F*) f_st->Get("2btag_topmatch/M_Semi");
  TH2F *h2_st_non = (TH2F*) f_st->Get("2btag_topmatch/M_Non");

  TH2F *h2_wjets = (TH2F*) f_wjets->Get("2btag_topmatch/M_Non");

  TH2F *h2_data = (TH2F*) f_data->Get("2btag_topmatch/M_Non");

  setCMSStyle();

  for (int i = 1; i <= 4; ++i) 
    {

      TH1D *h_tt_merged = h2_tt_merged->ProjectionY("tt_merged_py", i, i);
      h_tt_merged->SetFillColor(810);
      TH1D *h_tt_semi   = h2_tt_semi->ProjectionY("tt_semi_py", i, i);
      h_tt_semi->SetFillColor(808);
      TH1D *h_tt_non    = h2_tt_non->ProjectionY("tt_non_py", i, i);
      h_tt_non->SetFillColor(806);

      TH1D *h_st_merged = h2_st_merged->ProjectionY("st_merged_py", i, i);
      h_st_merged->SetFillColor(800);
      TH1D *h_st_semi   = h2_st_semi->ProjectionY("st_semi_py", i, i);
      h_st_semi->SetFillColor(791);
      TH1D *h_st_non    = h2_st_non->ProjectionY("st_non_py", i, i);
      h_st_non->SetFillColor(794);
      TH1D *h_wjets = h2_wjets->ProjectionY("wjets_py", i, i);
      h_wjets->SetFillColor(594);

      THStack *hs = new THStack("hs","");
      hs->Add(h_wjets);
      hs->Add(h_st_non);
      hs->Add(h_st_semi);
      hs->Add(h_st_merged);
      hs->Add(h_tt_non);
      hs->Add(h_tt_semi);
      hs->Add(h_tt_merged);

      TH1D *h_data = h2_data->ProjectionY("data_py", i, i);
      h_data->SetMarkerColor(1);
      h_data->SetMarkerStyle(8);

      TList *stackHists = hs->GetHists();
  
      TH1* tmpHist = (TH1*)stackHists->At(0)->Clone();
      tmpHist->Reset();
 
      for (int i=0;i<stackHists->GetSize();++i) {
	tmpHist->Add((TH1*)stackHists->At(i));
      }

      TGraphAsymmErrors *ratio_plot = new TGraphAsymmErrors(h_data, tmpHist, "pois");

      TCanvas* c = new TCanvas("cfit","cfit",600,600);
      // auto ratio_plot = new TRatioPlot(hs,h_data);
      TPad* pad_top = SetupRatioPadTop();
      TPad* pad_bot = SetupRatioPad();
      pad_top->Draw();
      pad_bot->Draw();
      pad_top->cd();
      
      hs->SetMaximum(max(hs->GetMaximum(), h_data->GetMaximum())*1.1);
      hs->Draw("HIST");
      h_data->Draw("SAME PE");

      TLegend* leg = new TLegend(0.6,0.5,0.95,0.9);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      leg->AddEntry(h_data, "Data", "p");
      leg->AddEntry(h_tt_merged, "ttbar merged", "f");
      leg->AddEntry(h_tt_semi, "ttbar semi-merged", "f");
      leg->AddEntry(h_tt_non, "ttbar non-merged", "f");
      leg->AddEntry(h_st_merged, "single top merged", "f");
      leg->AddEntry(h_st_semi, "single top semi-merged", "f");
      leg->AddEntry(h_st_non, "single top non-merged", "f");
      leg->AddEntry(h_wjets, "W+Jets", "f");
      leg->Draw();

      pad_bot->cd();
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

      TString hist_name = "hist_"+ to_string(i);
      c->Print(hist_name+".eps");
      c->Print(hist_name+".png");
      
    }

}
