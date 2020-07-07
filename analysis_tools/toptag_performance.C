#include "TSystem.h"
#include <vector>

#include <iostream>

#include "cosmetics.h"

vector<char*> get_filenames(const char* ext)
{
  const char* inDir = "/nfs/dust/cms/user/albrecha/uhh2_102X/HOTVRStudiesOutput/ROC/root";
  char* dir = gSystem->ExpandPathName(inDir);
  void* dirp = gSystem->OpenDirectory(dir);
  const char* entry;
  vector<char*> filenames;
  TString str;
  while((entry = (char*)gSystem->GetDirEntry(dirp))) {
    str = entry;
    if(str.Contains(ext))
      filenames.push_back(gSystem->ConcatFileName(dir, entry));
  }
  return filenames;
}

//commented out soft drop stuff and deltaR_max
void toptag_performance()
{
  Int_t n_ptbins = 7;
  Double_t pt_xbins[9] = {0, 200, 400, 600, 800, 1000, 1500, 2000};
  int n_points = 100;

  TH1F* a_ttbar_hotvr = new TH1F("A_ttbar_HOTVR", "P_{T}", n_ptbins, pt_xbins);
  vector<TH1F*> b_ttbar_hotvr;
  for (int i = 0; i < n_points; ++i)
    {
      string name = "B_ttbar_HOTVR_" + to_string(i);
      b_ttbar_hotvr.push_back(new TH1F(name.c_str(), "P_{T}", n_ptbins, pt_xbins));
    }

  TH1F* a_QCD_hotvr = new TH1F("A_QCD_HOTVR", "P_{T}", n_ptbins, pt_xbins);
  vector<TH1F*> b_QCD_hotvr;
  for (int i = 0; i < n_points; ++i)
    {
      string name = "B_QCD_HOTVR_" + to_string(i);
      b_QCD_hotvr.push_back(new TH1F(name.c_str(), "P_{T}", n_ptbins, pt_xbins));
    }

 //  TH1F* a_ttbar_Softdrop = new TH1F("A_ttbar_Softdrop", "P_{T}", n_ptbins, pt_xbins);
 //  vector<TH1F*> b_ttbar_Softdrop;
 //  for (int i = 0; i < n_points; ++i)
 //    {
 //      string name = "B_ttbar_Softdrop_" + to_string(i);
 //      b_ttbar_Softdrop.push_back(new TH1F(name.c_str(), "P_{T}", n_ptbins, pt_xbins));
 //    }
 //
 //  TH1F* a_QCD_Softdrop = new TH1F("A_QCD_Softdrop", "P_{T}", n_ptbins, pt_xbins);
 // vector<TH1F*> b_QCD_Softdrop;
 //  for (int i = 0; i < n_points; ++i)
 //    {
 //      string name = "B_QCD_Softdrop_" + to_string(i);
 //      b_QCD_Softdrop.push_back(new TH1F(name.c_str(), "P_{T}", n_ptbins, pt_xbins));
 //    }

  vector<char*> filenames = get_filenames(".root");

  for (int i = 0; i < filenames.size(); ++i) //loop over root files
    {
      TFile *f = new TFile(filenames.at(i));
      TString filename = filenames.at(i);
      if (filename.Contains(".ttbar")) // fill ttbar
	{
	  a_ttbar_hotvr->Add((TH1F*)f->Get("HOTVR_Pre/pt_gentop"));
	 // a_ttbar_Softdrop->Add((TH1F*)f->Get("Softdrop_Pre/pt_gen"));
	  for (int i = 0; i < n_points; ++i)
	    {
	      std::string name_hotvr = "HOTVR_Performance_" + to_string(i) + "/pt_top_matched";
	     // std::string name_Softdrop = "Softdrop_Performance_" + to_string(i) + "/pt_top_matched";
	      b_ttbar_hotvr.at(i)->Add((TH1F*)f->Get(name_hotvr.c_str()));
	     // b_ttbar_Softdrop.at(i)->Add((TH1F*)f->Get(name_Softdrop.c_str()));
	    }
	}
      if (filename.Contains(".QCD")) // fill QCD
	{
	  a_QCD_hotvr->Add((TH1F*)f->Get("HOTVR_Pre/pt_gen"));
	 // a_QCD_Softdrop->Add((TH1F*)f->Get("Softdrop_Pre/pt_gen"));
	  for (int i = 0; i < n_points; ++i)
	    {
	      std::string name_hotvr = "HOTVR_Performance_" + to_string(i) + "/pt_top_mismatched";
	    //  std::string name_Softdrop = "Softdrop_Performance_" + to_string(i) + "/pt_top_mismatched";
	      b_QCD_hotvr.at(i)->Add((TH1F*)f->Get(name_hotvr.c_str()));
	    //  b_QCD_Softdrop.at(i)->Add((TH1F*)f->Get(name_Softdrop.c_str()));
	    }
	}
    }

  gStyle->SetOptStat(0);

  TCanvas *c = new TCanvas("c", "c", 600, 600);

  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();
  gPad->SetTicks();
  for (int bin = 1; bin < n_ptbins+1; ++bin){

    double N_ttbar_HOTVR = a_ttbar_hotvr->GetBinContent(bin);
    cout << N_ttbar_HOTVR << endl;
  //  double N_ttbar_Softdrop = a_ttbar_Softdrop->GetBinContent(bin);

    double N_QCD_HOTVR = a_QCD_hotvr->GetBinContent(bin);
    cout << N_QCD_HOTVR << endl;
  //  double N_QCD_Softdrop = a_QCD_Softdrop->GetBinContent(bin);

    TGraph *HOTVR = new TGraph(n_points);
  //  TGraph *Softdrop = new TGraph(n_points);

    for (int i = 0; i < n_points; ++i)
      {

	       double x_hotvr = b_ttbar_hotvr.at(i)->GetBinContent(bin) / (double)N_ttbar_HOTVR;
	       
         //	double x_Softdrop = b_ttbar_Softdrop.at(i)->GetBinContent(bin) / N_ttbar_Softdrop;

	       double y_hotvr = b_QCD_hotvr.at(i)->GetBinContent(bin) / (double)N_QCD_HOTVR;
	       
         //	double y_Softdrop = b_QCD_Softdrop.at(i)->GetBinContent(bin) / N_QCD_Softdrop;

	       HOTVR->SetPoint(i, x_hotvr, y_hotvr);
         //	Softdrop->SetPoint(i, x_Softdrop, y_Softdrop);
      }

    HOTVR->SetLineWidth(2);
    HOTVR->SetLineColor(kRed);

  //  Softdrop->SetLineWidth(2);
  //  Softdrop->SetLineColor(kBlue);

    // TGraph *Softdrop_wp1 = new TGraph(1);
    // Softdrop_wp1->SetPoint(1,Softdrop->GetX()[19], Softdrop->GetY()[19]);
    // Softdrop_wp1->SetMarkerStyle(22); Softdrop_wp1->SetMarkerColor(1);
    // TGraph *Softdrop_wp2 = new TGraph(1);
    // Softdrop_wp2->SetPoint(1,Softdrop->GetX()[33], Softdrop->GetY()[33]);
    // Softdrop_wp2->SetMarkerStyle(22); Softdrop_wp2->SetMarkerColor(2);
    // TGraph *Softdrop_wp3 = new TGraph(1);
    // Softdrop_wp3->SetPoint(1,Softdrop->GetX()[43], Softdrop->GetY()[43]);
    // Softdrop_wp3->SetMarkerStyle(22); Softdrop_wp3->SetMarkerColor(3);
    // TGraph *Softdrop_wp4 = new TGraph(1);
    // Softdrop_wp4->SetPoint(1,Softdrop->GetX()[50], Softdrop->GetY()[50]);
    // Softdrop_wp4->SetMarkerStyle(22); Softdrop_wp4->SetMarkerColor(6);

    // TGraph *hotvr_wp = new TGraph(1);
    // hotvr_wp->SetPoint(1,HOTVR->GetX()[50], HOTVR->GetY()[50]);
    // hotvr_wp->SetMarkerStyle(23); hotvr_wp->SetMarkerColor(1);

    TMultiGraph *mg = new TMultiGraph();
    mg->Add(HOTVR);
    // mg->Add(Softdrop);
    // mg->Add(Softdrop_wp1, "P");
    // mg->Add(Softdrop_wp2, "P");
    // mg->Add(Softdrop_wp3, "P");
    // mg->Add(Softdrop_wp4, "P");
    // mg->Add(hotvr_wp, "P");
    mg->Draw("AL");

    gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn;
    mg->SetTitle();
    Hist_Cosmetics(mg, false);
    mg->GetXaxis()->SetTitle("#varepsilon_{S}");
    mg->GetYaxis()->SetTitle("#varepsilon_{B}");
    mg->SetMaximum(0.1);
    mg->SetMinimum(1e-4);
    mg->GetXaxis()->SetLimits(0,0.8);
    gPad->Modified(); gPad->Update(); // make sure it's really (re)drawn;


    std::string interval;
    interval = to_string((int)pt_xbins[bin-1]) + "GeV < p_{T} < " + to_string((int)pt_xbins[bin]) + "GeV";

    TLegend *leg = new TLegend(0.4,0.19,0.7,0.34);
    leg->SetFillColor(0);
    leg->SetLineColor(1);
    leg->SetBorderSize(0);
    leg->SetTextSize(.04);
    leg->SetFillStyle(0);
    leg->SetHeader(interval.c_str());
    leg->AddEntry(HOTVR, "HOTVR TopTagger", "l");
  //  leg->AddEntry(Softdrop, "SoftDrop TopTagger", "l");
    leg->Draw();

    TLegend *leg_wp = new TLegend(0.55,0.36,0.8,0.7);
    leg_wp->SetFillColor(0);
    leg_wp->SetLineColor(1);
    leg_wp->SetBorderSize(0);
    leg_wp->SetTextSize(.03);
    leg_wp->SetFillStyle(0);
    // leg_wp->AddEntry((TObject*)0, "SoftDrop TopTagger:", "");
    // leg_wp->AddEntry(Softdrop_wp1, "WP: #tau_{3/2}<0.81", "p");
    // leg_wp->AddEntry(Softdrop_wp2, "WP: #tau_{3/2}<0.67", "p");
    // leg_wp->AddEntry(Softdrop_wp3, "WP: #tau_{3/2}<0.57", "p");
    // leg_wp->AddEntry(Softdrop_wp4, "WP: #tau_{3/2}<0.50", "p");
    leg_wp->AddEntry((TObject*)0, "HOTVR TopTagger:", "");
    // leg_wp->AddEntry(hotvr_wp,  "WP: #tau_{3/2}<0.56", "p");
    leg_wp->Draw();

    pad->SetLogy();
    std::string name = "Performance_PtBin_" + to_string(bin) + ".eps";
    c->Print(name.c_str());
  }
}
