#include "TSystem.h"
#include <vector>

#include <iostream>
#include "TMath.h"

#include "cosmetics.h"

bool my_cmp(const char *c1, const char *c2)
{
    return strcmp(c1, c2) < 0;
}	

vector<char*> get_filenames(const char* ext)
{
  const char* inDir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Muon/all/NOMINAL";
  char* dir = gSystem->ExpandPathName(inDir);
  void* dirp = gSystem->OpenDirectory(dir);
  const char* entry;
  vector<char*> filenames;
  TString str;
  while((entry = (char*)gSystem->GetDirEntry(dirp))) {
    str = entry;
    if( str.Contains(ext) && !str.Contains("LH"))
      filenames.push_back(gSystem->ConcatFileName(dir, entry));
  }
  std::sort(filenames.begin(), filenames.end(), my_cmp );
  return filenames;
}

void fitTopMass()
{

  Double_t width = 400;
  Double_t height = 400;
  TCanvas * c1 = new TCanvas("c", "c", width, height);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();

  vector<char*> filenames = get_filenames("MC.BstarToTW");
  TH1F *combined = new TH1F("Top_reco_M", "M", 30, 0, 300);
  for (int i = 1; i < filenames.size(); ++i)
    {
      TFile *f = new TFile(filenames.at(i));
      combined->Add((TH1F*)f->Get("BstarToTWMatchedReco/Top_reco_M"));
    }
  
  TF1 *fitg = new TF1("fitg", "gaus", 150, 210);

  Hist_Cosmetics(combined, false);
  combined->GetXaxis()->SetTitle("M_{top} [GeV/c^{2}]");
  combined->GetYaxis()->SetTitle("Events");
  combined->Fit("fitg", "R");

  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(1);
  TPaveStats *st = (TPaveStats*)combined->GetListOfFunctions()->FindObject("stats");
  st->SetX1NDC(0.6); //new x start position
  st->SetX2NDC(0.92); //new x end position

  c1->Print("mTopFit.eps");
}


void fitPtGen()
{

  Double_t width = 400;
  Double_t height = 400;
  TCanvas * c1 = new TCanvas("c", "c", width, height);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();

  vector<char*> filenames = get_filenames("MC.BstarToTW");
  TH1F *combined = new TH1F("pt_rel_dev", "M", 40, -0.2, 0.2);
  for (int i = 1; i < filenames.size(); ++i)
    {
      TFile *f = new TFile(filenames.at(i));
      combined->Add((TH1F*)f->Get("PreSel_HOTVRPerformance/delta_pt_gen_reco"));
    }
  
  TF1 *fitg = new TF1("fitg", "gaus", -0.08, 0.1);

  Hist_Cosmetics(combined, false);
  combined->GetXaxis()->SetTitle("p_{T,topjet} - p_{T,gentop} / p_{T,gentop}");
  combined->GetYaxis()->SetTitle("Events");
  // combined->Fit("fitg", "R");
  combined->Draw();
  // gStyle->SetOptStat(0); 
  // gStyle->SetOptFit(1);
  // TPaveStats *st = (TPaveStats*)combined->GetListOfFunctions()->FindObject("stats");
  // st->SetX1NDC(0.6); //new x start position
  // st->SetX2NDC(0.92); //new x end position

  c1->Print("ptGenFit.eps");
}



Double_t lorentzian(Double_t *x, Double_t *par) {

  return par[2] * (par[1]*par[1])/((x[0]-par[0])*(x[0]-par[0])+(par[1]*par[1]));
}

void fitDeltaPhi()
{

  Double_t width = 800;
  Double_t height = 800;
  TCanvas * c1 = new TCanvas("c", "c", width, height);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();
  
  vector<char*> filenames = get_filenames(".BstarToTW");
  TH1F *combined = new TH1F("DeltaPhi_top_W", "Delta Phi", 50, 2.5, 3.5);
  for (int i = 1; i < filenames.size(); ++i)
    {
      TFile *f = new TFile(filenames.at(i));
      combined->Add((TH1F*)f->Get("BstarToTWMatchedReco/DeltaPhi_top_W"));
    }


  TF1 *lorentz = new TF1("lorentz", lorentzian, 2.6, 3.7, 3);
  lorentz->SetParLimits(0,2,4);
  lorentz->SetParLimits(1,1e-6,1);

  TF1 *fitg = new TF1("fitg", "gaus", 3.06, 3.22);
  // fitg->SetParLimits(2, 1e-10,0.1);

  Hist_Cosmetics(combined, false);
  combined->GetXaxis()->SetTitle("#Delta #phi_{top,W}");
  combined->GetYaxis()->SetTitle("Events");
  // combined->Fit("lorentz", "R");
  combined->Fit("fitg", "R");

  gPad->SetTicks();
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(1);
  TPaveStats *st = (TPaveStats*)combined->GetListOfFunctions()->FindObject("stats");
  st->SetX1NDC(0.25); //new x start position
  st->SetX2NDC(0.57); //new x end position

  c1->Print("deltaPhiFit.eps");
}

void fitDeltaPtOver()
{

  Double_t width = 400;
  Double_t height = 400;
  TCanvas * c1 = new TCanvas("c", "c", width, height);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();
  
  vector<char*> filenames = get_filenames(".BstarToTW");
  TH1F *combined = new TH1F("DeltaPt_top_W", "Delta Pt", 50, -0.5, 0.5);
  for (int i = 1; i < filenames.size(); ++i)
    {
      TFile *f = new TFile(filenames.at(i));
      combined->Add((TH1F*)f->Get("BstarToTWMatchedReco/DeltaPt_top_W_over_pt"));
    }
  
  // TF1 *lorentz = new TF1("lorentz", lorentzian, -0.15, 0.1, 3);
  // lorentz->SetParLimits(0,-1,1);
  // lorentz->SetParLimits(1,1e-10,0.1);
  TF1 *fitg = new TF1("fitg", "gaus", -0.09, 0.09);

  Hist_Cosmetics(combined, false);
  combined->GetXaxis()->SetTitle("#Delta p_{T top,W} / p_{T, top}");
  combined->GetYaxis()->SetTitle("Events");
  combined->Fit("fitg", "R");

  gPad->SetTicks();  

  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(1);
  TPaveStats *st = (TPaveStats*)combined->GetListOfFunctions()->FindObject("stats");
  st->SetX1NDC(0.63); //new x start position
  st->SetX2NDC(0.95); //new x end position

  c1->Print("deltaPtFit.eps");
}

void fitDeltaPt()
{

  Double_t width = 800;
  Double_t height = 800;
  TCanvas * c1 = new TCanvas("c", "c", width, height);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();
  
  vector<char*> filenames = get_filenames(".BstarToTW");
  TH1F *combined = new TH1F("DeltaPt_top_W", "p_{T balance}", 100, -1, 1);
  for (int i = 1; i < filenames.size(); ++i)
    {
      TFile *f = new TFile(filenames.at(i));
      combined->Add((TH1F*)f->Get("BstarToTWMatchedReco/PtBalance"));
    }
  
  TF1 *fitg = new TF1("fitg", "gaus", -0.08, 0.08);
  TF1 *fitl = new TF1("fitl", lorentzian, -1, 1, 3);

  combined->GetXaxis()->SetTitle("p_{T, balance}");
  combined->GetYaxis()->SetTitle("Events");
  combined->Fit("fitg", "R");
  // combined->Fit("fitl", "R");
  
  gPad->SetTicks();
  gStyle->SetOptStat(0); 
  gStyle->SetOptFit(1);
  // TPaveStats *st = (TPaveStats*)combined->GetListOfFunctions()->FindObject("stats");
  // st->SetX1NDC(0.6); //new x start position
  // st->SetX2NDC(0.92); //new x end position

  c1->Print("deltaPtFit.eps");
}
