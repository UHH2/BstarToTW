#include "../../cosmetics.h"


void draw_covariance_matrix() {

  TFile *f = new TFile("fitDiagnostics.root");
  TH2D *h_cov = (TH2D*) f->Get("covariance_fit_b");
  h_cov->GetXaxis()->SetTitle("");
  h_cov->GetYaxis()->SetTitle("");

  h_cov->GetXaxis()->SetLabelSize(0.01);
  h_cov->GetYaxis()->SetLabelSize(0.01);

  h_cov->SetMarkerSize(0.25);
  gStyle->SetPaintTextFormat("1.1f");
  gStyle->SetOptStat("0");

  TCanvas *c = new TCanvas("c", "c", 1200, 1200);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();

  h_cov->Draw("COLZ TEXT");
  c->Print("cov.eps");

}
