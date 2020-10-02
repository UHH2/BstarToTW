#include "../cosmetics.h"


void draw_covariance_matrix() {

  TFile *f = new TFile("fitDiagnostics.root");
  TH2D *h_cov = (TH2D*) f->Get("covariance_fit_b");
  TH2D *h_cov_reduced = (TH2D*) h_cov->Clone("covarianve_fit_b_reduced");
  h_cov_reduced->GetXaxis()->SetRange(1,48);
  h_cov_reduced->GetYaxis()->SetRange(h_cov_reduced->GetYaxis()->GetLast() - 47, h_cov_reduced->GetYaxis()->GetLast());

  h_cov->GetXaxis()->SetTitle("");
  h_cov->GetYaxis()->SetTitle("");

  h_cov->GetXaxis()->SetLabelSize(0.01);
  h_cov->GetYaxis()->SetLabelSize(0.01);

  h_cov_reduced->GetXaxis()->SetTitle("");
  h_cov_reduced->GetYaxis()->SetTitle("");

  h_cov_reduced->GetXaxis()->SetLabelSize(0.025);
  h_cov_reduced->GetYaxis()->SetLabelSize(0.025);

  h_cov->SetMarkerSize(0.25);
  h_cov_reduced->SetMarkerSize(0.5);

  gStyle->SetPaintTextFormat("1.1f");
  gStyle->SetOptStat("0");

  TCanvas *c = new TCanvas("c", "c", 1200, 1200);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();

  h_cov_reduced->Draw("COLZ TEXT");
  c->Print("cov.eps");

}
