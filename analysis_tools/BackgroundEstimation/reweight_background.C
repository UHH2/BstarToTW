#include "../cosmetics.h"
#include "../plot_styles.C"
#include <vector>
#include <cmath>

void fill_histogram_from_file(TH1F* hist, const TString &fName, const TString &region, const TString &histName, const int &weightSign = 1);
TGraphAsymmErrors* fit_ratio(TH1F* hSignal, TH1F* hControl);
void fill_reweighted(TH1F* hReweight, TF1* fitfun, const TString &fname, const int &weightSign = 1);
void fill_reweighted(TH1F* hReweight, TGraphAsymmErrors* ratio, const TString &fname, const int &weightSign = 1);
void draw_extrapolation(TH1F* hReweight, TH1F* hSignal);

// see math/mathcore/src/PdfFuncMathCore.cxx in ROOT 6.x
double crystalball_function(double x, double alpha, double n, double sigma, double mean) {
  // evaluate the crystal ball function
  if (sigma < 0.)     return 0.;
  double z = (x - mean)/sigma; 
  if (alpha < 0) z = -z; 
  double abs_alpha = std::abs(alpha);
  // double C = n/abs_alpha * 1./(n-1.) * std::exp(-alpha*alpha/2.);
  // double D = std::sqrt(M_PI/2.)*(1.+ROOT::Math::erf(abs_alpha/std::sqrt(2.)));
  // double N = 1./(sigma*(C+D));
  if (z  > - abs_alpha)
    return std::exp(- 0.5 * z * z);
  else {
    //double A = std::pow(n/abs_alpha,n) * std::exp(-0.5*abs_alpha*abs_alpha);
    double nDivAlpha = n/abs_alpha;
    double AA =  std::exp(-0.5*abs_alpha*abs_alpha);
    double B = nDivAlpha -abs_alpha;
    double arg = nDivAlpha/(B-z);
    return AA * std::pow(arg,n);
  }
}

double crystalball_function(const double *x, const double *p) {
  // if ((!x) || (!p)) return 0.; // just a precaution
  // [Constant] * ROOT::Math::crystalball_function(x, [Alpha], [N], [Sigma], [Mean])
  return (p[0] * crystalball_function(x[0], p[3], p[4], p[2], p[1]));
}

double landau_flattail(double x, double landau_scale, double landau_mpv, double landau_sigma, double erf_mean, double erf_sigma, double scale) {
  double L = landau_scale * TMath::Landau(x, landau_mpv, landau_sigma);
  double erf_arg = (x-erf_mean)/erf_sigma;
  double erf = (1 + TMath::Erf(erf_arg)) / 2;
  return ((1-erf) * L + scale*erf);
}

double landau_flattail(const double *x, const double *p) {
  return landau_flattail(x[0], p[0], p[1], p[2], p[3], p[4], p[5]);
}

void reweight_background(TString channel, const TString year = "2018", const TString signalRegion = "1btag1toptag20chi2") {

  const TString controlRegion = "0btag1toptag";

  vector<TString> inDirs;
  if (channel == "Combined")
    {
      inDirs.push_back("/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Electron/"+year+"/BACKGROUND/");
      inDirs.push_back("/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Muon/"+year+"/BACKGROUND/");
    }
  else
    {
      inDirs.push_back("/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/"+channel+"/"+year+"/BACKGROUND/");
    }
  const TString prefix = "uhh2.AnalysisModuleRunner.";
   
  const int nSamples = 4;
  const TString sampleNames[nSamples] = {"WJets", "DYJets", "Diboson", "QCD"};
  const int nTopSamples = 6 - nSamples;
  const TString topSampleNames[nTopSamples] = {"TT", "ST"};
  // read in histograms
  const int nbins = 17;
  const double xbins[nbins] = {0, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1700, 1900, 2100, 2400, 3500};
  const int nbins_rebin = 21;
  double xbins_rebin[nbins_rebin] = {500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3200};

  TH1F *hSignal_ST = new TH1F("signal_st", channel+"_"+signalRegion, nbins_rebin-1, xbins_rebin);
  TH1F *hControl_ST = new TH1F("control_st", channel+"_"+controlRegion, nbins_rebin-1, xbins_rebin);
  // TH1F *hSignal_ST2 = new TH1F("signal_st", channel+"_"+signalRegion, 50, 100, 3500);
  // TH1F *hControl_ST2 = new TH1F("control_st", channel+"_"+controlRegion, 50, 100, 3500);
  TH1F *hSignal_M = new TH1F("signal_m", channel+"_"+signalRegion, nbins-1, xbins);
  TH1F *hControl_M = new TH1F("control_m", channel+"_"+controlRegion, nbins-1, xbins);
  TH1F *hClosure_M = new TH1F("closure", channel+"_"+signalRegion+"_data", nbins-1, xbins);

  TH1F *hMCReweight_bin_M = new TH1F("MCreweight_bin", "_MC_reweight_bin", nbins-1, xbins);
  TH1F *hMCReweight_fit_M = new TH1F("MCreweight_fit", "_MC_reweight_fit", nbins-1, xbins);
  TH1F *hMCReweight_expo_M = new TH1F("MCreweight_expo", "_MC_reweight_expo", nbins-1, xbins);
  
  TH1F *hReweight_bin_M = new TH1F("reweight_bin", "_reweight_bin", nbins-1, xbins);
  TH1F *hReweight_fit_M = new TH1F("reweight_fit", "_reweight_fit", nbins-1, xbins);
  TH1F *hReweight_expo_M = new TH1F("reweight_expo", "_reweight_expo", nbins-1, xbins);

  hReweight_bin_M->Sumw2();
  hReweight_fit_M->Sumw2();
  hReweight_expo_M->Sumw2();
  vector<TString> fNames;
  for (TString inDir : inDirs)
    {
      for (int i = 0; i < nSamples; ++i)
	{
	  TString fName = inDir + prefix + "MC." + sampleNames[i] + ".root";
	  fNames.push_back(fName);
	  fill_histogram_from_file(hSignal_M, fName, signalRegion, "_reco/Bstar_reco_M");
	  fill_histogram_from_file(hControl_M, fName, controlRegion, "_reco/Bstar_reco_M");
	  if (sampleNames[i] == "QCD") continue;
	  // fill_histogram_from_file(hSignal_ST2, fName, signalRegion, "_reco/Bstar_reco_M_fine");
	  // fill_histogram_from_file(hControl_ST2, fName, controlRegion, "_reco/Bstar_reco_M_fine");
	  fill_histogram_from_file(hSignal_ST, fName, signalRegion, "_reco/Bstar_reco_M_rebin");
	  fill_histogram_from_file(hControl_ST, fName, controlRegion, "_reco/Bstar_reco_M_rebin");
	}
    }

  // fit_ratio
  TGraphAsymmErrors* ratio = fit_ratio(hSignal_ST, hControl_ST);
  
  // reweight control region
  TF1 *ffit = ratio->GetFunction("ffit");

  // hSignal_ST2->Draw();
  // TF1 *fexpo = ratio->GetFunction("fpol1");

  for (TString inDir : inDirs)
    {
      fill_reweighted(hReweight_bin_M, ratio, inDir + prefix+ "DATA.DATA.root");
      fill_reweighted(hReweight_fit_M, ffit, inDir + prefix+ "DATA.DATA.root");
      // fill_reweighted(hReweight_expo_M, fexpo, inDir + prefix+ "DATA.DATA.root");
      fill_histogram_from_file(hClosure_M, inDir + prefix+ "DATA.DATA.root", signalRegion, "_reco/Bstar_reco_M");
      for (int i = 0; i < nTopSamples; ++i)
  	{
  	  TString fName = inDir + prefix + "MC." + topSampleNames[i] + ".root";
  	  fill_reweighted(hReweight_bin_M, ratio, fName, -1);
  	  fill_reweighted(hReweight_fit_M, ffit, fName, -1);
  	  // fill_reweighted(hReweight_expo_M, fexpo, fName, -1);
  	  fill_histogram_from_file(hClosure_M, fName, signalRegion, "_reco/Bstar_reco_M", -1);
  	}
      for (int i = 0; i < nSamples; ++i)
  	{
  	  TString fName = inDir + prefix + "MC." + sampleNames[i] + ".root";
  	  fill_reweighted(hMCReweight_bin_M, ratio, fName);
  	  fill_reweighted(hMCReweight_fit_M, ffit, fName);
  	  // fill_reweighted(hMCReweight_expo_M, fexpo, fName);  	  
  	}
    }
  
  // create output file
  // hReweight_bin_M->Scale(1.,"width");
  // hMCReweight_bin_M->Scale(1.,"width");
  // hReweight_fit_M->Scale(1.,"width");
  // hMCReweight_fit_M->Scale(1.,"width");
  // hSignal_M->Scale(1.,"width");
  // hClosure_M->Scale(1.,"width");

  // hControl_M->Scale(1.,"width");
  // hControl_M->Draw();

  // draw_extrapolation(hReweight_bin_M, hSignal_M);
  // draw_extrapolation(hReweight_bin_M, hClosure_M);
  // draw_extrapolation(hMCReweight_bin_M, hSignal_M);
  // draw_extrapolation(hMCReweight_bin_M, hClosure_M);

  // draw_extrapolation(hReweight_fit_M, hSignal_M);
  // draw_extrapolation(hReweight_fit_M, hClosure_M);
  // draw_extrapolation(hMCReweight_fit_M, hSignal_M);
  // draw_extrapolation(hMCReweight_fit_M, hClosure_M);

  // draw_extrapolation(hReweight_expo_M, hSignal_M);
  // draw_extrapolation(hReweight_expo_M, hClosure_M);
  // draw_extrapolation(hMCReweight_expo_M, hSignal_M);
  // draw_extrapolation(hMCReweight_expo_M, hClosure_M);


  // const TString fout_name = prefix+"MC.Other.root";
  // TFile *fout = new TFile(fout_name, "update");
  // TDirectory* subdir = fout->mkdir(signalRegion+"_reco");
  // subdir->cd();
  // TH1F *hOut =(TH1F*) hReweight_fit_M->Clone("Bstar_reco_M");
  // hOut->Write();
  // fout->Close();
  // delete fout;

}

void fill_histogram_from_file(TH1F* hist, const TString &fName, const TString &region, const TString &histName, const int &weightSign) {

  TFile* f = new TFile(fName);
  hist->Add((TH1F*)  f->Get(region + histName), weightSign);
  f->Close();
  delete f;
}

TGraphAsymmErrors* fit_ratio(TH1F* hSignal, TH1F* hControl) {
  setCMSStyle();

  TGraphAsymmErrors* ratio = new TGraphAsymmErrors();
  TGraphErrors* ratio_int = new TGraphErrors();
  ratio->Divide(hSignal, hControl, "pois");

  
  // TF1 *ffit = new TF1("ffit", "TMath::Exp([0]*x/1000+[1])+[2]", 500, 5000);
  // TF1 *ffit = new TF1("ffit", "pol0", 500, 5000);
  // TF1 *ffit = new TF1("ffit", crystalball_function, 500, 5000,5);
  // ffit->SetParameters(0.04, 687.0, 276.0, -0.193, 0.17);
  
  // TF1 *ffit = new TF1("ffit", landau_flattail, 500, 5000, 6);
  // ffit->SetParameters(0.02, 1000.0, 400.0, 1500.0, 100.0, 0.002);
  // ffit->SetParLimits(3, 500.0, 5000.0);
  // ffit->SetParLimits(4, 1.0, 1000.0);

  TF1 *ffit = new TF1("ffit", "gaus(0)+pol0(3)", 500, 5000);
  ffit->SetParameters(0.01, 1450.0, 1500.0, 0.005);
  ffit->SetParLimits(0, 0.0, 0.1);
  ffit->SetParLimits(3, 0.0, 0.1);

  // TF1 *ffit = new TF1("ffit", "gaus", 500, 5000);
  
  ratio->Fit(ffit, "0", "", 500, 5000);
  int n_points = 10000;

  for (int i = 0; i<n_points; ++i)
    {
      ratio_int->SetPoint(i,500+3000-(3000*i/n_points),0);
    }
  ratio_int->SetFillColor(kGreen+1);
  (TVirtualFitter::GetFitter())->GetConfidenceIntervals(ratio_int,0.68);
  // Draw and safe
  TCanvas* c = new TCanvas("cfit","cfit",600,600);
  TPad* pad = SetupPad();

  gStyle->SetOptFit(1);
  pad->Draw();
  pad->cd();
  ratio_int->Draw("A3");
  ffit->SetLineColor(kAzure-3);
  ffit->SetLineWidth(2);
  ffit->Draw("same");
  ratio->Draw("P same");

  ratio_int->GetXaxis()->SetTitle("M_{tW} [GeV]");
  ratio_int->GetYaxis()->SetTitle("#alpha");
  ratio_int->GetYaxis()->CenterTitle();
  ratio_int->GetHistogram()->GetYaxis()->SetRangeUser(0, 0.3);
  TLegend* leg = new TLegend(0.25,0.66,0.5,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(ratio,"#alpha_{i}", "lpe");
  leg->AddEntry(ffit, "fit","l");
  leg->AddEntry(ratio_int,"68% CL.","f");
  
  leg->Draw();

  TString outputName; outputName.Form("ratiofit_%s.eps", hSignal->GetTitle());
  c->Print(outputName);

  return ratio;
}

void fill_reweighted(TH1F* hReweight, TF1* fitfun, const TString &fname, const int &weightSign) {
  TFile *f = new TFile(fname);
  TTreeReader tReader("AnalysisTree", f);
  TTreeReaderValue<float> preWeight(tReader, "final_weight");
  TTreeReaderValue<float> recoMass(tReader, "reco_mass");
  // TTreeReaderValue<double> st(tReader, "ST");
  while (tReader.Next()) 
    {
      float weight = fitfun->Eval(*recoMass);
      float mass = (*recoMass);
      weight = (weight < 0) ? 0 : weight * (*preWeight);
      hReweight->Fill(mass, weight * weightSign);
    }
  f->Close();
  delete f;
}

void fill_reweighted(TH1F* hReweight, TGraphAsymmErrors* ratio, const TString &fname, const int &weightSign) {
  TFile *f = new TFile(fname);
  TTreeReader tReader("AnalysisTree", f);
  TTreeReaderValue<float> preWeight(tReader, "final_weight");
  TTreeReaderValue<float> recoMass(tReader, "reco_mass");
  // TTreeReaderValue<double> st(tReader, "ST");
  while (tReader.Next()) {
    int bin = hReweight->FindBin(*recoMass);
    double weight = 0;
    double x = 0;
    float mass = (*recoMass) < 5000 ? *recoMass : 4999.9;
    ratio->GetPoint(bin-2, x, weight);
    weight = (weight < 0) ? 0 : weight * (*preWeight);
    // weight = weight*(*preWeight);
    hReweight->Fill(mass, weight * weightSign);
  }
  f->Close();
  delete f;
}

void draw_extrapolation(TH1F* hReweight, TH1F* hSignal) {
  setCMSStyle();
  gStyle->SetOptFit(0);
  TCanvas* c = new TCanvas("cmain","cmain",600,600);
  TPad* pad = SetupRatioPadTop();
  TPad* pad_ratio = SetupRatioPad();

  pad->Draw();
  pad_ratio->Draw();
  pad->cd();
  gPad->SetLogy();
  Hist_Cosmetics(hReweight);
  Hist_Cosmetics(hSignal);

  hReweight->SetFillColor(kGreen+1);
  hReweight->SetLineColor(kBlack);
  hReweight->SetLineWidth(2); 
  hReweight->SetMarkerStyle(1);
  hSignal->SetLineColor(kAzure);
  hSignal->SetLineWidth(2);
  hReweight->GetYaxis()->SetRangeUser(1.1e-4, hReweight->GetMaximum()*2.2);
  hReweight->Draw("H");
  hSignal->Draw("SAME P");

  TLegend* leg = new TLegend(0.5,0.7,0.8,0.85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hReweight,"extrapolation", "f");
  leg->AddEntry(hSignal, "MC prediction","lpe");
  leg->Draw();
  
  pad_ratio->cd();

  TGraphAsymmErrors *gRatio = new TGraphAsymmErrors(hSignal->GetNbinsX());
  gRatio->Divide(hSignal, hReweight, "pois");

  Hist_Cosmetics(gRatio, true); 
  gRatio->GetXaxis()->SetLimits(hReweight->GetXaxis()->GetXmin(), hReweight->GetXaxis()->GetXmax());
  double xmin = gRatio->GetXaxis()->GetXmin();
  double xmax = gRatio->GetXaxis()->GetXmax();
  TLine *line = new TLine(xmin, 1, xmax, 1);
  line->SetLineColor(kBlack);
  gRatio->GetYaxis()->CenterTitle();
  gRatio->GetYaxis()->SetRangeUser(0.3, 1.7);
  gRatio->GetYaxis()->SetTitle("VR/CR");
  gRatio->GetXaxis()->SetTitle("M_{tW} [GeV]");
  gRatio->GetYaxis()->SetNdivisions(505);
  gRatio->Draw("AP0");
  line->Draw();

  TString outputName; outputName.Form("closure_%s%s.eps", hSignal->GetTitle(), hReweight->GetTitle());
  c->Print(outputName);
}
