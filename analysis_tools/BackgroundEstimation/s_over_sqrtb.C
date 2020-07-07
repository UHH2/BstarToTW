#include "../cosmetics.h"
#include "../plot_styles.C"

void fill_histogram_from_file(TH1F* hist, const TString &fName, const TString &region, const int &weightSign = 1);

void s_over_sqrtb(TString channel = "Muon", const TString year = "all", const TString region = "1btag1toptag20chi2", const TString signal_name="BstarToTW1400LH") {

  const TString  in_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/"+channel+"/"+year+"/NOMINAL/";
  const TString prefix = "uhh2.AnalysisModuleRunner.MC.";
  const int n_backgrounds = 6;
  const TString background_names[n_backgrounds] = {"TT","ST","WJets", "DYJets", "Diboson", "QCD"};

  // const int n_signals = 20;
  // const TString signal_names[n_signals] = {"BstarToTW700LH", "BstarToTW800LH", "BstarToTW900LH", "BstarToTW1000LH", "BstarToTW1100LH", 
  // 					   "BstarToTW1200LH", "BstarToTW1400LH", "BstarToTW1600LH", "BstarToTW1800LH", "BstarToTW2000LH",
  // 					   "BstarToTW2200LH", "BstarToTW2400LH", "BstarToTW2600LH", "BstarToTW2800LH", "BstarToTW3000LH"};
  
  const int nbins = 23;
  const double xbins[nbins] = {0, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3000, 3500};

  setCMSStyle();

  TH1F *h_background = new TH1F("background", channel+"_"+region, nbins-1, xbins);
  TH1F* h_signal = new TH1F("signal", channel+"_"+region, nbins-1, xbins);
  TGraph* g_soverb = new TGraph();
  TCanvas* c = new TCanvas("c","c", 600, 600);
  TPad* pad = SetupPad();

  pad->Draw();
  pad->cd();


  for (int i = 0; i < n_backgrounds; ++i)
    {
      TString f_name = in_dir + prefix + background_names[i] + ".root";
      fill_histogram_from_file(h_background, f_name, region);
    }
  
  TString f_name = in_dir + prefix + signal_name + ".root";
  fill_histogram_from_file(h_signal, f_name, region);
  h_signal->Scale(0.8042);
  for (int point = 0; point < nbins; ++point)
    {
      double soverb = (h_background->GetBinContent(point+1) == 0) ? 0 : h_signal->GetBinContent(point+1) / std::sqrt(h_background->GetBinContent(point+1));
      g_soverb->SetPoint(point, xbins[point], soverb);
      std::cout << "N_bg= " << h_background->GetBinContent(point+1) << std::endl;
      std::cout << "N_sg= " << h_signal->GetBinContent(point+1) << std::endl;
      std::cout << "s/sqrt(b)= " << soverb << std::endl;
    } 

  // h_background->Draw();
  // h_signal->Draw("SAME");
  g_soverb->Draw("AP");  
  c->Print(region+"_"+signal_name+".eps");


}


void fill_histogram_from_file(TH1F* hist, const TString &fName, const TString &region, const int &weightSign) {

  TFile* f = new TFile(fName);
  hist->Add((TH1F*)  f->Get(region + "_reco/Bstar_reco_M_rebin"), weightSign);
  f->Close();
  delete f;
}

