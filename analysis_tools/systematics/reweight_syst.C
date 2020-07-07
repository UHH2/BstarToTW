#include "../cosmetics.h"
#include "../plot_styles.C"

TH1F* get_hist_from_file(const TString file_name, const TString hist_name);

TH1F* get_reweighted_hist(TH1F* h_num, TH1F* h_den, const TString name);
void write_hist_to_file(const TString fname, const TString subdir_name, TH1F* hist);

void reweight_syst(){
  TString year = "2016";
  // TString year = "2017";
  // TString year = "2018";
  // TString year = "all";
  
  TString signal_postfix = "_2016v3";
  // TString signal_postfix = "_2017v2";
  // TString signal_postfix = "_2018";
  // TString signal_postfix = "";

  TString in_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Electron/"+year+"/";
  // TString in_dir = "/nfs/dust/cms/user/froehlia/BstarToTW/Run2_102X_v1/Analysis/Muon/"+year+"/";
  
  TString processes[] = {"TT",  "ST",
			 "BstarToTW0700RH"+signal_postfix, "BstarToTW0800RH"+signal_postfix, "BstarToTW0900RH"+signal_postfix, "BstarToTW1000RH"+signal_postfix, "BstarToTW1100RH"+signal_postfix,
			 "BstarToTW1200RH"+signal_postfix, "BstarToTW1400RH"+signal_postfix, "BstarToTW1600RH"+signal_postfix, "BstarToTW1800RH"+signal_postfix,
			 "BstarToTW2000RH"+signal_postfix, "BstarToTW2200RH"+signal_postfix, "BstarToTW2400RH"+signal_postfix, "BstarToTW2600RH"+signal_postfix, 
			 "BstarToTW2800RH"+signal_postfix, "BstarToTW3000RH"+signal_postfix, "BstarToTW3200RH"+signal_postfix, "BstarToTW3400RH"+signal_postfix,
			 "BstarToTW3600RH"+signal_postfix, "BstarToTW3800RH"+signal_postfix, "BstarToTW4000RH"+signal_postfix, "BstarToTW4200RH"+signal_postfix,
			 "BstarToTW0700LH"+signal_postfix, "BstarToTW0800LH"+signal_postfix, "BstarToTW0900LH"+signal_postfix, "BstarToTW1000LH"+signal_postfix, "BstarToTW1100LH"+signal_postfix,
			 "BstarToTW1200LH"+signal_postfix, "BstarToTW1400LH"+signal_postfix, "BstarToTW1600LH"+signal_postfix, "BstarToTW1800LH"+signal_postfix,
			 "BstarToTW2000LH"+signal_postfix, "BstarToTW2200LH"+signal_postfix, "BstarToTW2400LH"+signal_postfix, "BstarToTW2600LH"+signal_postfix, 
			 "BstarToTW2800LH"+signal_postfix, "BstarToTW3000LH"+signal_postfix, "BstarToTW3200LH"+signal_postfix, "BstarToTW3400LH"+signal_postfix,
			 "BstarToTW3600LH"+signal_postfix, "BstarToTW3800LH"+signal_postfix, "BstarToTW4000LH"+signal_postfix, "BstarToTW4200LH"+signal_postfix,
			 "BstarToTW0700VL"+signal_postfix, "BstarToTW0800VL"+signal_postfix, "BstarToTW0900VL"+signal_postfix, "BstarToTW1000VL"+signal_postfix, "BstarToTW1100VL"+signal_postfix,
			 "BstarToTW1200VL"+signal_postfix, "BstarToTW1400VL"+signal_postfix, "BstarToTW1600VL"+signal_postfix, "BstarToTW1800VL"+signal_postfix,
 			 "BstarToTW2000VL"+signal_postfix, "BstarToTW2200VL"+signal_postfix, "BstarToTW2400VL"+signal_postfix, "BstarToTW2600VL"+signal_postfix, 
			 "BstarToTW2800VL"+signal_postfix, "BstarToTW3000VL"+signal_postfix, "BstarToTW3200VL"+signal_postfix, "BstarToTW3400VL"+signal_postfix,
			 "BstarToTW3600VL"+signal_postfix, "BstarToTW3800VL"+signal_postfix, "BstarToTW4000VL"+signal_postfix, "BstarToTW4200VL"+signal_postfix,
			 "BprimeTToTW0700RH_2016v3", "BprimeTToTW0800RH_2016v3", "BprimeTToTW0900RH_2016v3", "BprimeTToTW1000RH_2016v3", "BprimeTToTW1100RH_2016v3", "BprimeTToTW1200RH_2016v3",
			 "BprimeTToTW1300RH_2016v3", "BprimeTToTW1400RH_2016v3", "BprimeTToTW1500RH_2016v3", "BprimeTToTW1600RH_2016v3", "BprimeTToTW1700RH_2016v3", "BprimeTToTW1800RH_2016v3",
			 "BprimeTToTW0700LH_2016v3", "BprimeTToTW0800LH_2016v3", "BprimeTToTW0900LH_2016v3", "BprimeTToTW1000LH_2016v3", "BprimeTToTW1100LH_2016v3", "BprimeTToTW1200LH_2016v3",
			 "BprimeTToTW1300LH_2016v3", "BprimeTToTW1400LH_2016v3", "BprimeTToTW1500LH_2016v3", "BprimeTToTW1600LH_2016v3", "BprimeTToTW1700LH_2016v3", "BprimeTToTW1800LH_2016v3"};

  const TString prefix = "uhh2.AnalysisModuleRunner.MC.";
  
  const TString num_name = "JER";
  const TString den_name = "NOMINAL";

  const TString channels[] = {"1btag1toptag20chi2_reco", "2btag1toptag_reco"};
  const TString hist_name = "Bstar_reco_M_rebin";

  for (TString process : processes)
    {
      for (TString channel: channels)
	{
	  // get histogram of numerator
	  TH1F *h_num_up = get_hist_from_file(in_dir+num_name+"_up/"+prefix+process+".root", channel +"/"+ hist_name);
	  TH1F *h_num_down = get_hist_from_file(in_dir+num_name+"_down/"+prefix+process+".root", channel +"/"+ hist_name);
	  TH1F *h_den = get_hist_from_file(in_dir+den_name+"/"+prefix+process+".root", channel +"/"+ hist_name);
	  
	  TH1F *h_reweight_up = get_reweighted_hist(h_num_up, h_den, process+"_"+channel+"_up");
	  TH1F *h_reweight_down = get_reweighted_hist(h_num_down, h_den, process+"_"+channel+"_down");
	  
	  write_hist_to_file("up/" + prefix + process + ".root", channel, h_reweight_up);
	  write_hist_to_file("down/" + prefix + process + ".root", channel, h_reweight_down);
	}
    }
}


TH1F* get_hist_from_file(const TString file_name, const TString hist_name) {
  std::cout << "get histogram " << hist_name << "from file " << file_name << std::endl;

  TFile *f = new TFile(file_name);
  TH1F *hist = (TH1F*)( f->Get(hist_name))->Clone();
  hist->SetDirectory(0);
  f->Close();
  delete f;
  return hist;
}

TH1F* get_reweighted_hist(TH1F* h_num, TH1F* h_den, const TString name) {
  std::cout << "reweight hist ..." << std::endl;
  
  TH1F *h_reweight = (TH1F*) h_den->Clone();
  TGraphAsymmErrors* ratio = new TGraphAsymmErrors();
  ratio->Divide(h_num, h_den, "pois");
  TF1* ffit = new TF1("ffit", "pol1", 0., 4000.);
  ratio->Fit(ffit);
  
  for (int i = 0; i <= h_den->GetNbinsX(); ++i)
    {
      h_reweight->SetBinContent(i, h_den->GetBinContent(i) * ffit->Eval(h_den->GetBinCenter(i)));
    } 
  setCMSStyle();
  TCanvas* c = new TCanvas("c","c",600,600);
  TPad* pad = SetupPad();
  pad->Draw();
  pad->cd();
  ratio->GetHistogram()->GetYaxis()->SetRangeUser(0.,2.);
  ratio->Draw("AP");
  ffit->Draw("SAME");

  c->Print("plots/"+name+".eps");
  
  delete c;
  delete ffit;
  delete ratio;
  return  h_reweight;
}

void write_hist_to_file(const TString fname, const TString subdir_name, TH1F* hist) {
  std::cout << "write to file " << fname << std::endl;

  TFile *fout = new TFile(fname, "update");
  TDirectory* subdir = fout->mkdir(subdir_name);
  subdir->cd();
  TH1F *h_out = (TH1F*) hist->Clone("Bstar_reco_M_rebin");
  h_out->Write();
  fout->Close();
  delete fout;
}
