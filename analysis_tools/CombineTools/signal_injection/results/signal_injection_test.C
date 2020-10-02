
void draw_and_fit(TString file_name, Double_t signal_strength, TString plot_name);

void signal_injection_test(TString file_name, Double_t signal_strength) {

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  TString plot_name = file_name;
  plot_name.ReplaceAll(".root", "");
  draw_and_fit(file_name, signal_strength, plot_name);

}


void draw_and_fit(TString file_name, Double_t signal_strength, TString plot_name) {
  TFile* f = new TFile(file_name);
  TH1D* h = new TH1D("hist_r", "#hat{r} - r", 40, -4, 4);

  Double_t r, r_err;
  TTree* tree = (TTree*) f->Get("tree_fit_sb");
  tree->SetBranchAddress("r",&r);
  tree->SetBranchAddress("rErr",&r_err);

  for(int i=0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    h->Fill((r - signal_strength) / r_err);
  }
  
  TCanvas* c = new TCanvas("c","c", 600, 600);
  h->Fit("gaus");
  h->GetXaxis()->SetTitle(TString::Format("#hat{r} - %3.2f",signal_strength));
  
  h->Draw();

  c->Print(plot_name+".eps");
  f->Close();
}
