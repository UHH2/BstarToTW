
void goodnessoffit_analyzer() {

  TFile* f_test_stat = new TFile("higgsCombine_result_bonly_CRonly_toy.GoodnessOfFit.mH120.123456.root");

  TH1D* h_test_stat = new TH1D("Goodness of Fit", "test_stat", 90, 0, 90);

  double limit;
  TTree* tree = (TTree*) f_test_stat->Get("limit");
  tree->SetBranchAddress("limit",&limit);

  for(int i=0; i < tree->GetEntries(); ++i){

    tree->GetEntry(i);
    h_test_stat->Fill(limit);
  }

  h_test_stat->Scale(1/h_test_stat->Integral());

  TFile* f_data = new TFile("higgsCombine_result_bonly_CRonly.GoodnessOfFit.mH120.root");
  TTree* data_tree = (TTree*) f_data->Get("limit");
  double data_limit;
  data_tree->SetBranchAddress("limit",&data_limit);
  data_tree->GetEntry(0);
  
  TCanvas* c = new TCanvas("c","c",600,600);
  
  h_test_stat->Fit("gaus");
  TF1* fgaus = h_test_stat->GetFunction("gaus");

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_test_stat->Draw();



  TArrow* ar = new TArrow(data_limit,fgaus->GetParameter(0)*0.33,data_limit,0.0, 0.05, "|>");
  ar->SetLineWidth(2);
  ar->Draw();

  double p_value = fgaus->Integral(data_limit, TMath::Infinity());
  cout << "p-value of data: " << p_value << endl;

  TLatex *text = new TLatex();
  text->DrawTextNDC(0.15, 0.8, TString::Format("p-value: %f", p_value));
  
  c->Print("GoF_test.eps");
}
