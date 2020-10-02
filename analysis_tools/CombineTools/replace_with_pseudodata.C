

void replace_with_pseudodata(TString filename) {  
  // setup
  TString backgrounds[] = {"TT", "ST", "Other"};
  // TString regions[] = {"muon_region_1b1t", "muon_region_2b1t", "electron_region_1b1t", "electron_region_2b1t"};
  TString regions[] = {"combined_region_1b1t", "combined_region_2b1t"};

  
  TFile *f = new TFile(filename,"UPDATE");
  for (TString region : regions)
    {
      TH1F* h = (TH1F*) f->Get(region+"__data_obs")->Clone();
      h->Reset();
      for (TString background : backgrounds) 
	{
	  h->Add((TH1F*) f->Get(region+"__"+background));
	}
      h->Write(region+"__data_obs",TObject::kOverwrite);
    }
  f->Close();
}
