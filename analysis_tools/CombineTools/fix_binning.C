

void fix_binning() {

  // open file
  TFile *fout = new TFile("BstarToTW_combined_new.root", "recreate");
  TFile *f = new TFile("BstarToTW_combined.root");
  double xbins[22] = {500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 2800, 3000, 3500};
  
  // loop over all histograms in file
  TIter nextKey(gDirectory->GetListOfKeys());
  TKey* key;
  TObject* obj; 
  while ((key = (TKey*)nextKey()) )
    {
      obj = key->ReadObj();
      if (obj->IsA()->InheritsFrom("TH1F")) 
	{
	  TH1F *h = (TH1F*) obj;
	  std::cout << key->GetName() << std::endl;
	  std::cout << h->GetName() << std::endl;
	  // remove empty bin
	  TH1F *h_rebin = (TH1F*) h->Clone();
	  h_rebin->SetBins(21,xbins);
	  fout->cd();
	  h_rebin->Write(key->GetName());
	  f->cd();
	  delete h_rebin;
	}	
    }


}
