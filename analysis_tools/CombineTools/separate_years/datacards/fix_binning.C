

void fix_binning() {

  // open file
  TString filenames[] = {"Electron_2016","Electron_2017","Electron_2018","Muon_2016","Muon_2017","Muon_2018"};
  const int nbins = 20;
  double xbins[nbins] = {500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2200, 2400, 2600, 4000}; // new extended binning

  for (TString filename : filenames)
    {

      TFile *fout = new TFile("BstarToTW_"+filename+"_new.root", "recreate");
      TFile *f = new TFile("BstarToTW_"+filename+".root");
  
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
	      // std::cout << key->GetName() << std::endl;
	      // std::cout << h->GetName() << std::endl;
	      // remove empty bin
	      TH1F *h_rebin = (TH1F*) h->Rebin(nbins-1,"h_new",xbins);;
	      fout->cd();
	      h_rebin->Write(key->GetName());
	      f->cd();
	      delete h_rebin;
	    }	
	}

      fout->Close();
      f->Close();
    }


}
