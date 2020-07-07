#include "../cosmetics.h"

void draw_stat_uncert() {

  TFile *f = new TFile("BstarToTW_combined.root");

  TString regions[] = {"muon_region_1b1t__", "muon_region_2b1t__", "electron_region_1b1t__", "electron_region_2b1t__"};

  for (TString region : regions)
    {
      TH1F *h = (TH1F*) f->Get(region+"TT");
      h->Add((TH1F*) f->Get(region+"ST"));
      h->Add((TH1F*) f->Get(region+"Other"));
      
      cout << "statistical uncertainty in " << region << endl;

      for (int i = 1; i <= h->GetNbinsX(); ++i)
	{
	  cout << "bin " << i << ": " << h->GetBinError(i) / h->GetBinContent(i) << endl;
	}
    }

  // TCanvas *c = new TCanvas("c", "c", 1200, 1200);
  // TPad* pad = SetupPad();
  // pad->Draw();
  // pad->cd();


}
