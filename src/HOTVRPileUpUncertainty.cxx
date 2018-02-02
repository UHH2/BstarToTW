#include "UHH2/BstarToTW/include/HOTVRPileUpUncertainty.h"

using namespace uhh2;
using namespace std;

HOTVRPileUpUncertainty::HOTVRPileUpUncertainty(Context &ctx, TString path, TString sys_direction) {

  is_mc = ctx.get("dataset_type") == "MC";
  unique_ptr<TFile> file;
  file.reset(new TFile(path, "READ"));

  sf_hist.reset((TGraph*)file->Get("RhoRes"));

 if(sys_direction != "central" && sys_direction != "up" && sys_direction != "down") throw runtime_error("HOTVRPileUpUncertainty: Invalid sys_direction specified.");
  m_sys_direction = sys_direction;

}

bool HOTVRPileUpUncertainty::process(Event &event) {
  if (is_mc)
    {
      
      double rho = event.rho;
      double x,y;
      if (rho < 5) sf_hist->GetPoint(0, x, y);
      else if (rho < 10) sf_hist->GetPoint(1, x, y);
      else if (rho < 15) sf_hist->GetPoint(2, x, y);
      else if (rho < 20) sf_hist->GetPoint(3, x, y);
      else if (rho < 25) sf_hist->GetPoint(4, x, y);
      else if (rho < 30) sf_hist->GetPoint(5, x, y);
      else if (rho < 35) sf_hist->GetPoint(6, x, y);
      else if (rho < 40) sf_hist->GetPoint(7, x, y);
      else if (rho < 45) sf_hist->GetPoint(8, x, y);
      else if (rho < 50) sf_hist->GetPoint(9, x, y);
      else if (rho < 55) sf_hist->GetPoint(10, x, y);
      else if (rho > 55) sf_hist->GetPoint(11, x, y);

      double variation = abs(1 - y);

      if (m_sys_direction == "central")
	return true;
      else if (m_sys_direction == "up")
	{
	  for (auto & topjet : *event.topjets)
	    {
	      vector<Jet> new_subjets;
	      LorentzVector temp_jet;
	      vector<Jet> subjets = topjet.subjets();
	      for (auto & subjet : subjets)
		{
		  subjet.set_v4(subjet.v4() * (1+variation));
		  temp_jet += subjet.v4();
		  new_subjets.push_back(subjet);
		}
	      topjet.set_subjets(move(new_subjets));
	      topjet.set_v4(temp_jet);
	    }
	}
      else if (m_sys_direction == "down")
	{
	  for (auto & topjet : *event.topjets)
	    {
	      vector<Jet> new_subjets;
	      LorentzVector temp_jet;
	      vector<Jet> subjets = topjet.subjets();
	      for (auto & subjet : subjets)
		{
		  subjet.set_v4(subjet.v4() * (1 - variation));
		  temp_jet += subjet.v4();
		  new_subjets.push_back(subjet);
		}
	      topjet.set_subjets(move(new_subjets));
	      topjet.set_v4(temp_jet);
	    }
	}	
    }
  return true;
}
