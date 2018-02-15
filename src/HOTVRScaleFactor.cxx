#include "UHH2/BstarToTW/include/HOTVRScaleFactor.h"


using namespace uhh2;
using namespace std;

HOTVRScaleFactor::HOTVRScaleFactor(Context &ctx, string signal_name, TopJetId id_topjet, TString path, TString sys_direction) {
  string dataset_name = ctx.get("dataset_version");
  m_do_weight = (dataset_name.find("TTbar") != std::string::npos || dataset_name.find(signal_name) != std::string::npos);
  m_id_topjet = id_topjet;
  
  unique_ptr<TFile> file;
  file.reset(new TFile(path, "READ"));

  sf_hist.reset((TGraphAsymmErrors*)file->Get("ratio"));

  if(sys_direction != "central" && sys_direction != "up" && sys_direction != "down") throw runtime_error("HOTVRScaleFactor: Invalid sys_direction specified.");
  m_sys_direction = sys_direction;
  
}

bool HOTVRScaleFactor::process(Event &event) {
  if (m_do_weight)
    {
      for (TopJet topjet : *event.topjets)
	{
	  ///////////////////////////////////////////////////////////////////////
	  // double eta = topjet.eta();
	  // if (abs(eta) > 2.5) throw runtime_error("HOTVRScaleFactor: TopJet |eta| > 2.5 is not supported");
	  // if (m_id_topjet(topjet, event))
	  //   {
	  //     if (m_sys_direction == "up")
	  // 	{
	  // 	  if (eta < -1.479)      event.weight *= sf_hist->GetErrorYhigh(0);
	  // 	  else if (eta < 0)      event.weight *= sf_hist->GetErrorYhigh(1);
	  // 	  else if (eta < 1.479)  event.weight *= sf_hist->GetErrorYhigh(2);
	  // 	  else                   event.weight *= sf_hist->GetErrorYhigh(3);
	  // 	}
	  //     else if (m_sys_direction == "down")
	  // 	{
	  // 	  if (eta < -1.479)      event.weight *= sf_hist->GetErrorYlow(0);
	  // 	  else if (eta < 0)      event.weight *= sf_hist->GetErrorYlow(1);
	  // 	  else if (eta < 1.479)  event.weight *= sf_hist->GetErrorYlow(2);
	  // 	  else                   event.weight *= sf_hist->GetErrorYlow(3);
	  // 	}
	  //     else if (m_sys_direction == "central")
	  // 	{
	  // 	  double x,y;
	  // 	  if (eta < -1.479)
	  // 	    {
	  // 	      sf_hist->GetPoint(0, x, y); 
	  // 	      event.weight *= y;
	  // 	    }
	  // 	  else if (eta < 0)      
	  // 	    {
	  // 	      sf_hist->GetPoint(1, x, y); 
	  // 	      event.weight *= y;
	  // 	    }
	  // 	  else if (eta < 1.479)  
	  // 	    {
	  // 	      sf_hist->GetPoint(2, x, y); 
	  // 	      event.weight *= y;
	  // 	    }
	  // 	  else     
	  // 	    {
	  // 	      sf_hist->GetPoint(3, x, y); 
	  // 	      event.weight *= y;
	  // 	    }              
	  // 	}
	  //   }
	  //////////////////////////////////////////////////////////////////////////

	  double pt = topjet.pt();
	  if (m_id_topjet(topjet, event))
	    {
	      if (m_sys_direction == "central")
	  	{
	  	  double x,y;
		  if (pt < 300) 
	  	    {
	  	      sf_hist->GetPoint(0, x, y); 
	  	      event.weight *= y;
	  	    }
    	  	  else if (pt < 400)
	  	    {
	  	      sf_hist->GetPoint(1, x, y); 
	  	      event.weight *= y;
	  	    }
	  	  else
	  	    {
	  	      sf_hist->GetPoint(2, x, y); 
	  	      event.weight *= y;
	  	    }
		}
	      else if (m_sys_direction == "up")
		{
		  double x,y;
		  if (pt < 300) 
		    {
		      sf_hist->GetPoint(0, x, y); 
		      event.weight *= y + sf_hist->GetErrorYhigh(0);
		    }
		  else if (pt < 400)
		    {
		      sf_hist->GetPoint(1, x, y); 
		      event.weight *= y + sf_hist->GetErrorYhigh(1);
		    }
		  else
		    {
		      sf_hist->GetPoint(2, x, y); 
		      event.weight *= y +  sf_hist->GetErrorYhigh(2);
		    }
		}
	      else if (m_sys_direction == "down")
		{
		  double x,y;
		  if (pt < 300) 
		    {
		      sf_hist->GetPoint(0, x, y); 
		      event.weight *= y + sf_hist->GetErrorYlow(0);
		    }
		  else if (pt < 400)
		    {
		      sf_hist->GetPoint(1, x, y); 
		      event.weight *= y + sf_hist->GetErrorYlow(1);
		    }
		  else
		    {
		      sf_hist->GetPoint(2, x, y); 
		      event.weight *= y +  sf_hist->GetErrorYlow(2);
		    }
		}
	    }
	} 
    }
  return true;   
}

CMSTTScaleFactor::CMSTTScaleFactor(Context &ctx, string signal_name, TString sys_direction){

  string dataset_name = ctx.get("dataset_version");
  m_do_weight = (dataset_name.find("TTbar") != std::string::npos || dataset_name.find(signal_name) != std::string::npos);

  if(sys_direction != "central" && sys_direction != "up" && sys_direction != "down") throw runtime_error("HOTVRScaleFactor: Invalid sys_direction specified.");
  m_sys_direction = sys_direction;

}

bool CMSTTScaleFactor::process(Event &event) {
  if (m_do_weight)
    {
      if (m_sys_direction == "central")
	event.weight *= 1.07;
      else if (m_sys_direction == "up")
	event.weight *= 1.07+0.05;
      else if (m_sys_direction == "down")
	event.weight *= 1.07-0.03;
    }
  return true;
}
