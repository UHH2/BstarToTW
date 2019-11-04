#include "UHH2/BstarToTW/include/BstarToTWModules.h"

using namespace std;
using namespace uhh2;

ObjectTagger::ObjectTagger() {}

void ObjectTagger::init(Context &ctx) {
  fail_if_init();

  h_toptag = ctx.declare_event_output<vector<TopJet> >(m_toptag_name);
  h_btag_loose = ctx.declare_event_output<vector<Jet> >(m_btag_loose_name);
  h_btag_medium = ctx.declare_event_output<vector<Jet> >(m_btag_medium_name);
  h_btag_tight = ctx.declare_event_output<vector<Jet> >(m_btag_tight_name);

  b_init = true;
}

void ObjectTagger::fail_if_init() {
  if(b_init)
    {
      throw logic_error("ObjectTagger::init already called!");
    }
}

bool ObjectTagger::process(Event &event) {
  if(!b_init){
    throw runtime_error("ObjectTagger::init not called (has to be called in AnalysisModule constructor)!");
  }
  // tag and store
  vector<Jet> btag_loose, btag_medium, btag_tight;
  vector<TopJet> toptag;
  for (const Jet jet: *event.jets)
    {
      if (m_btag_tight_id(jet, event))
	{
	  btag_tight.push_back(jet);
	  btag_medium.push_back(jet);
	  btag_loose.push_back(jet);
	}
      else if (m_btag_medium_id(jet, event))
	{
	  btag_medium.push_back(jet);
	  btag_loose.push_back(jet);
	}
      else if (m_btag_loose_id(jet, event))
	{
	  btag_loose.push_back(jet);
	}
    }
  event.set(h_btag_tight, btag_tight);
  event.set(h_btag_medium, btag_medium);
  event.set(h_btag_loose, btag_loose);

  for (const TopJet topjet : *event.topjets)
    {
      if (m_toptag_id(topjet, event))
	toptag.push_back(topjet);
    }
  event.set(h_toptag, toptag);
  return true;
}

SystematicsModule::SystematicsModule(){}

bool SystematicsModule::process(Event &event){
  for (const auto &module : modules)
    {
      module->process(event);
    }
  return true;
}

ElectronTriggerWeights::ElectronTriggerWeights(Context & ctx, TString path_, TString SysDirection_): path(path_), SysDirection(SysDirection_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: ElectronTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }
  
  unique_ptr<TFile> file;												 
  file.reset(new TFile(path,"READ"));											 
  
  Eff_lowpt_MC.reset((TGraphAsymmErrors*)file->Get("gr_lowpt_eta_TTbar_eff"));					 
  Eff_highpt_MC.reset((TGraphAsymmErrors*)file->Get("gr_highpt_eta_TTbar_eff"));					 
  Eff_lowpt_DATA.reset((TGraphAsymmErrors*)file->Get("gr_lowpt_eta_DATA_eff"));					 
  Eff_highpt_DATA.reset((TGraphAsymmErrors*)file->Get("gr_highpt_eta_DATA_eff"));					 
  
  if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In BstarToTWModules.cxx, ElectronTriggerWeights.process(): Invalid SysDirection specified.");
  
}

bool ElectronTriggerWeights::process(Event & event){

  if(event.isRealData) return true;

  const auto ele = event.electrons->at(0);
  double eta = ele.eta();
  if(fabs(eta) > 2.4) throw runtime_error("In BstarToTWModules.cxx, ElectronTriggerWeights.process(): Ele-|eta| > 2.4 is not supported at the moment.");


  //find right bin in eta
  int idx = 0;
  bool lowpt = false;
  if(30 <= ele.pt() && ele.pt() < 120){
    lowpt = true;
    //lowpt trigger
    bool keep_going = true;
    while(keep_going){
      double x,y;
      Eff_lowpt_MC->GetPoint(idx,x,y);
      keep_going = eta > x + Eff_lowpt_MC->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }
  }
  else if(ele.pt() >= 120){
    //highpt trigger
    bool keep_going = true;
    while(keep_going){
      double x,y;
      Eff_highpt_MC->GetPoint(idx,x,y);
      keep_going = eta > x + Eff_highpt_MC->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }
  }
  else throw runtime_error("In BstarToTWModules.cxx, ElectronTriggerWeights.process(): Electron has pt<30. Clean electron collection before applying weights.");

  //access efficiencies for MC and DATA, possibly accout for systematics = statistical + add. 2% up/down
  double eff_data = -1, eff_mc = -1, dummy_x;
  double stat_data = -1, stat_mc = -1, tp = 0.02, total_syst_data = -1, total_syst_mc = -1;
  if(lowpt){
    Eff_lowpt_MC->GetPoint(idx,dummy_x,eff_mc);
    Eff_lowpt_DATA->GetPoint(idx,dummy_x,eff_data);

    if(SysDirection == "up"){		
      stat_mc = Eff_lowpt_MC->GetErrorYlow(idx);	
      stat_data = Eff_lowpt_DATA->GetErrorYhigh(idx);
      total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
      total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));

      eff_mc -= total_syst_mc;    	
      eff_data += total_syst_data;	
    }							
    else if(SysDirection == "down"){
      stat_mc = Eff_lowpt_MC->GetErrorYhigh(idx);	
      stat_data = Eff_lowpt_DATA->GetErrorYlow(idx);
      total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
      total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
			
      eff_mc += Eff_lowpt_MC->GetErrorYhigh(idx);    	
      eff_data -= Eff_lowpt_DATA->GetErrorYlow(idx);	
    }                                                 
  }
  else{
    Eff_highpt_MC->GetPoint(idx,dummy_x,eff_mc);
    Eff_highpt_DATA->GetPoint(idx,dummy_x,eff_data);

    if(SysDirection == "up"){	
      stat_mc = Eff_highpt_MC->GetErrorYlow(idx);	
      stat_data = Eff_highpt_DATA->GetErrorYhigh(idx);
      total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
      total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
			
      eff_mc -= Eff_highpt_MC->GetErrorYlow(idx);    	
      eff_data += Eff_highpt_DATA->GetErrorYhigh(idx);	
    }							
    else if(SysDirection == "down"){	
      stat_mc = Eff_highpt_MC->GetErrorYhigh(idx);	
      stat_data = Eff_highpt_DATA->GetErrorYlow(idx);
      total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
      total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
					
      eff_mc += Eff_highpt_MC->GetErrorYhigh(idx);    	
      eff_data -= Eff_highpt_DATA->GetErrorYlow(idx);	
    }                                                 
  }

  //Scale weight by (eff_data) / (eff_mc)
  double SF = eff_data/eff_mc;
  event.weight *= SF;

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

BstarToTWOutputModule::BstarToTWOutputModule(uhh2::Context &ctx, const std::string hyps_name){
  h_weight   = ctx.declare_event_output<float>("final_weight");
  h_topjets = ctx.declare_event_output<vector<vector<double> > >("topjets");
  h_recomass = ctx.declare_event_output<float>("reco_mass");
  h_hyps     = ctx.get_handle<std::vector<BstarToTWHypothesis> >(hyps_name);
}

bool BstarToTWOutputModule::process(Event &event) {
  // get mbstar
  std::vector<BstarToTWHypothesis> hyps = event.get(h_hyps);
  const BstarToTWHypothesis* hyp = get_best_hypothesis(hyps, "Chi2");
  double mbstar = 0;
  if((hyp->get_topjet() + hyp->get_w()).isTimelike())
    {    
      mbstar = (hyp->get_topjet() + hyp->get_w()).M();
    }
  else
    {
      mbstar = sqrt(-(hyp->get_topjet()+hyp->get_w()).mass2());
    }

  // get topjet fourvectors
  vector<TopJet> topjets = *event.topjets;
  vector<vector<double> > topjets_v4;
  for (TopJet &topjet : topjets)
    {
      vector<double> v4;
      v4.push_back(topjet.pt());
      v4.push_back(topjet.eta());
      v4.push_back(topjet.phi());
      v4.push_back(topjet.energy());
      topjets_v4.push_back(v4);      
    }
  
  event.set(h_topjets, topjets_v4);
  event.set(h_weight, event.weight);
  event.set(h_recomass, mbstar);
  return true;
}

