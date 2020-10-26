#include "UHH2/BstarToTW/include/BstarToTWModules.h"

#include "boost/algorithm/string.hpp"

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

ElectronTriggerWeights::ElectronTriggerWeights(Context & ctx, TString path_, TString syst_direction_): path(path_), syst_direction(syst_direction_) {

  h_ele_weight      = ctx.declare_event_output<float>("weight_sfelec_trigger");
  h_ele_weight_up   = ctx.declare_event_output<float>("weight_sfelec_trigger_up");
  h_ele_weight_down = ctx.declare_event_output<float>("weight_sfelec_trigger_down");

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: ElectronTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }
  year = extract_year(ctx);
  TString yeartag = "2016";
  if(year == Year::is2017v1 || year == Year::is2017v2) yeartag = "2017";
  else if(year == Year::is2018) yeartag = "2018";
  unique_ptr<TFile> file_pt1, file_pt2;
  if(yeartag == "2016"){
    file_pt1.reset(new TFile(path+"/" + yeartag + "/ElectronTriggerScaleFactors_eta_ele_binned_official_pt30to175_withsyst.root","READ"));
    file_pt2.reset(new TFile(path+"/" + yeartag + "/ElectronTriggerScaleFactors_eta_ele_binned_official_pt175toInf.root","READ"));
  }
  else if(yeartag == "2017" || yeartag == "2018"){
    file_pt1.reset(new TFile(path+"/" + yeartag + "/ElectronTriggerScaleFactors_eta_ele_binned_official_pt30to200_withsyst.root","READ"));
    file_pt2.reset(new TFile(path+"/" + yeartag + "/ElectronTriggerScaleFactors_eta_ele_binned_official_pt200toInf.root","READ"));
  }
  else throw runtime_error("invalid year");

  if(yeartag == "2016") pt_bins = {30, 175};
  else if(yeartag == "2017" || yeartag == "2018") pt_bins = {30, 200};

  g_sf_pt1.reset((TGraphAsymmErrors*)file_pt1->Get("ScaleFactors"));
  g_sf_pt2.reset((TGraphAsymmErrors*)file_pt2->Get("ScaleFactors"));
}

bool ElectronTriggerWeights::process(Event & event){

  if(event.isRealData){
    event.set(h_ele_weight, 1.);
    event.set(h_ele_weight_up, 1.);
    event.set(h_ele_weight_down, 1.);
    return true;
  }

  const Electron ele = event.electrons->at(0);
  double eta = ele.eta();
  double pt = ele.pt();
  if(fabs(eta) > 2.4) throw runtime_error("In BstarToTWModules.cxx, ElectronTriggerWeights.process(): |eta| > 2.4 is not supported.");
  if(pt < pt_bins[0]) throw runtime_error("In BstarToTWModules.cxx, ElectronTriggerWeights.process(): too small ele-pt for this year, not supported.");

  // find number of correct eta bin
  int idx = 0;
  if(pt < pt_bins[1]){
    bool keep_going = true;
    while(keep_going){
      double x,y;
      g_sf_pt1->GetPoint(idx,x,y);
      keep_going = eta > x + g_sf_pt1->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }
  }
  else {
    bool keep_going = true;
    while(keep_going){
      double x,y;
      g_sf_pt2->GetPoint(idx,x,y);
      keep_going = eta > x + g_sf_pt2->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }
  } 
  //access scale factors, add 2% t&p systematic uncertainty
  double sf, sf_up, sf_down, dummy_x;
  double stat_up = -1., stat_down = -1., tp = 0.0, total_up = -1., total_down = -1.;
  if(pt < pt_bins[1]){
    g_sf_pt1->GetPoint(idx,dummy_x,sf);

    stat_up = g_sf_pt1->GetErrorYhigh(idx);
    stat_down = g_sf_pt1->GetErrorYlow(idx);
    total_up = sqrt(pow(stat_up,2) + pow(tp,2));
    total_down = sqrt(pow(stat_down,2) + pow(tp,2));

    sf_up = sf + total_up;
    sf_down = sf - total_down;
  }
  else {
    g_sf_pt2->GetPoint(idx,dummy_x,sf);

    stat_up = g_sf_pt2->GetErrorYhigh(idx);
    stat_down = g_sf_pt2->GetErrorYlow(idx);
    total_up = sqrt(pow(stat_up,2) + pow(tp,2));
    total_down = sqrt(pow(stat_down,2) + pow(tp,2));

    sf_up = sf + total_up;
    sf_down = sf - total_down;
  }

  event.set(h_ele_weight, sf);
  event.set(h_ele_weight_up, sf_up);
  event.set(h_ele_weight_down, sf_down);

  if (syst_direction == "up") {
    event.weight *= sf_up;
  } else if (syst_direction == "down") {
    event.weight *= sf_down;
  } else {
    event.weight *= sf;
  }

  return true;
}


MuonTriggerWeights::MuonTriggerWeights(Context & ctx, TString path_, TString syst_direction_): path(path_), syst_direction(syst_direction_) {

  h_muo_weight      = ctx.declare_event_output<float>("weight_sfmuon_trigger");
  h_muo_weight_up   = ctx.declare_event_output<float>("weight_sfmuon_trigger_up");
  h_muo_weight_down = ctx.declare_event_output<float>("weight_sfmuon_trigger_down");

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: MuonTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }
  year = extract_year(ctx);
  TString yeartag = "2016";
  if(year == Year::is2017v1 || year == Year::is2017v2) yeartag = "2017";
  else if(year == Year::is2018) yeartag = "2018";

  unique_ptr<TFile> file_0to0p9, file_0p9to1p2, file_1p2to2p1, file_2p1to2p4;
  file_0to0p9.reset(new TFile(path+"/" + yeartag + "/MuonTriggerScaleFactors_pt_mu_binned_official_eta0to0p9_withsyst.root","READ"));
  file_0p9to1p2.reset(new TFile(path+"/" + yeartag + "/MuonTriggerScaleFactors_pt_mu_binned_official_eta0p9to1p2_withsyst.root","READ"));
  file_1p2to2p1.reset(new TFile(path+"/" + yeartag + "/MuonTriggerScaleFactors_pt_mu_binned_official_eta1p2to2p1_withsyst.root","READ"));
  file_2p1to2p4.reset(new TFile(path+"/" + yeartag + "/MuonTriggerScaleFactors_pt_mu_binned_official_rebin_eta2p1to2p4_withsyst.root","READ"));


  g_sf_0to0p9.reset((TGraphAsymmErrors*)file_0to0p9->Get("ScaleFactors"));
  g_sf_0p9to1p2.reset((TGraphAsymmErrors*)file_0p9to1p2->Get("ScaleFactors"));
  g_sf_1p2to2p1.reset((TGraphAsymmErrors*)file_1p2to2p1->Get("ScaleFactors"));
  g_sf_2p1to2p4.reset((TGraphAsymmErrors*)file_2p1to2p4->Get("ScaleFactors"));
}

bool MuonTriggerWeights::process(Event & event){

  if(event.isRealData){
    event.set(h_muo_weight, 1.);
    event.set(h_muo_weight_up, 1.);
    event.set(h_muo_weight_down, 1.);
    return true;
  }

  const Muon muon = event.muons->at(0);
  double eta = fabs(muon.eta());
  double pt = (muon.pt() >= 1000.0 ? 999.0 : muon.pt());
  if(eta > 2.4) throw runtime_error("In BstarToTWModules.cxx, MuonTriggerWeights.process(): Muon-|eta| > 2.4 is not supported at the moment.");
  if(pt < 30.) throw runtime_error("In BstarToTWModules.cxx, MuonTriggerWeights.process(): Muon-pt < 30 is not supported at the moment.");

  // find number of correct pt bin
  int idx = 0;
  if(eta < 0.9) {
    bool keep_going = true;
    while(keep_going) {
      double x,y;
      g_sf_0to0p9->GetPoint(idx,x,y);
      keep_going = (pt > x + g_sf_0to0p9->GetErrorXhigh(idx));
      if(keep_going) idx++;
    }
  }
  else if(eta < 1.2){
    bool keep_going = true;
    while(keep_going){
      double x,y;
      g_sf_0p9to1p2->GetPoint(idx,x,y);
      keep_going =  (pt > x + g_sf_0p9to1p2->GetErrorXhigh(idx));
      if(keep_going) idx++;
    }
  }
  else if(eta < 2.1){
    bool keep_going = true;
    while(keep_going){
      double x,y;
      g_sf_1p2to2p1->GetPoint(idx,x,y);
      keep_going =  (pt > x + g_sf_1p2to2p1->GetErrorXhigh(idx));
      if(keep_going) idx++;
    }
  }
  else {
    bool keep_going = true;
    while(keep_going){
      double x,y;
      g_sf_2p1to2p4->GetPoint(idx,x,y);
      keep_going =  (pt > x + g_sf_2p1to2p4->GetErrorXhigh(idx));
      if(keep_going) idx++;
    }
  }

  //access scale factors, add 2% t&p systematic uncertainty
  double sf, sf_up, sf_down, dummy_x;
  double stat_up = -1., stat_down = -1., tp = 0.0, total_up = -1., total_down = -1.;
  if(eta < 0.9){
    g_sf_0to0p9->GetPoint(idx,dummy_x,sf);

    stat_up = g_sf_0to0p9->GetErrorYhigh(idx);
    stat_down = g_sf_0to0p9->GetErrorYlow(idx);
    total_up = sqrt(pow(stat_up,2) + pow(tp,2));
    total_down = sqrt(pow(stat_down,2) + pow(tp,2));

    sf_up = sf + total_up;
    sf_down = sf - total_down;
  }
  else if(eta < 1.2){
    g_sf_0p9to1p2->GetPoint(idx,dummy_x,sf);

    stat_up = g_sf_0p9to1p2->GetErrorYhigh(idx);
    stat_down = g_sf_0p9to1p2->GetErrorYlow(idx);
    total_up = sqrt(pow(stat_up,2) + pow(tp,2));
    total_down = sqrt(pow(stat_down,2) + pow(tp,2));

    sf_up = sf + total_up;
    sf_down = sf - total_down;
  }
  else if(eta < 2.1){
    g_sf_1p2to2p1->GetPoint(idx,dummy_x,sf);

    stat_up = g_sf_1p2to2p1->GetErrorYhigh(idx);
    stat_down = g_sf_1p2to2p1->GetErrorYlow(idx);
    total_up = sqrt(pow(stat_up,2) + pow(tp,2));
    total_down = sqrt(pow(stat_down,2) + pow(tp,2));

    sf_up = sf + total_up;
    sf_down = sf - total_down;
  }
  else {
    g_sf_2p1to2p4->GetPoint(idx,dummy_x,sf);

    stat_up = g_sf_2p1to2p4->GetErrorYhigh(idx);
    stat_down = g_sf_2p1to2p4->GetErrorYlow(idx);
    total_up = sqrt(pow(stat_up,2) + pow(tp,2));
    total_down = sqrt(pow(stat_down,2) + pow(tp,2));

    sf_up = sf + total_up;
    sf_down = sf - total_down;
  }

  if (syst_direction == "up") 
    {
      event.weight *= sf_up;
    } 
  else if (syst_direction == "down") 
    {
      event.weight *= sf_down;
    } 
  else 
    {
      event.weight *= sf;
    }

  event.set(h_muo_weight, sf);
  event.set(h_muo_weight_up, sf_up);
  event.set(h_muo_weight_down, sf_down);
  return true;
}


// reimplementation from common modules TopPtReweight
TopPtReweighting::TopPtReweighting(uhh2::Context& ctx,
				   float a, float b,
				   const std::string& syst_a,
				   const std::string& syst_b,
				   const std::string& ttgen_name,
				   const std::string& weight_name):
  a_(a), b_(b),
  ttgen_name_(ttgen_name){
  
  weight_name_ = weight_name;
  if(!weight_name_.empty())
    h_weight_= ctx.get_handle<float>(weight_name);
  version_ = ctx.get("dataset_version", "");
  boost::algorithm::to_lower(version_);
  if(!ttgen_name_.empty()){
    h_ttbargen_ = ctx.get_handle<TTbarGen>(ttgen_name);
  }

  if (syst_a == "up")
    a_ *= 1.5;
  else if (syst_a == "down")
    a_ *= 0.5;

  if (syst_b == "up")
    b_ *= 1.5;
  else if (syst_b == "down")
    b_ *= 0.5;
}

bool TopPtReweighting::process(uhh2::Event& event){
  if (event.isRealData || (!boost::algorithm::contains(version_,"ttbar") && !boost::algorithm::contains(version_,"ttjets") && !boost::algorithm::starts_with(version_,"tt")) ) {
    return true;
  }
  const TTbarGen& ttbargen = !ttgen_name_.empty() ? event.get(h_ttbargen_) : TTbarGen(*event.genparticles,false);
  float wgt = 1.;
  if (ttbargen.DecayChannel() != TTbarGen::e_notfound) {
    float tpt1 = ttbargen.Top().v4().Pt();
    float tpt2 = ttbargen.Antitop().v4().Pt();
    wgt = sqrt(exp(a_+b_*tpt1)*exp(a_+b_*tpt2));
  }

  if(!weight_name_.empty())
    event.set(h_weight_, wgt);

  event.weight *= wgt;
  return true;
}
