#include "UHH2/core/include/Utils.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/BstarToTW/include/BstarToTWHists.h"


#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>


using namespace std;
using namespace uhh2;


bool operator<(const run_lumi & rl1, const run_lumi & rl2){
    if(rl1.run == rl2.run){
        return rl1.lumiblock < rl2.lumiblock;
    }
    else{
        return rl1.run < rl2.run;
    }
}

BstarToTWHists::BstarToTWHists(Context & ctx, const string & dirname): 
  Hists(ctx, dirname){

  lumi_per_bin = string2double(ctx.get("lumihists_lumi_per_bin", "50.0"));
  if(lumi_per_bin <= 0.0) {
    throw runtime_error("lumihists_lumi_per_bin is <= 0.0; this is not allowed");
  }
    
  string lumifile = ctx.get("lumi_file");
  std::unique_ptr<TFile> file(TFile::Open(lumifile.c_str(), "read"));
  TTree * tree = dynamic_cast<TTree*>(file->Get("AnalysisTree"));
  if(!tree){
    throw runtime_error("LuminosityHists: Did not find TTree 'AnalysisTree' in file ;" + lumifile + "'");
  }
  // only fetch branches we really need:
  TBranch * brun = tree->GetBranch("run");
  TBranch * blumiblock = tree->GetBranch("luminosityBlock");
  TBranch * bilumi = tree->GetBranch("intgRecLumi");
  run_lumi rl;
  double ilumi;
  brun->SetAddress(&rl.run);
  blumiblock->SetAddress(&rl.lumiblock);
  bilumi->SetAddress(&ilumi);
    
  // use the file to build a map from run/lumi --> integrated lumi in pb.
  // Inserting into a map sorts by run and lumi.
  std::map<run_lumi, double> rl2lumi;
  double total_lumi = 0.0; // in   1/pb
  auto ientries = tree->GetEntries();
  for(auto ientry = 0l; ientry < ientries; ientry++){
    for(auto b : {brun, blumiblock, bilumi}){
      b->GetEntry(ientry);
    }
    double ilumi_pb = ilumi * 1e-6; // convert units in file (microbarn) to pb.
    total_lumi += ilumi_pb;
    auto it_inserted = rl2lumi.insert(make_pair(rl, ilumi_pb));
    if(!it_inserted.second){
      throw runtime_error("Duplicate run/lumi entry in lumi file '" + lumifile + "'");
    }
  }
  //cout << "LuminosityHists: read total lumi " << total_lumi << " from lumi file " << lumifile << endl;
    
  // Save the bin borders to find out the number of bins to use and for later assigning each event to a bin.
  int nbins_estimated = int(total_lumi / lumi_per_bin + 1);
  if(nbins_estimated >= 20000){
    throw runtime_error("LuminosityHists misconfiguration: would take more than 20000 bins. Please increase lumi_per_bin");
  }
  upper_binborders.reserve(nbins_estimated);
  double ilumi_current_bin = 0.0;
  for(const auto & rl : rl2lumi){
    ilumi_current_bin += rl.second;
    if(ilumi_current_bin >= lumi_per_bin){
      upper_binborders.push_back(rl.first);
      ilumi_current_bin = ilumi_current_bin - lumi_per_bin;
    }
  }
  int nbins = upper_binborders.size() + 1; // add one for the partial bin

  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  // book all histograms here

  // Event variables
  MET = book<TH1F>("MET", "#slash{E}_{T} [GeV/c]", 20, 0, 1000);
  HT_lep = book<TH1F>("HT_lep", "HT_{lep} [GeV/c]", 70, 0, 3500);
  HT_jet = book<TH1F>("HT_jet", "HT_{jet} [GeV/c]", 70, 0, 3500);
  ST = book<TH1F>("ST", "S_{T} [GeV/c]", 80, 0, 4000);
  rho = book<TH1F>("rho", "#rho [GeV/c]", 12, 0, 65);
  deltaPhi_blep = book<TH1F>("deltaPhi_blep", "#Delta #phi_{b,l}", 35, 0, 3.5);
  deltaPhi_btop = book<TH1F>("deltaPhi_btop", "#Delta #phi_{b,t}", 35, 0, 3.5);

  LumiBlock_vs_NPV = book<TH2F>("LumiBlock_vs_NPV", "LumiBlock_vs_NPV",
				nbins, 0, ( int(total_lumi / lumi_per_bin) + 1)*lumi_per_bin,
				40, 0, 40);

    }


void BstarToTWHists::fill(const Event & event){
  double weight = event.weight;
  vector<Jet> jets = *event.jets;
  vector<Electron> electrons = *event.electrons;
  vector<Muon> muons = *event.muons;
  vector<TopJet> topjets = *event.topjets;
  
  
  double ht_lep = 0;
  double ht_jet = 0;

  const Particle &primlep = event.get(h_primlep);

  for (Electron ele : electrons)
    {
      ht_lep += ele.v4().pt();
    }

  for (Muon muo : muons)
    {
      ht_lep += muo.v4().pt();
    }
  MET->Fill(event.met->pt(), weight);
  ht_lep += event.met->pt();

  for (Jet jet : jets)
    {
      ht_jet += jet.v4().pt();
      if(btag_loose(jet, event))
	{
	  deltaPhi_blep->Fill(deltaPhi(jet.v4(),primlep.v4()), weight);
	  for (TopJet topjet : topjets)
	    {
	      deltaPhi_btop->Fill(deltaPhi(jet.v4(),topjet.v4()), weight);
	    }
	}
    }

  HT_lep->Fill(ht_lep, weight);
  HT_jet->Fill(ht_jet, weight);
  ST->Fill(ht_lep + ht_jet, weight);

  if(event.isRealData)
    {
      run_lumi rl{event.run, event.luminosityBlock};
      auto it = upper_bound(upper_binborders.begin(), upper_binborders.end(), rl);
      int ibin = distance(upper_binborders.begin(), it);
      LumiBlock_vs_NPV->Fill(ibin*lumi_per_bin, event.pvs->size(), weight);
    }
  rho->Fill(event.rho, weight);

}

BstarToTWHists::~BstarToTWHists(){}


BstarToTWBackgroundHists::BstarToTWBackgroundHists(Context & ctx, const string & dirname, const string & hyps_name, const TString & path):
  Hists(ctx, dirname){
  h_hyps = ctx.get_handle<std::vector<BstarToTWHypothesis>>(hyps_name);
  TFile* f = new TFile(path);
  TH1F* ratio = (TH1F*)(f->Get("ratio"))->Clone("ratio");
  TF1* linfit = ratio->GetFunction("linfit");
  m_p0 = linfit->GetParameter(0);
  m_p1 = linfit->GetParameter(1);
  //  m_p2 = linfit->GetParameter(2);
  m_p2 = 0;

  double xbins[17] = {0, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1700, 1900, 2100, 2400, 5000};
  double xbins_fitbin[13] = {0,500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100, 2300, 2500, 5000};
  Bstar_reco_M_unbinned = book<TH1F>("Bstar_reco_M_unbinned", "M_{tW} [GeV/c^{2}]", 25, 100, 5100);
  Bstar_reco_M_rebin = book<TH1F>("Bstar_reco_M_rebin", "M_{tW} [GeV/c^{2}]", 12, xbins_fitbin);
  Bstar_reco_M = book<TH1F>("Bstar_reco_M", "M_{tW} [GeV/c^{2}]", 16, xbins);
}

void BstarToTWBackgroundHists::fill(const Event & event){
  
  std::vector<BstarToTWHypothesis> hyps = event.get(h_hyps);
  const BstarToTWHypothesis* hyp = get_best_hypothesis( hyps, "Chi2" );
  if (!hyp)
    {
      //cout << "WARNING: " + m_hyps_name  +": " + m_discriminator_name + " No hypothesis was valid!" << endl;
      return;
    }
  
  double mbstar = 0;
  if((hyp->get_topjet() + hyp->get_w()).isTimelike())
    {    
      mbstar = (hyp->get_topjet() + hyp->get_w()).M();
    }
  else
    {
      mbstar = sqrt(-(hyp->get_topjet()+hyp->get_w()).mass2());
    }
  //  double background_weight = m_p0 + m_p1 * mbstar + m_p2 * mbstar * mbstar;
  double background_weight = exp(m_p0 + m_p1 * mbstar);
  if ( background_weight < 0) 
    background_weight = 0;
  const double event_weight = event.weight;
  background_weight *= event_weight;
  if (mbstar < 5000.) 
    {
      Bstar_reco_M->Fill(mbstar, background_weight);
      Bstar_reco_M_rebin->Fill(mbstar, background_weight);
      Bstar_reco_M_unbinned->Fill(mbstar, background_weight);
    }
  else 
    {
      Bstar_reco_M->Fill(4999., background_weight);
      Bstar_reco_M_rebin->Fill(4999., background_weight);
      Bstar_reco_M_unbinned->Fill(4999., background_weight);
    }

}

BstarToTWBackgroundHists::~BstarToTWBackgroundHists(){}
