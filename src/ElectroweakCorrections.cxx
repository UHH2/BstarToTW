#include "UHH2/BstarToTW/include/ElectroweakCorrections.h"

using namespace std;
using namespace uhh2;

ElectroweakCorrections::ElectroweakCorrections(Context &ctx) {
  h_ewk_weight = ctx.declare_event_output<float>("weight_ewk");
  TString dataset_name = ctx.get("dataset_version");
  TString file_path = ctx.get("ElectroWeakCorrections","/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/BstarToTW/data/EWKCorrections.root");
  TFile *f = new TFile(file_path,"READ");
  
  // determine V (W or Z)
  if (dataset_name.Contains("WJets"))
    {
      m_pdg_id = 24;
      m_fun = (TF1*)(f->Get("w_ewkcorr/w_ewkcorr_func")->Clone());;
    }
  else if (dataset_name.Contains("DYJets"))
    {
      m_pdg_id = 23;
      m_fun = (TF1*)(f->Get("w_ewkcorr/w_ewkcorr_func")->Clone());;      
    }
  f->Close();
}

bool ElectroweakCorrections::process(Event &event) {
  if (m_pdg_id == 0)
    {
      event.set(h_ewk_weight, 1.0);
      return true;
    }

  GenParticle v = get_genp(event);

  // check if V was actually found
  if ( abs(v.pdgId()) != m_pdg_id ) 
    {
      event.set(h_ewk_weight, 1.0);
      return true;
    }
  // reweight
  double ewk_correction_weight = m_fun->Eval(v.pt());

  event.weight *= ewk_correction_weight;
  event.set(h_ewk_weight,ewk_correction_weight);
  return true;
}

GenParticle ElectroweakCorrections::get_genp(Event &event) {
  GenParticle v_cand;
  for (const auto & genp: *event.genparticles)
    {
      // get genparticle with matching pdg ID of V
      if (abs(genp.pdgId()) == m_pdg_id) v_cand = genp;
    }
  return v_cand;
}
