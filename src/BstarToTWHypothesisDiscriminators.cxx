#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/common/include/Utils.h"

#include <set>

using namespace uhh2;
using namespace std;


namespace {
    
  // invariant mass of a lorentzVector, but safe for timelike / spacelike vectors
  float inv_mass(const LorentzVector& p4){ return p4.isTimelike() ? p4.mass() : -sqrt(-p4.mass2()); }

}

const BstarToTWHypothesis * get_best_hypothesis(const std::vector<BstarToTWHypothesis> & hyps, const std::string & label) {

  const BstarToTWHypothesis * best = nullptr;
  float current_best_disc = numeric_limits<float>::infinity();

  for(const auto & hyp : hyps)
    {
      if(!hyp.has_discriminator(label)) continue;
      auto disc = hyp.get_discriminator(label);
      if(disc < current_best_disc)
	{
	  best = &hyp;
	  current_best_disc = disc;
        }
    }
  if(std::isfinite(current_best_disc))
    {
      return best;
    }
  else
    {
      return nullptr;
    }
}

const LeptonicTopHypothesis * get_best_hypothesis(const std::vector<LeptonicTopHypothesis> & hyps, const std::string & label) {

  const LeptonicTopHypothesis * best = nullptr;
  float current_best_disc = numeric_limits<float>::infinity();

  for(const auto & hyp : hyps)
    {
      if(!hyp.has_discriminator(label)) continue;
      auto disc = hyp.get_discriminator(label);
      if(disc < current_best_disc)
	{
	  best = &hyp;
	  current_best_disc = disc;
        }
    }
  if(std::isfinite(current_best_disc))
    {
      return best;
    }
  else
    {
      return nullptr;
    }
}


BstarToTWChi2Discriminator::BstarToTWChi2Discriminator(Context & ctx, const std::string & rechyps_name) {

  h_hyps = ctx.get_handle<vector<BstarToTWHypothesis>>(rechyps_name);
  h_bjets = ctx.get_handle<vector<Jet> >("btag_medium");
  // from fits to matched distribution
  m_mtop_mean  = 179.9;
  m_mtop_sigma =  17.7;

  // m_deltaPhi_mean = M_PI;
  // m_deltaPhi_sigma = 0.105;

  // m_deltaPt_mean = -19.988;
  // m_deltaPt_sigma = 54.090;
  // over p_T
  // m_deltaPt_mean = -0.028;
  // m_deltaPt_sigma = 0.081;

  m_deltaPhi_mean = M_PI;
  m_deltaPhi_sigma = 0.064;

  // m_deltaPt_mean = -19.36;
  // m_deltaPt_sigma = 53.89;

  // m_deltaPt_mean = 0.05;
  // m_deltaPt_sigma = 0.06;
  m_deltaPt_mean = 0.0;
  m_deltaPt_sigma = 0.046;
}

bool BstarToTWChi2Discriminator::process(uhh2::Event& event) {

  auto& hyps = event.get(h_hyps);

  for(auto& hyp : hyps)
    {
      double mtop_reco = inv_mass(hyp.get_topjet());
      const double chi2_mtop = pow((mtop_reco - m_mtop_mean) / m_mtop_sigma, 2);

      double deltaPhi_reco = (hyp.get_w().phi() - hyp.get_topjet().phi());
      if (deltaPhi_reco < 0) deltaPhi_reco += 2*M_PI;
      const double chi2_deltaPhi = pow((deltaPhi_reco - m_deltaPhi_mean) / m_deltaPhi_sigma, 2);

      // double deltaPt_reco = (hyp.get_topjet().pt() - hyp.get_w().pt());
      // const double chi2_deltaPt = pow((deltaPt_reco - m_deltaPt_mean) / m_deltaPt_sigma, 2);
      double deltaPt_reco = (hyp.get_topjet().pt() - hyp.get_w().pt()) / (hyp.get_topjet().pt() + hyp.get_w().pt());
      const double chi2_deltaPt = pow((deltaPt_reco - m_deltaPt_mean) / m_deltaPt_sigma, 2);

      double deltaPt_toplep = std::numeric_limits<double>::infinity();
      double dRmin = std::numeric_limits<double>::infinity();

      vector<Jet> jet_collection = (event.get(h_bjets).size() > 1) ? event.get(h_bjets) : *event.jets; // decide weather to take bjet or jet collection
      Jet *closest_jet = 0;
      for (Jet &j : jet_collection)
	{
	  double dR = deltaR(hyp.get_w(), j);
	  if (dR < dRmin)
	    {
	      dRmin = dR;
	      closest_jet = &j;
	    }
	}
      if (dRmin < 2.0) // make sure chosen jet does not overlap with top region
	{
	  deltaPt_toplep = (hyp.get_topjet().pt() - (hyp.get_w() + closest_jet->v4()).pt()) / (hyp.get_topjet().pt() + (hyp.get_w() + closest_jet->v4()).pt());
	}

      double deltaZ_lnu = abs(hyp.get_lepton().pz() - hyp.get_neutrino().pz());

      hyp.set_discriminator("Chi2_top", chi2_mtop);
      hyp.set_discriminator("Chi2_deltaPhi", chi2_deltaPhi);
      hyp.set_discriminator("Chi2_deltaPt", chi2_deltaPt);
      hyp.set_discriminator("Chi2", chi2_deltaPhi + chi2_deltaPt);
      hyp.set_discriminator("Chi2_with_mass", chi2_mtop + chi2_deltaPhi + chi2_deltaPt);
      hyp.set_discriminator("deltaPt_W",deltaPt_reco);
      hyp.set_discriminator("deltaPt_toplep",deltaPt_toplep);
      hyp.set_discriminator("closest_nu",deltaZ_lnu);
    }

  return true;
}


LeptonicTopChi2Discriminator::LeptonicTopChi2Discriminator(Context & ctx, const std::string & rechyps_name) {

  h_hyps = ctx.get_handle<vector<LeptonicTopHypothesis>>(rechyps_name);

  // from fits to matched distribution
  m_mtop_mean  = 172.5;
  m_mtop_sigma =  5.0;

}

bool LeptonicTopChi2Discriminator::process(uhh2::Event& event) {

  auto& hyps = event.get(h_hyps);

  for(auto& hyp : hyps)
    {
      // double mtophad = inv_mass(hyp.get_tophad());
      double mtoplep = inv_mass(hyp.get_toplep());
      const double chi2_mtoplep = pow((mtoplep - m_mtop_mean) / m_mtop_sigma, 2);

      hyp.set_discriminator("Chi2", chi2_mtoplep);
    }

  return true;
}


BstarToTWMatchDiscriminator::BstarToTWMatchDiscriminator(Context & ctx, const std::string & rechyps_name, const cfg & config_): config(config_){
  h_hyps = ctx.get_handle<vector<BstarToTWHypothesis>>(rechyps_name);
  h_bstartotwgen = ctx.get_handle<BstarToTWGen>(config.bstartotwgen_name);
}


bool BstarToTWMatchDiscriminator::process(uhh2::Event & event){
  auto & hyps = event.get(h_hyps);
  const auto & bstartotwgen = event.get(h_bstartotwgen);
  if(!bstartotwgen.IsSemiLeptonicDecay() || !bstartotwgen.IsTopHadronicDecay())
    {
      for(auto & hyp: hyps)
	{
	  hyp.set_discriminator(config.discriminator_label, infinity);
	}
      return true;
    }
  for(auto & hyp: hyps)
    {
      float deltaR_top = deltaR(bstartotwgen.tbstar(), hyp.get_topjet());
      if (deltaR_top >= 0.4)
	{
	  hyp.set_discriminator(config.discriminator_label, infinity);
	  continue;
	}

      float deltaR_lep = deltaR(bstartotwgen.ChargedLepton(), hyp.get_lepton());
      if (deltaR_lep >= 0.4)
	{
	  hyp.set_discriminator(config.discriminator_label, infinity);
	  continue;
	}

      float deltaR_neutrino = deltaR(bstartotwgen.Neutrino(), hyp.get_neutrino());
      if (deltaR_neutrino >= 0.4)
	{
	  hyp.set_discriminator(config.discriminator_label, infinity);
	  continue;
	}

      float correct_dr = deltaR_top + deltaR_lep + deltaR_neutrino;
      hyp.set_discriminator(config.discriminator_label, correct_dr);
    }
  return true;
}
