#include "UHH2/BstarToTW/include/BstarToTWHypothesisDiscriminators.h"
#include "UHH2/core/include/Utils.h"

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


BstarToTWChi2Discriminator::BstarToTWChi2Discriminator(Context & ctx, const std::string & rechyps_name) {

  h_hyps = ctx.get_handle<vector<BstarToTWHypothesis>>(rechyps_name);

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
  m_deltaPhi_sigma = 0.05458;

  // m_deltaPt_mean = -19.36;
  // m_deltaPt_sigma = 53.89;

  m_deltaPt_mean = -0.002863;
  m_deltaPt_sigma = 0.06937;

}

bool BstarToTWChi2Discriminator::process(uhh2::Event& event) {

  auto& hyps = event.get(h_hyps);

  for(auto& hyp : hyps)
    {
      double mtop_reco = inv_mass(hyp.get_topjet());
      const double chi2_mtop = pow((mtop_reco - m_mtop_mean) / m_mtop_sigma, 2);

      double deltaPhi_reco = (hyp.get_topjet().phi() - hyp.get_w().phi());
      if (deltaPhi_reco < 0) deltaPhi_reco += 2*M_PI;
      const double chi2_deltaPhi = pow((deltaPhi_reco - m_deltaPhi_mean) / m_deltaPhi_sigma, 2);

      // double deltaPt_reco = (hyp.get_topjet().pt() - hyp.get_w().pt());
      // const double chi2_deltaPt = pow((deltaPt_reco - m_deltaPt_mean) / m_deltaPt_sigma, 2);
      double deltaPt_reco = (hyp.get_topjet().pt() - hyp.get_w().pt()) / hyp.get_topjet().pt();
      const double chi2_deltaPt = pow((deltaPt_reco - m_deltaPt_mean) / m_deltaPt_sigma, 2);

      hyp.set_discriminator("Chi2_top", chi2_mtop);
      hyp.set_discriminator("Chi2_deltaPhi", chi2_deltaPhi);
      hyp.set_discriminator("Chi2_deltaPt", chi2_deltaPt);
      hyp.set_discriminator("Chi2", chi2_deltaPhi + chi2_deltaPt);
      hyp.set_discriminator("Chi2_with_mass", chi2_mtop + chi2_deltaPhi + chi2_deltaPt);
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
      if (deltaR_top >= 0.2)
	{
	  hyp.set_discriminator(config.discriminator_label, infinity);
	  continue;
	}

      float deltaR_lep = deltaR(bstartotwgen.ChargedLepton(), hyp.get_lepton());
      if (deltaR_lep >= 0.2)
	{
	  hyp.set_discriminator(config.discriminator_label, infinity);
	  continue;
	}

      float correct_dr = deltaR_top + deltaR_lep;

      //add deltaR between reconstructed and true neutrino
      correct_dr += deltaR(bstartotwgen.Neutrino(), hyp.get_neutrino());
      hyp.set_discriminator(config.discriminator_label, correct_dr);
    }
  return true;
}
