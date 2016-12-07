#include "UniversalJetCluster.h"

using namespace std;
using namespace uhh2;
using namespace fastjet;

UniversalJetCluster::UniversalJetCluster(const vector<GenParticle> &genparticles)
{
  for(vector<GenParticle>::iterator it = genparticles.begin(); it != genparticles.end(); ++it) 
    {
      _psj.push_back(ConvertGenToPsj(*it));
    }
  ClusterHOTVR();
}

UniversalJetCluster::UniversalJetCluster(const vector<PFParticle> &pfparticles)
{
  for(vector<GenParticle>::iterator it = pfparticles.begin(); it != pfparticles.end(); ++it) 
    {
      _psj.push_back(ConvertPFToPsj(*it));
    }
  ClusterHOTVR();
}

void UniversalJetCluster::ClusterHOTVR()
{
  double mu(30.),                 // massjump threshold
    theta(0.7),                   // massjump parameter
    max_r(1.5),                   // maximum allowed distance R
    min_r(0.1),                   // minimum allowed distance R
    rho(600.),                    // cone shrinking parameter
    hotvr_pt_min(30.),            // minimum pT of subjets
    ptfraction_min(0.8),          // minimum pt,1/pt
    m_lower(140.),                // lower limit of mass window
    m_upper(220.),                // upper limit of mass window
    m_pairwise_min(50.);          // minimum pairwise mass
  unsigned int n_subjets_min = 3; // minimum number of subjets

  HOTVR hotvr_plugin(mu, theta, min_r, max_r, rho, hotvr_pt_min, HOTVR::CALIKE); 
  JetDefinition jet_def(&hotvr_plugin);
  ClusterSequence cs(_psj, jet_def);
  vector<PseudoJet> hotvr_jets = hotvr_plugin.get_jets();

  for (unsigned int i = 0; i < hotvr_jets.size(); ++i)
    {
      HOTVRinfo hi = hotvr_jets[i].user_info<HOTVRinfo>();
      vector<PseudoJet> subjets = hi.subjets();
      if (subjets.size() >= n_subjets_min)
	{
	  if ( (hi.ptfraction(1) < ptfraction_min) &&
	       (hotvr_jets[i].m() > m_lower && hotvr_jets[i].m() < m_upper) &&
	       (hi.mmin() >m_pairwise_min) )
	    {
	      _hotvrTop.push_back(ConvertPsjToTopJet(hotvr_jets[i]), subjets);
	    }
	}
    }
}

// ---------------------------------------------------------------
// Converters
// TODO: Write another converter for stable psj

PseudoJet UniversalJetCluster::ConvertGenToPsj(const GenParticle & genp);
{
  PseudoJet psj(genp.X(), genp.Y(), genp.Z(), genp.T());
  return psj;
}

PseudoJet UniversalJetCluster::ConvertPFToPsj(const PFParticle & pfp);
{
  PseudoJet psj(pfp.X(), pfp.Y(), pfp.Z(), pfp.T());
  return psj;
}

Jet UniversalJetCluster::ConvertPsjToJet(const PseudoJet & psj)
{
  Jet jet;
  jet.set_pt(psj.pt());
  jet.set_eta(psj.eta());
  jet.set_phi(psj.phi());
  jet.set_energy(psj.E());
  return jet;
}

TopJet UniversalJetCluster::ConvertPsjToTopJet(const PseudoJet & psj, const vector<PseudoJet> &subpsj)
{
  Topjet topjet;
  topjet.set_pt(psj.pt());
  topjet.set_eta(psj.eta());
  topjet.set_phi(psj.phi());
  topjet.set_energy(psj.E());
  for (vector<PseudoJet>::iterator it = subpsj.begin(); it != subpsj.end(); ++it) 
    {
      topjet.add_subjet(ConvertPsjToJet(*it));
    }
  return topjet;
}
