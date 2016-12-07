#include "UHH2/BstarToTW/include/GenJetCluster.h"

using namespace std;
using namespace uhh2;
using namespace fastjet;
using namespace contrib;

/* ------------------------------------------------------------------
 * GenJetCluster - Clusters GenParticles to PseudoJets with different
 * jet-algorithms. It can be used on it's own or in Analysis Module
 * with GenJetProducer.
 */
GenJetCluster::GenJetCluster(const Event & event) 
{
  vector<GenParticle> allGenp = *event.genparticles;

  // ------------------------------------------------------------------
  // get all stable particles except neutrinos and prepare for
  // clustering
  for(unsigned int i = 0; i < allGenp.size(); ++i)
    {
      const auto & genp = allGenp[i];
      if( (genp.status() == 1) && (abs(genp.pdgId()) != 12 && abs(genp.pdgId()) != 14 && abs(genp.pdgId()) != 16) )
	{
	  _stablePsj.push_back(convertGenParticle(genp));
	}
      // including 'hard' leptons results in double-counting, but excluding them results in less genjets then expected
      // else if ( genp.status() == 23 && (abs(genp.pdgId()) == 11 || abs(genp.pdgId()) == 13 || abs(genp.pdgId()) == 15 ) )
      // 	{ 
      // 	  _stableGenp.push_back(genp);
      // 	  _stablePsj.push_back(convertGenParticle(genp));
      // 	}
    }

  // ------------------------------------------------------------------
  // cluster anit-kt jets
  double ak4_pt_min = 10;  // ak4 jets are required to have pt > 10 GeV
  double ak8_pt_min = 150; // ak8 jets are required to have pt > 150 GeV

  // ak4 jets
  JetDefinition jet_def_ak4(antikt_algorithm, 0.4);
  ClusterSequence cs_ak4(_stablePsj, jet_def_ak4);
  _ak4jets = sorted_by_pt(cs_ak4.inclusive_jets(ak4_pt_min));

  // ak8 jets
  JetDefinition jet_def_ak8(antikt_algorithm, 0.8);
  ClusterSequence cs_ak8(_stablePsj, jet_def_ak8);
  _ak8jets = sorted_by_pt(cs_ak8.inclusive_jets(ak8_pt_min));

  // ------------------------------------------------------------------
  // Top Tagger

  // Top tagger using the soft-drop algorithm. Settings are as described in CMS PAS B2G-15-003
  double y_max = 2.4;          // maximum rapidity of top-candidate
  double beta_softdrop  = 0.0;
  double beta_nsubj = 1.0;
  double z_cut = 0.1;
  double R0 = 0.8;
  double cms_pt_min = 400;

  SoftDrop sd(beta_softdrop, z_cut, R0);
  NsubjettinessRatio   nSub32(3,2, OnePass_WTA_KT_Axes(), UnnormalizedMeasure(beta_nsubj));
  for (unsigned int i = 0; i < _ak8jets.size(); ++i)
    {
      if (_ak8jets[i].pt() > cms_pt_min && abs(_ak8jets[i].rap()) < y_max)
	{
	  
	  PseudoJet sdjet = sd(_ak8jets[i]);
	  double tau32 = nSub32(sdjet);
	  if ((sdjet.m() > 110 && sdjet.m() < 210) && (tau32 < 0.69))
	    {
	      _cmstoptagged.push_back(sdjet);
	    }
	}
    }

  // Top tagger using the HOTVR algorithm. Settings are as described in arXiv:1606.04961v1
  double mu(30.),      // massjump threshold
    theta(0.7),        // massjump parameter
    max_r(1.5),        // maximum allowed distance R
    min_r(0.1),        // minimum allowed distance R
    rho(600.),         // cone shrinking parameter
    hotvr_pt_min(30.); // minimum pT of subjets

  HOTVR hotvr_plugin(mu, theta, min_r, max_r, rho, hotvr_pt_min, HOTVR::CALIKE); 
  JetDefinition jet_def(&hotvr_plugin);
  ClusterSequence cs(_stablePsj, jet_def);
  vector<PseudoJet> hotvr_jets = hotvr_plugin.get_jets();
  for (unsigned int i = 0; i < hotvr_jets.size(); ++i)
    {
      HOTVRinfo hi = hotvr_jets[i].user_info<HOTVRinfo>();
      vector<PseudoJet> subjets = hi.subjets();
      if ( (hi.ptfraction(1) < 0.8) && (subjets.size() >= 3) && (hotvr_jets[i].m() > 140 && hotvr_jets[i].m() < 220) && (hi.mmin() > 50))
	{
	  _hotvrtoptagged.push_back(hotvr_jets[i]);
	}
	
    }
}

// get all stable genparticles except neutrinos as pseudojets.
vector<PseudoJet> GenJetCluster::getPseudojets() const
{
  return _stablePsj;
}

// get ak4 Jets
vector<PseudoJet> GenJetCluster::getAK4Jets() const
{
  return _ak4jets;
}

// get ak8 Jets
vector<PseudoJet> GenJetCluster::getAK8Jets() const
{
  return _ak8jets;
}

// get CMS top-tagged Jets
vector<PseudoJet> GenJetCluster::getCMSTopTagged() const
{
  return _cmstoptagged;
}

// get HOTVR top-tagged Jets
vector<PseudoJet> GenJetCluster::getHOTVRTopTagged() const
{
  return _hotvrtoptagged;
}

// convert genparticle to pseudojet.
PseudoJet GenJetCluster::convertGenParticle(const GenParticle & genp)
{
  LorentzVector v = genp.v4();

  return PseudoJet(v.Px(), v.Py(), v.Pz(), v.E());
}

/* ------------------------------------------------------------------
 * GenJetProducer - Creates handle for PseudoJets in event.
 */
GenJetProducer::GenJetProducer(Context & ctx, const string & name)
{ 
  h_genjetcluster = ctx.get_handle<GenJetCluster>(name);
}

bool GenJetProducer::process(Event & event)
{
  event.set(h_genjetcluster, GenJetCluster(event));
  return true;
}
