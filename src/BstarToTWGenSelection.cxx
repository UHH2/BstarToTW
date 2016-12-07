#include "UHH2/BstarToTW/include/BstarToTWGenSelection.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/Particle.h"


using namespace uhh2;
using namespace std;
using namespace fastjet;

/* ------------------------------------------------------------------
 * IsSemiLeptonic - checks if the event was SemiLeptonic, i.e. either
 * the top or the W decayed leptonicly
 */
IsSemiLeptonic::IsSemiLeptonic(Context & ctx):
  h_BstarToTWGen(ctx.get_handle<BstarToTWGen>("BstarToTWgen")){}

bool IsSemiLeptonic::passes(const Event & event){
  const auto & BstarToTWGen = event.get(h_BstarToTWGen);
  return (BstarToTWGen.IsSemiLeptonicDecay());
}

/* ------------------------------------------------------------------
 * IsTopHad - checks if the top decayed hadronicly
 */
IsTopHad::IsTopHad(Context & ctx):
  h_BstarToTWGen(ctx.get_handle<BstarToTWGen>("BstarToTWgen")){}

bool IsTopHad::passes(const Event & event){
  const auto & BstarToTWGen = event.get(h_BstarToTWGen);
  return (BstarToTWGen.IsTopHadronicDecay());
}

// NMuo::NMuo(Context & ctx, unsigned int n_min, unsigned int n_max, double pt_min):
//   h_genprod(ctx.get_handle<BstarToTWGen>("BstarToTWgen")),
//   _n_min(n_min),
//   _n_max(n_max),
//   _pt_min(pt_min){}

// bool NMuo::passes(const Event & event)
// {
//   unsigned int n = 0;
//   const auto & genproducer = event.get(h_genprod);
//   vector<GenParticle> stable = genproducer.getStabel();
	    
/* ------------------------------------------------------------------
 * NJet - checks if there are n >= n_min ak4-jets with pt > pt_min in
 * the event.
 * 
 * n_min: minimum number of ak4-jets required.
 * pt_min: minimum pt of each ak4-jet.
 */
NJet::NJet(Context & ctx, unsigned int n_min, double pt_min):
  h_genjetcluster(ctx.get_handle<GenJetCluster>("BstarToTWgenjet")),
  _n_min(n_min),
  _pt_min(pt_min){}

bool NJet::passes(const Event & event)
{
  const auto & genjetcluster = event.get(h_genjetcluster);
  unsigned int n = 0;
  vector<PseudoJet> jets = genjetcluster.getAK4Jets();
  for (unsigned int i = 0; i < jets.size(); ++i)
    {
      if (jets[i].pt() > _pt_min) ++n;
    }
  
  return (n >= _n_min);
}

/* ------------------------------------------------------------------
 * NCMSTopJet - checks if there are n_min <= n <= n_max jets, that are
 * top tagged by the CMS top tagger in the event.
 * 
 * n_min: minimum number of top tagged jets
 * n_max: maximum number of top tagged jets
 */
NCMSTopJet::NCMSTopJet(Context & ctx, unsigned int n_min, unsigned int n_max):
  h_genjetcluster(ctx.get_handle<GenJetCluster>("BstarToTWgenjet")),
  _n_min(n_min),
  _n_max(n_max){}

bool NCMSTopJet::passes(const Event & event)
{
  const auto & genjetcluster = event.get(h_genjetcluster);
  vector<PseudoJet> jets = genjetcluster.getCMSTopTagged();
  
  return (jets.size() >= _n_min && jets.size() <= _n_max);
}

/* ------------------------------------------------------------------
 * NHOTVRTopJet - checks if there are n_min <= n <= n_max jets, that
 * are top tagged by the HOTVR top tagger in the event.
 * 
 * n_min: minimum number of top tagged jets
 * n_max: maximum number of top tagged jets
 */
NHOTVRTopJet::NHOTVRTopJet(Context & ctx, unsigned int n_min, unsigned int n_max):
  h_genjetcluster(ctx.get_handle<GenJetCluster>("BstarToTWgenjet")),
  _n_min(n_min),
  _n_max(n_max){}

bool NHOTVRTopJet::passes(const Event & event)
{
  const auto & genjetcluster = event.get(h_genjetcluster);
  vector<PseudoJet> jets = genjetcluster.getHOTVRTopTagged();
  
  return (jets.size() >= _n_min && jets.size() <= _n_max);
}
