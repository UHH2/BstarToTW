#include "UHH2/BstarToTW/include/BstarToTWGenJetHists.h"
#include "TH1F.h"
#include "TH2F.h"
#include <vector>

using namespace uhh2;
using namespace std;
using namespace fastjet;

BstarToTWGenJetHists::BstarToTWGenJetHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  // Book histograms
  N_ak4Jets = book<TH1F>( "N_ak4Jets", "N_{jet, ak4}", 20, 0, 20 );
  N_ak8Jets = book<TH1F>( "N_ak8Jets", "N_{jet, ak8}", 20, 0, 20 );
  N_topJets_CMS = book<TH1F>( "N_topJets", "N_{jet, cms top}", 20, 0, 20 );
  N_topJets_HOTVR = book<TH1F>( "N_topJets_HOTVR", "N_{jet, hotvr top}", 20, 0, 20 );

  M_top_CMS = book<TH1F>( "M_top_CMS", "M_{top, CMS} [GeV/c^{2}]", 30, 100, 250 );
  M_top_HOTVR = book<TH1F>( "M_top_HOTVR", "M_{top, HOTVR} [GeV/c^{2}]", 30, 100, 250 );

  // Get container BstarToTWgen
  h_genjetcluster = ctx.get_handle<GenJetCluster>("BstarToTWgenjet");

}

void BstarToTWGenJetHists::fill(const uhh2::Event & event){
  // Do not fill histograms if BstarToTWgen information has not been filled
  if(!event.is_valid(h_genjetcluster)){
    return;
  }

  double eventweight = event.weight;

  // Get handles from BstarToTWGen
  const auto & genjetcluster = event.get(h_genjetcluster);

  vector<PseudoJet> ak4Genjets    = genjetcluster.getAK4Jets();
  vector<PseudoJet> ak8Genjets    = genjetcluster.getAK8Jets();
  vector<PseudoJet> cmsTopJets    = genjetcluster.getCMSTopTagged();
  vector<PseudoJet> hotvrTopJets  = genjetcluster.getHOTVRTopTagged();
  double mCMS;
  double mHOTVR;
  PseudoJet cmsTop;
  PseudoJet hotvrTop;
  for(unsigned int i = 0; i < cmsTopJets.size(); ++i)
    {
      cmsTop += cmsTopJets[i];
    }
  mCMS = cmsTop.m();

  for(unsigned int i = 0; i < hotvrTopJets.size(); ++i)
    {
      hotvrTop += hotvrTopJets[i];
    }
  mHOTVR = hotvrTop.m();

  N_ak4Jets->Fill(ak4Genjets.size(), eventweight);
  N_ak8Jets->Fill(ak8Genjets.size(), eventweight);
  N_topJets_CMS->Fill(cmsTopJets.size(), eventweight);
  N_topJets_HOTVR->Fill(hotvrTopJets.size(), eventweight);

  M_top_CMS->Fill(mCMS, eventweight);
  M_top_HOTVR->Fill(mHOTVR, eventweight); 
}
