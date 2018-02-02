#include "UHH2/BstarToTW/include/GenCleaningModules.h"
#include "UHH2/common/include/Utils.h"

using namespace uhh2;
using namespace std;

GenTopJetCleaner::GenTopJetCleaner(uhh2::Context & ctx, GenTopJetId gentopjetid):
  m_gentopjetid(gentopjetid) { }

bool GenTopJetCleaner::process(Event &event) {
  
  clean_collection(*event.gentopjets, event, m_gentopjetid);
  return true;
}

GenJetCleaner::GenJetCleaner(uhh2::Context & ctx, double pt, double eta):
  m_pt(pt),
  m_eta(eta) { }

bool GenJetCleaner::process(Event &event) {
  vector<Particle> & genjets = *event.genjets;
  vector<Particle> result;
  for (Particle & genjet : genjets)
    {
      if (genjet.pt() > m_pt && genjet.eta() < m_eta) result.push_back(genjet);
    }
  std::swap(result, genjets);
    
  return true;
}
