#include "UHH2/BstarToTW/include/BstarToTWCleaningModules.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/common/include/Utils.h"

#include <stdexcept>
#include <vector>

using namespace uhh2;
using namespace std;

HOTVRTopCleaner::HOTVRTopCleaner(unsigned int nsub_min_, double fpt_max_, double m_min_, double m_max_, double mpairwise_min_): 
  nsub_min(nsub_min_), 
  fpt_max(fpt_max_),
  m_min(m_min_),
  m_max(m_max_),
  mpairwise_min(mpairwise_min_){}

bool HOTVRTopCleaner::process(Event & event) {
  assert(event.topjets); // check if topjets are properly read in
  vector<TopJet> & topjets = *event.topjets;
  vector<TopJet> results;
  for (TopJet topjet : topjets)
    {
      // double m = topjet.v4().M();
      vector<Jet> subjets = topjet.subjets();
      if (subjets.size() < nsub_min) continue;
      double fpt = subjets.at(0).v4().pt()/topjet.v4().pt();
      double mpairwise = min(min((subjets.at(0).v4() + subjets.at(1).v4()).M(), (subjets.at(0).v4() + subjets.at(2).v4()).M()), (subjets.at(1).v4() + subjets.at(2).v4()).M());
      if (  fpt < fpt_max && mpairwise > mpairwise_min)	
	{
	  results.push_back(topjet);
	}
    }
  swap(results, *event.topjets);
  return true;
}

HOTVRGenTopCleaner::HOTVRGenTopCleaner(unsigned int nsub_min_, double fpt_max_, double m_min_, double m_max_, double mpairwise_min_): 
  nsub_min(nsub_min_), 
  fpt_max(fpt_max_),
  m_min(m_min_),
  m_max(m_max_),
  mpairwise_min(mpairwise_min_){}

bool HOTVRGenTopCleaner::process(Event & event) {
  assert(event.gentopjets); // check if gentopjets are properly read in
  vector<GenTopJet> & topjets = *event.gentopjets;
  vector<GenTopJet> results;
  for (GenTopJet topjet : topjets)
    {
      vector<Particle> subjets = topjet.subjets();
      if (subjets.size() < nsub_min) continue;
      double fpt = subjets.at(0).v4().pt()/topjet.v4().pt();
      double mpairwise = min(min((subjets.at(0).v4() + subjets.at(1).v4()).M(), (subjets.at(0).v4() + subjets.at(2).v4()).M()), (subjets.at(1).v4() + subjets.at(2).v4()).M());
      if (  fpt < fpt_max && mpairwise > mpairwise_min)	
	{
	  results.push_back(topjet);
	}
    }
  swap(results, *event.gentopjets);
  return true;
}
