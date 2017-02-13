#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/common/include/Utils.h"

#include <stdexcept>
#include <vector>

using namespace uhh2;
using namespace std;

NHotvrTopSelection::NHotvrTopSelection(unsigned int n_): n(n_) {}

bool NHotvrTopSelection::passes(const Event & event) {
  assert(event.topjets); // check if topjets are properly read in
  return (event.topjets->size() == n);
}

MHotvrTopSelection::MHotvrTopSelection(double m_min_, double m_max_):
  m_min(m_min_),
  m_max(m_max_) {}

bool MHotvrTopSelection::passes(const Event & event) {
  assert(event.topjets);
  double m = event.topjets->at(0).v4().M();
  return (m_min < m && m < m_max);
}

NHotvrGenTopSelection::NHotvrGenTopSelection(unsigned int n_): n(n_) {}

bool NHotvrGenTopSelection::passes(const Event & event) {
  assert(event.gentopjets); // check if topjets are properly read in
  return (event.gentopjets->size() == n);
}

