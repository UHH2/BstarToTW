#include "UHH2/BstarToTW/include/HOTVRPerformanceSelections.h"
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

ptHotvrTopSelection::ptHotvrTopSelection(double pt_min_): pt_min(pt_min_) {}

bool ptHotvrTopSelection::passes(const Event & event) {
  assert(event.topjets); // check if topjets are properly read in
  if (event.topjets->at(0).pt() > pt_min) return true;

  return false;
}

NHotvrGenTopSelection::NHotvrGenTopSelection(unsigned int n_): n(n_) {}

bool NHotvrGenTopSelection::passes(const Event & event) {
  assert(event.gentopjets); // check if topjets are properly read in
  return (event.gentopjets->size() == n);
}
