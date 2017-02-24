#include "UHH2/BstarToTW/include/BstarToTWSelections.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/Jet.h"
#include "UHH2/common/include/Utils.h"

#include <stdexcept>
#include <vector>

using namespace uhh2;
using namespace std;

NHotvrSelection::NHotvrSelection(unsigned int n_min_, unsigned int n_max_, double pt_min_, double eta_max_) {
  n_min   = n_min_;
  n_max   = n_max_;
  pt_min  = pt_min_;
  eta_max = eta_max_;
}

bool NHotvrSelection::passes(const Event &event) {
  vector<TopJet> topjets = *event.topjets;
  unsigned int n = 0;
  for (TopJet topjet : topjets)
    { 
      double pt = topjet.v4().pt();
      double eta = abs(topjet.v4().eta());
      if ( (pt_min < pt) && (eta < eta_max) ) ++n;
    }
  return (n_min <= n && n <= n_max);
}

METSelection::METSelection(double met_min_) {
  met_min = met_min_;
}

bool METSelection::passes(const Event &event) {
  return (met_min < event.met->pt());
}
