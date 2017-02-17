#include "UHH2/BstarToTW/include/AndHists.h"

using namespace std;
using namespace uhh2;
AndHists::AndHists() {}

void AndHists::fill(const Event & event) {  
  for (Hists *hist : hists_vector)
    {
      hist->fill(event);
    }
}

void AndHists::add_hist(Hists *hist) {
  hists_vector.push_back(hist);
}

AndHists::~AndHists() {
  hists_vector.clear();
}
