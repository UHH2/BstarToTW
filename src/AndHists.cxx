#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/HOTVRHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/MuonHists.h"


using namespace std;
using namespace uhh2;
AndHists::AndHists(Context &ctx, const string & dirname) {
  // Add common hists to vector
  hists_vector.push_back(new HOTVRHists(ctx, dirname + "_HOTVR"));
  hists_vector.push_back(new EventHists(ctx, dirname + "_Event"));
  hists_vector.push_back(new MuonHists(ctx, dirname + "_Muon"));
}

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
