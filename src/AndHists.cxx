#include "UHH2/BstarToTW/include/AndHists.h"
#include "UHH2/BstarToTW/include/BstarToTWHists.h"
#include "UHH2/BstarToTW/include/HOTVRHists.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/JetHists.h"


using namespace std;
using namespace uhh2;
AndHists::AndHists(Context &ctx, const string & dirname, const boost::optional<TopJetId> &id_topjet):
  Hists(ctx, dirname+"_Counter")
{
  // Counting Hist
  nevt = book<TH1F>("NEvt", "", 1,0,1);

  // Add common hists to vector
  hists_vector.push_back(new BstarToTWHists(ctx, dirname + "_BstarToTW"));
  if (id_topjet)  hists_vector.push_back(new HOTVRHists(ctx, dirname + "_HOTVR", *id_topjet));
  else hists_vector.push_back(new HOTVRHists(ctx, dirname + "_HOTVR"));
  hists_vector.push_back(new LuminosityHists(ctx, dirname + "_Lumi"));
  hists_vector.push_back(new EventHists(ctx, dirname + "_Event"));
  hists_vector.push_back(new MuonHists(ctx, dirname + "_Muon"));
  hists_vector.push_back(new ElectronHists(ctx, dirname + "_Electron"));
  hists_vector.push_back(new JetHists(ctx, dirname + "_Jet"));
}

void AndHists::fill(const Event & event) {  
  nevt->Fill(0., event.weight);
  for (Hists *hist : hists_vector)
    {
      hist->fill(event);
    }
}

void AndHists::add_hist(Hists *hist) { 
  hists_vector.push_back(hist);
}

AndHists::~AndHists() {
  for (unsigned int i = 0; i < hists_vector.size(); ++i)
    {
      delete hists_vector[i];
    }
  hists_vector.clear();
}
