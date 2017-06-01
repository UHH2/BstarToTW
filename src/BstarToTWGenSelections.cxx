#include "UHH2/BstarToTW/include/BstarToTWGenSelections.h"

using namespace uhh2;

SemiLepSelection::SemiLepSelection(Context &ctx) {
  h_BstarToTWGen = ctx.get_handle<BstarToTWGen>("BstarToTWgen");
}
bool SemiLepSelection::passes(const Event &event) {
  BstarToTWGen gen = event.get(h_BstarToTWGen);
  return gen.IsSemiLeptonicDecay();
}

MuonChannelSelection::MuonChannelSelection(Context &ctx) {
  h_BstarToTWGen = ctx.get_handle<BstarToTWGen>("BstarToTWgen");
}
bool MuonChannelSelection::passes(const Event &event) {
  BstarToTWGen gen = event.get(h_BstarToTWGen);
  return (gen.IsMuonDecay() && gen.IsTopHadronicDecay());
}
