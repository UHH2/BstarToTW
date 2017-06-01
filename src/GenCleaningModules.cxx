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
