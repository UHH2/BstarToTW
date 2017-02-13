#include "UHH2/BstarToTW/include/TopTagIndexer.h"

using namespace std;
using namespace uhh2;

TopTagIndexer::TopTagIndexer(Event &event) {
  _index.clear();
  assert(event.topjets); // check if topjets are properly read in
  vector<TopJet> *topjets = event.topjets;
  for (unsigned int i = 0; i < topjets->size(); ++i)
    {
      _index.push_back(i);
    }
}

vector<int> TopTagIndexer::GetIndex() const {
  return _index;
}

void TopTagIndexer::SetIndex(vector<int> ind) {
  swap(_index, ind);
}

/*
 * Producer module for TopTagIndexer
 */
TopTagIndexerProducer::TopTagIndexerProducer(Context &ctx, const string name) {
  h_TopTagIndexer = ctx.get_handle<TopTagIndexer>(name);
}

bool TopTagIndexerProducer::process(Event &event) {
  event.set(h_TopTagIndexer, TopTagIndexer(event));
  return true;
}

