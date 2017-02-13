#pragma once

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/TopJet.h"

#include <vector>

class TopTagIndexer {
 public:
  explicit TopTagIndexer(uhh2::Event &event);
  
  std::vector<int> GetIndex() const;
  void SetIndex(std::vector<int>);

 private:
  std::vector<int> _index;
};

class TopTagIndexerProducer : public uhh2::AnalysisModule {
 public:
  explicit TopTagIndexerProducer(uhh2::Context &ctx, const std::string name = "TopTagIndexer");
  virtual bool process(uhh2::Event &event) override;

 private:
  uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
};
