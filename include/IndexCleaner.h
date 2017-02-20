#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Selection.h"

#include "UHH2/BstarToTW/include/TopTagIndexer.h"

namespace uhh2 {
  
  class PtTopIndexCleaner: public uhh2::AnalysisModule {
  public:
    PtTopIndexCleaner(uhh2::Context &ctx, double pt_min_);
    virtual bool process(uhh2::Event &event);
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double pt_min;
  };

  class MTopIndexCleaner: public uhh2::AnalysisModule {
  public:
    MTopIndexCleaner(uhh2::Context &ctx, double m_min_, double m_max_);
    virtual bool process(uhh2::Event &event);
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double m_min, m_max;
  };

  class EtaTopIndexCleaner: public uhh2::AnalysisModule {
  public:
    EtaTopIndexCleaner(uhh2::Context &ctx, double eta_max_);
    virtual bool process(uhh2::Event &event);
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double eta_max;
  };

  class NSubTopIndexCleaner: public uhh2::AnalysisModule {
  public:
    NSubTopIndexCleaner(uhh2::Context &ctx, unsigned int nsub_min_);
    virtual bool process(uhh2::Event &event);
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    unsigned int nsub_min;
  };

  class FptTopIndexCleaner: public uhh2::AnalysisModule {
  public:
    FptTopIndexCleaner(uhh2::Context &ctx, double fpt_max_);
    virtual bool process(uhh2::Event &event);
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double fpt_max;
  };

  class MpairTopIndexCleaner: public uhh2::AnalysisModule {
  public:
    MpairTopIndexCleaner(uhh2::Context &ctx, double mpair_min_);
    virtual bool process(uhh2::Event &event);
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double mpair_min;
  };

class Tau32TopIndexCleaner: public uhh2::AnalysisModule {
  public:
    Tau32TopIndexCleaner(uhh2::Context &ctx, double t32_max_);
    virtual bool process(uhh2::Event &event);
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double t32_max;
  };

class HotvrTopIndexCleaner: public uhh2::AnalysisModule {
  public:
    HotvrTopIndexCleaner(uhh2::Context &ctx, double pt_min_, double eta_max_, unsigned int nsub_min_, double fpt_max_, double mpair_min_, double tau32_max_);
    virtual bool process(uhh2::Event &event);
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double pt_min;
    double eta_max;
    unsigned int nsub_min;
    double fpt_max;
    double mpair_min;
    double tau32_max;
  };

  class NTopIndexSelection: public uhh2::Selection {
  public:
    NTopIndexSelection(uhh2::Context &ctx, unsigned int n_min_, unsigned int n_max_ = -1);
    virtual bool passes(const uhh2::Event &event);
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    unsigned int n_min, n_max;
  };

}
