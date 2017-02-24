#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Selection.h"

#include "UHH2/BstarToTW/include/TopTagIndexer.h"

namespace uhh2 {
  // ------- Cleaning Modules
  /*
   * The cleaning modules remove indices from TopTagIndexer, if they
   * don't fulfil the given requirements. For details on the
   * requirements see the respective module. If no jets survive the
   * cleaning the module returns false.
   */

  class PtEtaTopIndexCleaner: public uhh2::AnalysisModule {
  public:
    PtEtaTopIndexCleaner(uhh2::Context &ctx, double pt_min_, double eta_max_);
    virtual bool process(uhh2::Event &event);
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double pt_min, eta_max;
  };

  class MTopIndexCleaner: public uhh2::AnalysisModule {
  public:
    MTopIndexCleaner(uhh2::Context &ctx, double m_min_, double m_max_);
    virtual bool process(uhh2::Event &event);
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double m_min, m_max;
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

  class DeltaPhiTopIndexCleaner: public uhh2::AnalysisModule {
  public:
    DeltaPhiTopIndexCleaner(uhh2::Context &ctx, double dphi_min_);
    virtual bool process(uhh2::Event &event);
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double dphi_min;
  };

  // Selection modules
  /*
   * The selection modules work the same way as the respective cleaning
   * modules, except that they don't change the indices from
   * TopJetIndexer. This makes the modules only applicable for studies
   * of those cuts, since they cannot be combined. (except
   * NTopIndexSelection)
   */
  class NTopIndexSelection: public uhh2::Selection {
  public:
    NTopIndexSelection(uhh2::Context &ctx, unsigned int n_min_, unsigned int n_max_ = -1);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    unsigned int n_min, n_max;
  };

  class PtEtaTopIndexSelection: public uhh2::Selection {
  public:
    PtEtaTopIndexSelection(uhh2::Context &ctx, double pt_min_, double eta_max_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double pt_min, eta_max;
  };

  class MTopIndexSelection: public uhh2::Selection {
  public:
    MTopIndexSelection(uhh2::Context &ctx, double m_min_, double m_max_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double m_min, m_max;
  };

  class NSubTopIndexSelection: public uhh2::Selection {
  public:
    NSubTopIndexSelection(uhh2::Context &ctx, unsigned int nsub_min_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    unsigned int nsub_min;
  };

  class FptTopIndexSelection: public uhh2::Selection {
  public:
    FptTopIndexSelection(uhh2::Context &ctx, double fpt_max_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double fpt_max;
  };

  class MpairTopIndexSelection: public uhh2::Selection {
  public:
    MpairTopIndexSelection(uhh2::Context &ctx, double mpair_min_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double mpair_min;
  };

  class Tau32TopIndexSelection: public uhh2::Selection {
  public:
    Tau32TopIndexSelection(uhh2::Context &ctx, double t32_max_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double t32_max;
  };

  class DeltaPhiTopIndexSelection: public uhh2::Selection {
  public:
    DeltaPhiTopIndexSelection(uhh2::Context &ctx, double dphi_min_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    uhh2::Event::Handle<TopTagIndexer> h_TopTagIndexer;
    double dphi_min;
  };

}
