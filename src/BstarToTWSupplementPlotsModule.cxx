#include "UHH2/BstarToTW/include/BstarToTWSupplementPlotsModule.h"
#include "UHH2/HOTVR/include/HOTVRHists.h"


BstarToTWSupplementPlotsModule::BstarToTWSupplementPlotsModule(Context & ctx) {

  double top_fpt_max   = 0.8;                                    // maximum pt fraction of leading subjet
  double top_m_min     = 0.0;                                    // minimum topjet mass
  double top_m_max     = std::numeric_limits<float>::infinity(); // maximum topjet mass
  double top_mpair_min = 50.;                                    // minimum pairwise mass of first three subjets
  double top_tau32_max = 0.56;                                   // maximum nsubjetiness tau_3/2
  BTag::algo btag_algo = BTag::DEEPJET; // b tag algortihm
  BTag::wp btag_wp = BTag::WP_MEDIUM;   // b tag working point
  // b tag id
  JetId id_btag = BTag(btag_algo, btag_wp);
  // top tag id without mass window cut
  TopJetId id_toptag = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), Tau32Groomed(top_tau32_max));

  // --- b tag categories
  // 0 b tag
  sel_0btag.reset(new NJetSelection(0, 0, id_btag));
  hist_0btag.reset(new AndHists(ctx, "0btag"));
  hist_0btag->add_hist(new HOTVRMatchingHists(ctx, "0btag_topmatch"));
  hist_0btag->add_hist(new HOTVRHists(ctx, "0btag_HOTVR_tagged", id_toptag));
  // 1 b tag
  sel_1btag.reset(new NJetSelection(1, 1, id_btag));
  hist_1btag.reset(new AndHists(ctx, "1btag"));
  hist_1btag->add_hist(new HOTVRMatchingHists(ctx, "1btag_topmatch"));
  hist_1btag->add_hist(new HOTVRHists(ctx, "1btag_HOTVR_tagged", id_toptag));
  // 1+ btag
  sel_1plusbtag.reset(new NJetSelection(1, -1, id_btag));
  hist_1plusbtag.reset(new AndHists(ctx, "1plusbtag"));
  hist_1plusbtag->add_hist(new HOTVRMatchingHists(ctx, "1plusbtag_topmatch"));
  hist_1plusbtag->add_hist(new HOTVRHists(ctx, "1plusbtag_HOTVR_tagged", id_toptag));
  // 2  btag
  sel_2btag.reset(new NJetSelection(2, -1, id_btag));
  hist_2btag.reset(new AndHists(ctx, "2btag"));
  hist_2btag->add_hist(new HOTVRMatchingHists(ctx, "2btag_topmatch"));
  hist_2btag->add_hist(new HOTVRHists(ctx, "2btag_HOTVR_tagged", id_toptag));

  // --- top tag (without mass window)
  sel_1toptag.reset(new NTopJetSelection(1, -1, id_toptag));
  hist_1toptag.reset(new AndHists(ctx, "1toptag_no_masswindow"));
  hist_1toptag->add_hist(new HOTVRMatchingHists(ctx, "1toptag_no_masswindow_topmatch"));
  hist_1toptag->add_hist(new HOTVRHists(ctx, "1toptag_no_masswindow_HOTVR_tagged", id_toptag));


}

bool BstarToTWSupplementPlotsModule::process(Event & event) {
    

  // --- b tag categories
  if (sel_0btag->passes(event))
    hist_0btag->fill(event);
  else if (sel_1plusbtag->passes(event))
    {
      hist_1plusbtag->fill(event);
      if (sel_1btag->passes(event))
	hist_1btag->fill(event);
      else if (sel_2btag->passes(event))
	hist_2btag->fill(event);
    }
    
  // --- top tag (without mass window)
  if (sel_1toptag->passes(event))
    hist_1toptag->fill(event);

  return true;
}

