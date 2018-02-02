#pragma once

#include "UHH2/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/JetCorrections.h"


class FactorizedJetCorrector;

namespace JERFiles {

  extern const std::vector<std::string> Summer16_23Sep2016_V4_L1_AK4PFchs_MC;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L1_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L1_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L1_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L1_AK4PFchs_DATA;

  extern const std::vector<std::string> Summer16_23Sep2016_V4_L23_AK4PFchs_MC;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_BCD_L23_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_EF_L23_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_G_L23_AK4PFchs_DATA;
  extern const std::vector<std::string> Summer16_23Sep2016_V4_H_L23_AK4PFchs_DATA;
}

class HOTVRJetCorrector: public uhh2::AnalysisModule {
 public:
  explicit HOTVRJetCorrector(uhh2::Context &ctx, const std::vector<std::string> & filenames);

  virtual bool process(uhh2::Event &event) override;

  virtual ~HOTVRJetCorrector();

 private:
  std::unique_ptr<FactorizedJetCorrector> corrector;
  JetCorrectionUncertainty* jec_uncertainty;
  int direction = 0; // -1 = down, +1 = up, 0 = nominal

  double par[12][3]; // values with function parameters in 12 eta bins
  double pt_bins[7] = {0, 80, 130, 180, 250, 350, 500};
  double eta_bins[13] = {-4, -1.5, -1.0, -0.7, -0.4, -0.2, 0.0, 0.2, 0.4, 0.7, 1.0, 1.5, 4};

  double get_factor(double pt, double eta);
};
