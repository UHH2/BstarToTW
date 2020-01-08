#include "UHH2/BstarToTW/include/BstarToTWSystematics.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/Utils.h"

#include <TRandomGen.h>

using namespace std;
using namespace uhh2;

MuonScaleFactors2016::MuonScaleFactors2016(Context &ctx) {

  m_sf_trigger.reset(new MCMuonScaleFactor(ctx,
					   ctx.get("SF_MUO_TRIGGER_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root"),
					   ctx.get("SF_MUO_TRIGGER_FILE", "IsoMu24_OR_IsoTkMu24_PtEtaBins"),
					   0.5,
					   "trigger",
					   false,
					   ctx.get("SYS_MUO_TRIGGER", "nominal")));
  
  m_sf_id.reset(new MCMuonScaleFactor(ctx,
				      ctx.get("SF_MUO_ID_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonID_EfficienciesAndSF_average_RunBtoH.root"),
				      ctx.get("SF_MUO_ID_FILE", "NUM_TightID_DEN_genTracks_eta_pt"),
				      1.0,
				      "tight_id",
				      false,
				      ctx.get("SYS_MUO_ID", "nominal")));

  m_sf_iso.reset(new MCMuonScaleFactor(ctx,
				       ctx.get("SF_MUO_ISO_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2016/MuonIso_EfficienciesAndSF_average_RunBtoH.root"),
				       ctx.get("SF_MUO_ISO_FILE", "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt"),
				       1.0,
				       "isolation",
				       false,
				       ctx.get("SYS_MUO_ISO", "nominal")));

}

bool MuonScaleFactors2016::process(Event &event) {

  m_sf_trigger->process(event);
  m_sf_id->process(event);
  m_sf_iso->process(event);

  return true;
}


MuonScaleFactors2017::MuonScaleFactors2017(Context &ctx) {

  m_sf_trigger.reset(new MCMuonScaleFactor(ctx,
					   ctx.get("SF_MUO_TRIGGER_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2017/MuonTrigger_EfficienciesAndSF_RunBtoF_Nov17Nov2017.root"),
					   ctx.get("SF_MUO_TRIGGER_FILE", "IsoMu27_PtEtaBins"),
					   0.5,
					   "trigger",
					   false,
					   ctx.get("SYS_MUO_TRIGGER", "nominal")));
  
  m_sf_id.reset(new MCMuonScaleFactor(ctx,
				      ctx.get("SF_MUO_ID_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2017/MuonID_94X_RunBCDEF_SF_ID.root"),
				      ctx.get("SF_MUO_ID_FILE", "NUM_TightID_DEN_genTracks_pt_abseta"),
				      1.0,
				      "tight_id",
				      true,
				      ctx.get("SYS_MUO_ID", "nominal")));

  m_sf_iso.reset(new MCMuonScaleFactor(ctx,
				       ctx.get("SF_MUO_ISO_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2017/MuonIso_94X_RunBCDEF_SF_ISO.root"),
				       ctx.get("SF_MUO_ISO_FILE", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"),
				       1.0,
				       "isolation",
				       true,
				       ctx.get("SYS_MUO_ISO", "nominal")));

}

bool MuonScaleFactors2017::process(Event &event) {

  m_sf_trigger->process(event);
  m_sf_id->process(event);
  m_sf_iso->process(event);

  return true;
}


MuonScaleFactors2018::MuonScaleFactors2018(Context &ctx, long int seed) {

  m_seed = seed;
  m_rng = new TRandomMixMax();
  m_rng->SetSeed(m_seed);

  m_sf_trigger_before.reset(new MCMuonScaleFactor(ctx,
						  ctx.get("SF_MUO_TRIGGER_BEFORE_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2018/Muon_Trigger_Eff_SF_BeforeMuonHLTUpdate.root"),
						  ctx.get("SF_MUO_TRIGGER_BEFORE_FILE", "IsoMu24_PtEtaBins"),
						  0.5,
						  "trigger",
						  false,
						  ctx.get("SYS_MUO_TRIGGER", "nominal")));

  m_sf_trigger_after.reset(new MCMuonScaleFactor(ctx,
						 ctx.get("SF_MUO_TRIGGER_AFTER_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2018/Muon_Trigger_Eff_SF_AfterMuonHLTUpdate.root"),
						 ctx.get("SF_MUO_TRIGGER_AFTER_FILE", "IsoMu24_PtEtaBins"),
						 0.5,
						 "trigger",
						 false,
						 ctx.get("SYS_MUO_TRIGGER", "nominal")));

  m_sf_id.reset(new MCMuonScaleFactor(ctx,
				      ctx.get("SF_MUO_ID_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2018/Muon_ID_SF_RunABCD.root"),
				      ctx.get("SF_MUO_ID_FILE", "NUM_TightID_DEN_TrackerMuons_pt_abseta"),
				      1.0,
				      "tight_id",
				      true,
				      ctx.get("SYS_MUO_ID", "nominal")));

  m_sf_iso.reset(new MCMuonScaleFactor(ctx,
				       ctx.get("SF_MUO_ISO_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2018/Muon_Iso_SF_RunABCD.root"),
				       ctx.get("SF_MUO_ISO_FILE", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"),
				       1.0,
				       "isolation",
				       true,
				       ctx.get("SYS_MUO_ISO", "nominal")));
}

bool MuonScaleFactors2018::process(Event &event) {

  if ( m_rng->Uniform() < m_lumi_fraction )
    {
      m_sf_trigger_before->process(event);
    }
  else
    {
      m_sf_trigger_after->process(event);
    }
  m_sf_id->process(event);
  m_sf_iso->process(event);

  return true;
}


ElectronScaleFactors2016::ElectronScaleFactors2016(Context &ctx) {

  m_sf_reco.reset(new MCElecScaleFactor(ctx,
				      ctx.get("SF_ELE_RECO_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2016/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root"),
				      1.0,
				      "reco",
				      ctx.get("SYS_ELE_RECO", "nominal")));

  m_sf_id.reset(new MCElecScaleFactor(ctx,
				      ctx.get("SF_ELE_ID_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2016/2016LegacyReReco_ElectronTight_Fall17V2.root"),
				      1.0,
				      "tight_id",
				      ctx.get("SYS_ELE_ID", "nominal")));


}

bool ElectronScaleFactors2016::process(Event &event) {

  m_sf_reco->process(event);
  m_sf_id->process(event);

  return true;
}

ElectronScaleFactors2017::ElectronScaleFactors2017(Context &ctx) {

  m_sf_reco.reset(new MCElecScaleFactor(ctx,
				      ctx.get("SF_ELE_RECO_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2017/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root"),
				      1.0,
				      "reco",
				      ctx.get("SYS_ELE_RECO", "nominal")));

  m_sf_id.reset(new MCElecScaleFactor(ctx,
				      ctx.get("SF_ELE_ID_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2017/2017_ElectronTight.root"),
				      1.0,
				      "tight_id",
				      ctx.get("SYS_ELE_ID", "nominal")));


}

bool ElectronScaleFactors2017::process(Event &event) {

  m_sf_reco->process(event);
  m_sf_id->process(event);

  return true;
}


ElectronScaleFactors2018::ElectronScaleFactors2018(Context &ctx) {

  m_sf_reco.reset(new MCElecScaleFactor(ctx,
				      ctx.get("SF_ELE_RECO_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2018/Electron_reco_SF_RunABCD.root"),
				      1.0,
				      "reco",
				      ctx.get("SYS_ELE_RECO", "nominal")));

  m_sf_id.reset(new MCElecScaleFactor(ctx,
				      ctx.get("SF_ELE_ID_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/common/data/2018/Electron_ID_tight_SF_RunABCD.root"),
				      1.0,
				      "tight_id",
				      ctx.get("SYS_ELE_ID", "nominal")));


}

bool ElectronScaleFactors2018::process(Event &event) {

  m_sf_reco->process(event);
  m_sf_id->process(event);

  return true;
}

LeptonScaleFactors::LeptonScaleFactors(Context &ctx) {
  m_sf_lepton.reset(new YearSwitcher(ctx));

  if (ctx.get("analysis_channel") == "ELECTRON")
    {
      m_sf_lepton->setup2016(std::make_shared<ElectronScaleFactors2016>(ctx));
      m_sf_lepton->setup2017(std::make_shared<ElectronScaleFactors2017>(ctx));
      m_sf_lepton->setup2018(std::make_shared<ElectronScaleFactors2018>(ctx));
    }
  else if (ctx.get("analysis_channel") == "MUON")
    {
      m_sf_lepton->setup2016(std::make_shared<MuonScaleFactors2016>(ctx));
      m_sf_lepton->setup2017(std::make_shared<MuonScaleFactors2017>(ctx));
      m_sf_lepton->setup2018(std::make_shared<MuonScaleFactors2018>(ctx));
 }
  else
    {
      throw invalid_argument("LeptonScaleFactorss::LeptonScaleFactors : no analysis channel specified in xml file (MUON or ELECTRON)");
    }

}

bool LeptonScaleFactors::process(Event &event) {
  return m_sf_lepton->process(event);
}
