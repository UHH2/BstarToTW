#include "UHH2/BstarToTW/include/BstarToTWSystematics.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/Utils.h"

using namespace std;
using namespace uhh2;

MuonScaleFactors2018::MuonScaleFactors2018(Context &ctx) {

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
				      ctx.get("SF_MUO_ID", "NUM_TightID_DEN_TrackerMuons_pt_abseta"),
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

  if (event.run < m_hlt_runnr)
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
