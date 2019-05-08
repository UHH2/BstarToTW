#include "UHH2/BstarToTW/include/BstarToTWSystematics.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/Utils.h"

using namespace std;
using namespace uhh2;

MuonScaleFactors2018::MuonScaleFactors2018(Context &ctx) {

  m_sf_trigger_before.reset(new MCMuonScaleFactor(ctx,
						  ctx.get("SF_MUO_TRIGGER_BEFORE_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/BstarToTW/data/EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root"),
						  ctx.get("SF_MUO_TRIGGER_BEFORE_FILE", "IsoMu24_PtEtaBins"),
						  0.5,
						  "trigger",
						  false,
						  ctx.get("SYS_MUO_TRIGGER", "nominal")));

  m_sf_trigger_after.reset(new MCMuonScaleFactor(ctx,
						 ctx.get("SF_MUO_TRIGGER_AFTER_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/BstarToTW/data/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root"),
						 ctx.get("SF_MUO_TRIGGER_AFTER_FILE", "IsoMu24_PtEtaBins"),
						 0.5,
						 "trigger",
						 false,
						 ctx.get("SYS_MUO_TRIGGER", "nominal")));

  m_sf_id.reset(new MCMuonScaleFactor(ctx,
				      ctx.get("SF_MUO_ID_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/BstarToTW/data/RunABCD_SF_ID.root"),
				      ctx.get("SF_MUO_ID", "NUM_TightID_DEN_TrackerMuons_pt_abseta"),
				      1.0,
				      "tight_id",
				      true,
				      ctx.get("SYS_MUO_ID", "nominal")));

  m_sf_iso.reset(new MCMuonScaleFactor(ctx,
				       ctx.get("SF_MUO_ISO_PATH", "/nfs/dust/cms/user/froehlia/CMSSW_10_2_10/src/UHH2/BstarToTW/data/RunABCD_SF_ISO.root"),
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

