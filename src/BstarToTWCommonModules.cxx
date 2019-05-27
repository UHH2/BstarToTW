#include "UHH2/BstarToTW/include/BstarToTWCommonModules.h"

using namespace std;
using namespace uhh2;

JEC2016Module::JEC2016Module(Context &ctx) {
  is_mc = ctx.get("dataset_type") == "MC";
  // JetCorrectors
  if(is_mc)
    {  
      jec_ak4_mc.reset(new JetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_L123_AK4PFchs_MC));
      // jet_resolution_smearer.reset(new JetResolutionSmearer(ctx));
      jlc_ak4_mc.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Summer16_07Aug2017_V11_L123_AK4PFchs_MC));
      jec_hotvr_mc.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_L123_AK4PFPuppi_MC));
    }  
  else
    { 
      // -- ak4
      // JEC
      jec_ak4_B.reset(new JetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_B_L123_AK4PFchs_DATA));
      jec_ak4_C.reset(new JetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_C_L123_AK4PFchs_DATA));
      jec_ak4_D.reset(new JetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_D_L123_AK4PFchs_DATA));
      jec_ak4_E.reset(new JetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_E_L123_AK4PFchs_DATA));
      jec_ak4_F.reset(new JetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_F_L123_AK4PFchs_DATA));
      jec_ak4_G.reset(new JetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_G_L123_AK4PFchs_DATA));
      jec_ak4_H.reset(new JetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_H_L123_AK4PFchs_DATA));
      // Jet Lepton Cleaning
      jlc_ak4_B.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Summer16_07Aug2017_V11_B_L123_AK4PFchs_DATA));
      jlc_ak4_C.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Summer16_07Aug2017_V11_C_L123_AK4PFchs_DATA));
      jlc_ak4_D.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Summer16_07Aug2017_V11_D_L123_AK4PFchs_DATA));
      jlc_ak4_E.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Summer16_07Aug2017_V11_E_L123_AK4PFchs_DATA));
      jlc_ak4_F.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Summer16_07Aug2017_V11_F_L123_AK4PFchs_DATA));
      jlc_ak4_G.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Summer16_07Aug2017_V11_G_L123_AK4PFchs_DATA));
      jlc_ak4_H.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Summer16_07Aug2017_V11_H_L123_AK4PFchs_DATA));

      // -- HOTVR
      // JEC
      jec_hotvr_B.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_B_L123_AK4PFPuppi_DATA));
      jec_hotvr_C.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_C_L123_AK4PFPuppi_DATA));
      jec_hotvr_D.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_D_L123_AK4PFPuppi_DATA));
      jec_hotvr_E.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_E_L123_AK4PFPuppi_DATA));
      jec_hotvr_F.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_F_L123_AK4PFPuppi_DATA));
      jec_hotvr_G.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_G_L123_AK4PFPuppi_DATA));
      jec_hotvr_H.reset(new HOTVRJetCorrector(ctx, JERFiles::Summer16_07Aug2017_V11_H_L123_AK4PFPuppi_DATA));
    }  
  jlc_hotvr.reset(new HOTVRJetLeptonCleaner());

}

bool JEC2016Module::process(Event &event) {

  if (b_jetlep)
    jlc_hotvr->process(event);
  if(is_mc) 
    {
      if (b_jetlep)
	jlc_ak4_mc->process(event);
      jec_ak4_mc->process(event);
      jec_ak4_mc->correct_met(event);
      // jet_resolution_smearer->process(event);
      jec_hotvr_mc->process(event);
    }
  else{
	if(event.run <= runnr_B)
	  {
	    if (b_jetlep)
	      jlc_ak4_B->process(event);
	    jec_ak4_B->process(event);
	    jec_hotvr_B->process(event);
	  }
	else if(event.run <= runnr_C)
	  {
	    if (b_jetlep)
	      jlc_ak4_C->process(event);
	    jec_ak4_C->process(event);
	    jec_hotvr_C->process(event);
	  }

	else if(event.run <= runnr_D)
	  {
	    if (b_jetlep)
	      jlc_ak4_D->process(event);
	    jec_ak4_D->process(event);
	    jec_hotvr_D->process(event);
	  }

	else if(event.run <= runnr_E)
	  {
	    if (b_jetlep)
	      jlc_ak4_E->process(event);
	    jec_ak4_E->process(event);
	    jec_hotvr_E->process(event);
	  }

	else if(event.run < runnr_F) // "<" is correct since checking for Fearly
	  {
	    if (b_jetlep)
	      jlc_ak4_F->process(event);
	    jec_ak4_F->process(event);
	    jec_hotvr_F->process(event);
	  }

	else if(event.run <= runnr_G)
	  {
	    if (b_jetlep)
	      jlc_ak4_G->process(event);
	    jec_ak4_G->process(event);
	    jec_hotvr_G->process(event);
	  }
	else if(event.run <= runnr_H)
	  {
	    if (b_jetlep)
	      jlc_ak4_H->process(event);
	    jec_ak4_H->process(event);
	    jec_hotvr_H->process(event);
	  }
	else runtime_error("ObjectCleaner.cxx: run number not covered by if-statements in process-routine.");
  }
  return true;
}

JEC2018Module::JEC2018Module(Context &ctx) {
  is_mc = ctx.get("dataset_type") == "MC";
  // JetCorrectors
  if(is_mc)
    {  
      jec_ak4_mc.reset(new JetCorrector(ctx, JERFiles::Autumn18_V8_L123_AK4PFchs_MC));
      // jet_resolution_smearer.reset(new JetResolutionSmearer(ctx));
      jlc_ak4_mc.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Autumn18_V8_L123_AK4PFchs_MC));
      jec_hotvr_mc.reset(new HOTVRJetCorrector(ctx, JERFiles::Autumn18_V8_L123_AK4PFPuppi_MC));
    }  
  else
    { 
      // -- ak4
      // JEC
      jec_ak4_A.reset(new JetCorrector(ctx, JERFiles::Autumn18_V8_A_L123_AK4PFchs_DATA));
      jec_ak4_B.reset(new JetCorrector(ctx, JERFiles::Autumn18_V8_B_L123_AK4PFchs_DATA));
      jec_ak4_C.reset(new JetCorrector(ctx, JERFiles::Autumn18_V8_C_L123_AK4PFchs_DATA));
      jec_ak4_D.reset(new JetCorrector(ctx, JERFiles::Autumn18_V8_D_L123_AK4PFchs_DATA));
      // Jet Lepton Cleaning
      jlc_ak4_A.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Autumn18_V8_A_L123_AK4PFchs_DATA));
      jlc_ak4_B.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Autumn18_V8_B_L123_AK4PFchs_DATA));
      jlc_ak4_C.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Autumn18_V8_C_L123_AK4PFchs_DATA));
      jlc_ak4_D.reset(new JetLeptonCleaner_by_KEYmatching(ctx, JERFiles::Autumn18_V8_D_L123_AK4PFchs_DATA));

      // -- HOTVR
      // JEC
      jec_hotvr_A.reset(new HOTVRJetCorrector(ctx, JERFiles::Autumn18_V8_A_L123_AK4PFPuppi_DATA));
      jec_hotvr_B.reset(new HOTVRJetCorrector(ctx, JERFiles::Autumn18_V8_B_L123_AK4PFPuppi_DATA));
      jec_hotvr_C.reset(new HOTVRJetCorrector(ctx, JERFiles::Autumn18_V8_C_L123_AK4PFPuppi_DATA));
      jec_hotvr_D.reset(new HOTVRJetCorrector(ctx, JERFiles::Autumn18_V8_D_L123_AK4PFPuppi_DATA));
    }  
  jlc_hotvr.reset(new HOTVRJetLeptonCleaner());

}

bool JEC2018Module::process(Event &event) {

  if (b_jetlep)
    jlc_hotvr->process(event);
  if(is_mc) 
    {
      if (b_jetlep)
	jlc_ak4_mc->process(event);
      jec_ak4_mc->process(event);
      jec_ak4_mc->correct_met(event);
      // jet_resolution_smearer->process(event);
      jec_hotvr_mc->process(event);
    }
  else{
	if(event.run <= runnr_A)
	  {
	    if (b_jetlep)
	      jlc_ak4_A->process(event);
	    jec_ak4_A->process(event);
	    jec_ak4_A->correct_met(event);
	    jec_hotvr_A->process(event);
	  }
	else if(event.run <= runnr_B)
	  {
	    if (b_jetlep)
	      jlc_ak4_B->process(event);
	    jec_ak4_B->process(event);
	    jec_ak4_B->correct_met(event);
	    jec_hotvr_B->process(event);
	  }

	else if(event.run <= runnr_C)
	  {
	    if (b_jetlep)
	      jlc_ak4_C->process(event);
	    jec_ak4_C->process(event);
	    jec_ak4_C->correct_met(event);
	    jec_hotvr_C->process(event);
	  }

	else if(event.run <= runnr_D)
	  {
	    if (b_jetlep)
	      jlc_ak4_D->process(event);
	    jec_ak4_D->process(event);
	    jec_ak4_D->correct_met(event);
	    jec_hotvr_D->process(event);
	  }

	else runtime_error("ObjectCleaner.cxx: run number not covered by if-statements in process-routine.");
  }
  return true;
}


ObjectSetup::ObjectSetup(Context & ctx){

  is_mc = ctx.get("dataset_type") == "MC";
  year = extract_year(ctx);
  PrimaryVertexId id_pv = StandardPrimaryVertexId();
  MuonId id_muo_veto = AndId<Muon>(MuonID(Muon::Selector::CutBasedIdLoose), PtEtaCut(lepveto_pt_min, lep_eta_max), MuonIso(muo_iso_max)); // muon veto ID
  ElectronId id_ele_veto; // to be set era dependent
  JetId id_jetpfid = JetPFID(JetPFID::WP_TIGHT_CHS);
  JetId id_jet = PtEtaCut(jet_pt_min, jet_eta_max); // jet ID
  TopJetId id_topjet =  PtEtaCut(top_pt_min, top_eta_max);

  // Set era dependend stuff (JEC, IDs, whatnot)  
  if (year == Year::is2018) {
    id_ele_veto = AndId<Electron>(ElectronID_Fall17_veto, PtEtaCut(lepveto_pt_min, lep_eta_max)); // electron veto ID
    jec_module.reset(new JEC2018Module(ctx));
  }
  // else if (year == Year::is2017v2) jec_module.reset(new JEC2017Module(ctx));
  else if (year == Year::is2016v3) {
    id_ele_veto = AndId<Electron>(ElectronID_Summer16_veto, PtEtaCut(lepveto_pt_min, lep_eta_max)); // electron veto ID
    jec_module.reset(new JEC2016Module(ctx));
  }

  // Cleaner
  cl_pv.reset(new PrimaryVertexCleaner(id_pv));
  cl_muo.reset(new MuonCleaner(id_muo_veto));
  cl_ele.reset(new ElectronCleaner(id_ele_veto));
  cl_jetpfid.reset(new JetCleaner(ctx, id_jetpfid));
  cl_jet.reset(new JetCleaner(ctx, id_jet));
  cl_topjet.reset(new TopJetCleaner(ctx, id_topjet));

  // Metfilters
  metfilters_selection.reset(new AndSelection(ctx, "metfilters"));
  metfilters_selection->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
  metfilters_selection->add<TriggerSelection>("globalSuperTightHalo2016Filter", "Flag_globalSuperTightHalo2016Filter");
  metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
  if (!is_mc) metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");
  // metfilters_selection->add<TriggerSelection>("BadChargedCandidateFilter", "Flag_BadChargedCandidateFilter"); // Not recommended, under review. Separate module in ntuple_generator for 2016v2
  if (year != Year::is2016v2) {
    metfilters_selection->add<TriggerSelection>("BadPFMuonFilter", "Flag_BadPFMuonFilter");
  } else {
    metfilters_selection->add<TriggerSelection>("BadPFMuonFilter", "Extra_BadPFMuonFilter");
  }
  metfilters_selection->add<TriggerSelection>("goodVertices", "Flag_goodVertices");
  metfilters_selection->add<EcalBadCalibSelection>("EcalBadCalibSelection"); // Use this instead of Flag_ecalBadCalibFilter, uses ecalBadCalibReducedMINIAODFilter in ntuple_generator
  metfilters_selection->add<NPVSelection>("1 good PV",1,-1,id_pv);

}

bool ObjectSetup::process(Event & event){
  
  cl_pv->process(event);
  if(!metfilters_selection->passes(event)) return false;
  // lepton cleaning
  if (b_lepclean)
    {
      cl_muo->process(event);
      cl_ele->process(event);
    }

  // jet energy corrections and cleaning
  cl_jetpfid->process(event);
  jec_module->process(event);
  
  cl_jet->process(event);
  cl_topjet->process(event);

  return true;
}
