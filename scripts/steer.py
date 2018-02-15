import subprocess
import os
import time

from ModuleRunner import *
from MacroRunner import *
from constants import *


"""This macro and guide describes and steers the complete progress of the LQ->t mu search."""


"""
First of all, the user has to set up the Preselection for the Signal Region. For this, basically two files have to be modified.
1) LQToTopMu/src/LQToTopMuPreselectionModule.cxx
2) LQToTopMu/config/LQToTopMuPreselection.xml
"""

#Shortcuts to options
do_JECJER = True

#Modify these items, everything else will work. The files given here must already be adapted to new settings in case there are any.
#path_LQDIR = '/nfs/dust/cms/user/reimersa/CMSSWTEST/CMSSW_8_0_24_patch1/src/UHH2/LQToTopMu'
path_LQDIR = '/nfs/dust/cms/user/reimersa/CMSSW_8_0_24_patch1/src/UHH2/LQToTopMu'
path_TagProbeMacros = '/nfs/dust/cms/user/reimersa/LQToTopMu/Run2_80X_v3/CompleteAnalysis/macros/TagProbe/src'
name_PRESEL_TAGPROBE = 'LQToTopMuTagProbePreselection.xml'
name_PRESEL_FAKERATE = 'LQToTopMuFakeRatePreselection.xml'
name_PRESEL_FAKERATE_MU = 'LQToTopMuFakeRatePreselection_Muon.xml'
name_PRESEL_SR = 'LQToTopMuPreselection.xml'
name_PRESEL_CR = 'LQToTopMuSidebandPreselection_EE.xml'
name_SEL_TAGPROBE = 'LQToTopMuTagProbeAnalysis.xml'
name_SEL_FAKERATE_DIBOSON = 'LQToTopMuFakeRateAnalysis_Diboson.xml'
module_SEL_FAKERATE_DIBOSON = 'LQToTopMuFakeRateAnalysisModule.cxx'
name_SEL_FAKERATE_ELE = 'LQToTopMuFakeRateAnalysis_Electron.xml'
name_SEL_FAKERATE_MU = 'LQToTopMuFakeRateAnalysis_Muon.xml'
name_SEL_SR = 'LQToTopMuAnalysis.xml'
name_SEL_CR = 'LQToTopMuSidebandAnalysis_EE.xml'
closure_postfix = 'Closure'
#systematics = ['NOMINAL','PU', 'ALPHA', 'BTAG_bc', 'BTAG_udsg', 'ELEID', 'ELEREC', 'ELETR', 'MUID', 'MUISO', 'MUTR', 'PDF', 'SCALE', 'JEC', 'JER', 'RATE_ttbar', 'RATE_dy', 'RATE_st', 'RATE_wj', 'RATE_qcd', 'RATE_db', 'RATE_ttv', 'MUTRK', 'ELEFAKE', 'MUFAKE']
#systematics = ['NOMINAL','PU', 'ALPHA', 'BTAG_bc', 'BTAG_udsg', 'ELEID', 'ELEREC', 'ELETR', 'MUID', 'MUISO', 'MUTR', 'PDF', 'SCALE_ttbar', 'SCALE_dy', 'SCALE_st', 'SCALE_wj', 'SCALE_db', 'SCALE_ttv', 'JEC', 'JER', 'RATE_ttbar', 'RATE_dy', 'RATE_st', 'RATE_wj', 'RATE_qcd', 'RATE_db', 'RATE_ttv', 'MUTRK', 'ELEFAKE', 'MUFAKE']
systematics = ['SCALE_ttbar', 'SCALE_dy', 'SCALE_st', 'SCALE_wj', 'SCALE_db', 'SCALE_ttv']
#systematics = ['SCALE_db']
systematics_pre = ['NOMINAL', 'JEC', 'JER']
ModuleRunner = ModuleRunner(path_LQDIR)
TagProbeRunner = MacroRunner(path_TagProbeMacros)



###      PRESEL FOR TagAndProbe
###      ======================

#ModuleRunner.CompileModules()
#ModuleRunner.BuildBatch(name_PRESEL_TAGPROBE)
#ModuleRunner.SubmitBatch(name_PRESEL_TAGPROBE)
#ModuleRunner.RetrieveLogBatch(name_PRESEL_TAGPROBE, kill_proc=False)
#ModuleRunner.SubmitMissingFiles(name_PRESEL_TAGPROBE)
#ModuleRunner.ForceAddBatch(name_PRESEL_TAGPROBE)
#while isThisRunning('hadd'): time.sleep(20)
#ModuleRunner.AddFiles(name_PRESEL_TAGPROBE, True)
#ModuleRunner.CleanUp(name_PRESEL_TAGPROBE)


###      FULL SEL FOR TagAndProbe
###      ========================

#ModuleRunner.CompileModules()
#ModuleRunner.ModifyEntity(name_SEL_TAGPROBE, 'CLOSURETEST', 'true', closure_postfix)
#ModuleRunner.MakeOutDir(name_SEL_TAGPROBE)
#ModuleRunner.RunLocal(name_SEL_TAGPROBE)
#ModuleRunner.AddFiles(name_SEL_TAGPROBE, True)

### Calculate TriggerSF with macros
#TagProbeRunner.SetInputDir(out_dir(ModuleRunner.path, name_SEL_TAGPROBE), 'TriggerTool')
#TagProbeRunner.CompileMacros()
#TagProbeRunner.ExecuteMacros()



###      PRESEL FOR FAKERATE AND DIBOSON, ELE & MUON
###      ===========================================
#ModuleRunner.CompileModules()
#ModuleRunner.BuildBatch(name_PRESEL_FAKERATE)
#ModuleRunner.SubmitBatch(name_PRESEL_FAKERATE)
#ModuleRunner.RetrieveLogBatch(name_PRESEL_FAKERATE, kill_proc=False)
#ModuleRunner.SubmitMissingFiles(name_PRESEL_FAKERATE)
#ModuleRunner.ForceAddBatch(name_PRESEL_FAKERATE)
#while isThisRunning('hadd'): time.sleep(20)
#ModuleRunner.AddFiles(name_PRESEL_FAKERATE, False)
#ModuleRunner.CleanUp(name_PRESEL_FAKERATE)


#ModuleRunner.CompileModules()
#ModuleRunner.BuildBatch(name_PRESEL_FAKERATE_MU)
#ModuleRunner.SubmitBatch(name_PRESEL_FAKERATE_MU)
#ModuleRunner.RetrieveLogBatch(name_PRESEL_FAKERATE_MU, kill_proc=False)
#ModuleRunner.SubmitMissingFiles(name_PRESEL_FAKERATE_MU)
#ModuleRunner.ForceAddBatch(name_PRESEL_FAKERATE_MU)
#while isThisRunning('hadd'): time.sleep(20)
#ModuleRunner.AddFiles(name_PRESEL_FAKERATE_MU, False)
#ModuleRunner.CleanUp(name_PRESEL_FAKERATE_MU)


###      FULL SEL FOR DIBOSON
###      ====================

#1st time --> For BTag SF
#ModuleRunner.PrepareBTagSF(module_SEL_FAKERATE_DIBOSON, True)
#ModuleRunner.CompileModules()




###      FULL SEL FOR ELECTRON
###      =====================
#currently only a dummy, create a proper cycle for this!
#ModuleRunner.SubmitMissingFiles(name_SEL_FAKERATE_ELE)
#ModuleRunner.AddFiles(name_SEL_FAKERATE_ELE, False)



###      PRESEL for SR
###      =============

#ModuleRunner.CompileModules()
#ModuleRunner.create_JECJER_xmlfiles(name_PRESEL_SR)
"""
if do_JECJER:
    for sys in systematics_pre:
        if sys != 'NOMINAL':
            for direction in ['up', 'down']:
                ModuleRunner.BuildBatch(name_PRESEL_SR, sys, direction)
                ModuleRunner.SubmitBatch(name_PRESEL_SR, sys, direction)
                ModuleRunner.RetrieveLogBatch(name_PRESEL_SR, sys, direction, kill_proc=False)
                ModuleRunner.ResubmitBatch(name_PRESEL_SR, sys, direction)
                ModuleRunner.RetrieveLogBatch(name_PRESEL_SR, sys, direction, kill_proc=False)
                ModuleRunner.SubmitMissingFiles(name_PRESEL_SR, sys, direction)
                ModuleRunner.RetrieveLogBatch(name_PRESEL_SR, sys, direction, kill_proc=True)
                ModuleRunner.AddBatch(name_PRESEL_SR, sys, direction)
                while isThisRunning('hadd'): time.sleep(20)
                ModuleRunner.CleanUp(name_PRESEL_SR, sys, direction)
        else: 
            ModuleRunner.BuildBatch(name_PRESEL_SR)
            ModuleRunner.SubmitBatch(name_PRESEL_SR)
            ModuleRunner.RetrieveLogBatch(name_PRESEL_SR, kill_proc=False)
            ModuleRunner.ResubmitBatch(name_PRESEL_SR)
            ModuleRunner.RetrieveLogBatch(name_PRESEL_SR, kill_proc=False)
            ModuleRunner.SubmitMissingFiles(name_PRESEL_SR)
            ModuleRunner.RetrieveLogBatch(name_PRESEL_SR, kill_proc=True)
            ModuleRunner.AddBatch(name_PRESEL_SR)
            while isThisRunning('hadd'): time.sleep(20)
            ModuleRunner.CleanUp(name_PRESEL_SR)
else: 
    ModuleRunner.BuildBatch(name_PRESEL_SR)
    ModuleRunner.SubmitBatch(name_PRESEL_SR)
    ModuleRunner.RetrieveLogBatch(name_PRESEL_SR, kill_proc=False)
    ModuleRunner.ResubmitBatch(name_PRESEL_SR)
    ModuleRunner.RetrieveLogBatch(name_PRESEL_SR, kill_proc=False)
    ModuleRunner.SubmitMissingFiles(name_PRESEL_SR)
    ModuleRunner.RetrieveLogBatch(name_PRESEL_SR, kill_proc=True)
    ModuleRunner.AddBatch(name_PRESEL_SR)
    while isThisRunning('hadd'): time.sleep(20)
    ModuleRunner.CleanUp(name_PRESEL_SR)
"""


###      PRESEL FOR CR
###      =============

#ModuleRunner.CompileModules()
#ModuleRunner.create_JECJER_xmlfiles(name_PRESEL_CR)
"""
if do_JECJER:
    for sys in systematics_pre:
        if sys != 'NOMINAL':
            for direction in ['up', 'down']:
                ModuleRunner.BuildBatch(name_PRESEL_CR, sys, direction)
                ModuleRunner.SubmitBatch(name_PRESEL_CR, sys, direction)
                ModuleRunner.RetrieveLogBatch(name_PRESEL_CR, sys, direction, kill_proc=False)
                ModuleRunner.ResubmitBatch(name_PRESEL_CR, sys, direction)
                ModuleRunner.RetrieveLogBatch(name_PRESEL_CR, sys, direction, kill_proc=False)
                ModuleRunner.SubmitMissingFiles(name_PRESEL_CR, sys, direction)
                ModuleRunner.RetrieveLogBatch(name_PRESEL_CR, sys, direction, kill_proc=True)
                ModuleRunner.AddBatch(name_PRESEL_CR, sys, direction)
                while isThisRunning('hadd'): time.sleep(20)
        else: 
            ModuleRunner.BuildBatch(name_PRESEL_CR)
            ModuleRunner.SubmitBatch(name_PRESEL_CR)
            ModuleRunner.RetrieveLogBatch(name_PRESEL_CR, kill_proc=False)
            ModuleRunner.ResubmitBatch(name_PRESEL_CR, sys, direction)
            ModuleRunner.RetrieveLogBatch(name_PRESEL_CR, sys, direction, kill_proc=False)
            ModuleRunner.SubmitMissingFiles(name_PRESEL_CR)
            ModuleRunner.RetrieveLogBatch(name_PRESEL_CR, kill_proc=True)
            ModuleRunner.AddBatch(name_PRESEL_CR)
            while isThisRunning('hadd'): time.sleep(20)
    #ModuleRunner.AddFilesWithSyst(name_PRESEL_CR, False, systematics_pre)
    ModuleRunner.CleanUp(name_PRESEL_CR, sys, direction)
else: 
    ModuleRunner.BuildBatch(name_PRESEL_CR)
    ModuleRunner.SubmitBatch(name_PRESEL_CR)
    ModuleRunner.RetrieveLogBatch(name_PRESEL_CR, kill_proc=False)
    ModuleRunner.ResubmitBatch(name_PRESEL_CR)
    ModuleRunner.RetrieveLogBatch(name_PRESEL_CR, kill_proc=False)
    ModuleRunner.SubmitMissingFiles(name_PRESEL_CR)
    ModuleRunner.RetrieveLogBatch(name_PRESEL_CR, kill_proc=True)
    ModuleRunner.AddBatch(name_PRESEL_CR)
    while isThisRunning('hadd'): time.sleep(20)
    #ModuleRunner.AddFiles(name_PRESEL_CR, False)
    ModuleRunner.CleanUp(name_PRESEL_CR)
"""
       



### FULL SEL FOR SR

#ModuleRunner.CompileModules()
#configfiles = ModuleRunner.CreateSystXMLFiles(name_SEL_SR, systematics, systematics_dict_SR, 'AllSystsSR')
#ModuleRunner.RunSystJobsLocal(configfiles)
#ModuleRunner.AddFilesWithSyst(name_SEL_SR, True, systematics)
#ModuleRunner.CreateSCALEDirs(name_SEL_SR, systematics)


###FULL SEL FOR CR

#ModuleRunner.CompileModules()
#configfiles = ModuleRunner.CreateSystXMLFiles(name_SEL_CR, systematics, systematics_dict_CR, 'AllSystsCR')
#ModuleRunner.RunSystJobsLocal(configfiles)
#ModuleRunner.AddFilesWithSyst(name_SEL_CR, False, systematics)
#ModuleRunner.CreateSCALEDirs(name_SEL_CR, systematics)


###MACROS FOR CR


### CR 2nd TIME


###REMAINING MACROS


###SYSTEMATICS SR


###SYSTEMATICS CR, 1st TIME


###MACROS CR


###SYSTEMATICS CR, 2nd TIME


###REMAINING MACROS


###LIMITS

