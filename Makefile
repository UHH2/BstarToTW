LIBRARY := SUHH2BstarToTW
USERLDFLAGS := -lSUHH2core -lSUHH2common -lGenVector -Wl,-rpath,/nfs/dust/cms/user/froehlia/CMSSW_8_0_20/src/FastJet/lib -lm -L/nfs/dust/cms/user/froehlia/CMSSW_8_0_20/src/FastJet/lib -lfastjettools -lfastjet -lHOTVR -lNsubjettiness -lRecursiveTools
USERCXXFLAGS := -I/nfs/dust/cms/user/froehlia/CMSSW_8_0_20/src/FastJet/include

# enable par creation; this is necessary for all packages containing AnalysisModules
# to be loaded from by AnalysisModuleRunner.
PAR := 1
include ../Makefile.common
