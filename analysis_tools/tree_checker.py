#!/usr/bin/env python

import ROOT
import sys
from glob import glob

def check_TreeExists(filename,treename):
     rootfile = ROOT.TFile.Open(filename)
     if treename!='AnalysisTree':
          raise ValueError('treename is not \'AnalysisTree\'')
     #print filename
     try:
          #rootTree = rootfile.Get(treename)
          entries = rootfile.AnalysisTree.GetEntriesFast()
          print filename,' entries: ',entries,entries>0
          return  entries>0
     except:
          #print 'False'
          print 'tree '+treename+' in file '+filename+' has ==0 entries: '
          return False
          #entries = rootTree.GetEntriesFast()
          #if rootTree: return False
          #return True


if __name__ == "__main__":
     for arg in sys.argv[2:]:
          if '*' in arg:
               for itfile in glob(pattern):
                     check_TreeExists(itfile,sys.argv[1])
          else:
               check_TreeExists(arg,sys.argv[1])
