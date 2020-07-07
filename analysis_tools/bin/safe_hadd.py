#!/usr/bin/env python

import argparse
from glob import glob
import ROOT
import subprocess

def check_tree(filename):
    """
    check if AnalysisTree is in root file and has entries
    """
    rootfile = ROOT.TFile.Open(filename)
    try:
        entries = rootfile.AnalysisTree.GetEntriesFast()
        print filename,' entries: ',entries,entries>0
        return (entries > 0)
    except:
        print 'AnalysisTree in file '+filename+' has 0 entries: '
        return False

def spawn_hadd(file_out_name, filenames, nice='10'):
    """
    spawn hadd using correct order of files to avoid empty AnalysisTree
    """
    command = ['nice', '-n', nice, 'hadd', '-f', file_out_name] + filenames
    proc = subprocess.Popen(command)
    proc.wait()

if __name__ == "__main__":
    # parse args
    parser = argparse.ArgumentParser(description='Check files for valid AnalysisTree and change hadd order to avoid empty AnalysisTree of hadded root file.')
    parser.add_argument('outputname', help='path to output file. If this exists, it will be overwritten')
    parser.add_argument('filenames', help='path to the files to hadd.', nargs='+')
    parser.add_argument('--nice', '-n', default=10, type=int, help='level of niceness')
    args = parser.parse_args()

    filenames = args.filenames
    file_out_name = args.outputname
    for i in range(0, len(filenames)):
        if (check_tree(filenames[i])):
            # swap elements
            filenames[i], filenames[0] = filenames[0], filenames[i]
            break
    
    spawn_hadd(file_out_name, filenames, str(args.nice))
