#!/bin/env python

import ROOT
import array
from optparse import OptionParser
from pprint import pprint
import os.path as p
import os, sys
from subprocess import check_call as call

runBinning = [0, 197495, 203767, 210000]
h2 = ROOT.TH2F("pu_2D", "pu_2D", 100, 0, 100, 3, array.array('d',runBinning))

parser = OptionParser(usage="usage: %prog [options] EOS_source_directory\nrun with --help to get list of options")
parser.add_option("--putBack",    action="store_true", default=False, help="Put back merged file in source directory [default: %default].")
parser.add_option("--putHadoop",  action="store_true", default=False, help="Put back merged file in source directory [default: %default].")

(options, args) = parser.parse_args()

if len(args) != 1:
    parser.print_usage()
    raise RuntimeError("Need exactly one directory to work with.")

options.inDir = args[0]
options.inDirName = p.basename(options.inDir)
options.inDirBase = p.dirname(options.inDir)
options.eosLs = "/afs/cern.ch/project/eos/installation/pro/bin/eos root://eoscms ls /eos/cms"

if not options.inDirName:
    raise RuntimeError("Empty target directory name (which defines the sample name). Check path.")

try:
    ld_path = os.getenv("LD_LIBRARY_PATH")
except:
    ld_path = ""
if not "/afs/cern.ch/project/eos/installation/pro/lib64/" in ld_path:
    os.putenv("LD_LIBRARY_PATH", "%s:%s" % ( ld_path, "/afs/cern.ch/project/eos/installation/pro/lib64/" ) )
    
if options.putHadoop:
    call( """ls %(inDir)s | awk '/root$/ { print \"%(inDir)s/\"$1  }' | sed 's/\?.*$//' > %(inDirName)s.files.txt"""
      % vars(options), shell=True)
else: 
    call( """%(eosLs)s%(inDir)s | awk '/root$/ { print \"root://eoscms//eos/cms%(inDir)s/\"$1  }' | sed 's/\?.*$//' > %(inDirName)s.files.txt"""
      % vars(options), shell=True)
call( """rm -f %(inDirName)s.pileup.root %(inDirName)s.pileup.root.log""" % vars(options), shell=True)

file = open("%(inDirName)s.files.txt"%vars(options))
lines = file.readlines()
file.close()

for l in lines:
    rootFile = ROOT.TFile.Open(l.split("\n")[0])
    tree = rootFile.Get("event")
    tree.SetBranchStatus("*", 0)
    tree.SetBranchStatus("pu_n", 1)
    tree.SetBranchStatus("run", 1)

    entries = tree.GetEntries()
    for z in xrange(entries):
        tree.GetEntry(z)
        h2.Fill(tree.pu_n, tree.run)
        h1.Fill(tree.run)
        
out = ROOT.TFile("%(inDirName)s.pileup.root"%vars(options), "recreate")
h2.Write()
out.Close()

if options.putHadoop:
    cmd = """cp %(inDirName)s.pileup.root %(inDirBase)s/%(inDirName)s.pileup.root"""
    print "Copying back root"
    call(cmd % vars(options), shell=True)
    
if options.putBack:
    cmd = """cmsStage -f %(inDirName)s.pileup.root %(inDirBase)s/%(inDirName)s.pileup.root"""
    print "Copying back root"
    call(cmd % vars(options), shell=True)

