#!/usr/bin/env python

import ROOT
from sys import argv
from python.configProducer import *
from python.lumi import *
import sys

ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("../libLoopAll.so");

ROOT.gROOT.SetBatch(1)

ROOT.gBenchmark.Start("Analysis");

config_file="datafiles.dat"
if len(sys.argv) > 1:
    config_file = sys.argv[1]

njobs=-1
jobId=0
if len(argv) > 2:
    if len(argv) == 3:
        print "Usage: fitter.py <config_file> [<njobs> <job_id>]" 
        exit(1)
    njobs = int(argv[2])
    jobId = int(argv[3])

ut = ROOT.LoopAll();
cfg = configProducer(ut,config_file,0,njobs,jobId)
  
ut.LoopAndFillHistos();
ROOT.gBenchmark.Show("Analysis");

ut.WriteFits();  

jsonname = "%s.json" % ( ut.histFileName.rsplit(".root")[0] )

if "/" in str(ut.histFileName):
    path,name = str(ut.histFileName).rsplit("/",1)
    ut.histFileName=path+"/histograms_"+name
else:
    ut.histFileName="histograms_"+ut.histFileName
ut.WriteHist();  

print "Producing JSON file for lumi calculation.."
runLines = dumpLumi([ut.histFileName])

print "Processed lumi sections: "
print "{" + ", ".join(runLines) + "}"

json = open( "%s" % jsonname, "w+" )
print >>json, "{" + ", ".join(runLines) + "}"
json.close()

ut.WriteCounters();  

ROOT.gROOT.Reset()
