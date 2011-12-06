#!/usr/bin/env python

import ROOT
from sys import argv
from python.configProducer import *
from python.lumi import *
import sys
from python.runOptions import *
(options,args)=parser.parse_args()
if (int(options.nJobs) > 0) and (int(options.jobId) >= int(options.nJobs)):
  sys.exit("Job id's must run from 0 -> %d when splitting into %d jobs"%(int(options.nJobs)-1,int(options.nJobs)))

ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("../libLoopAll.so");
#ROOT.gSystem.Load("../../../../../lib/slc5_amd64_gcc434/libHiggsAnalysisHiggsToGammaGamma.so")

ROOT.gROOT.SetBatch(1)

ROOT.gBenchmark.Start("Analysis");

config_file="datafiles.dat"
if options.inputDat:
    config_file = options.inputDat

ut = ROOT.LoopAll();
cfg = configProducer(ut,config_file,0,int(options.nJobs),int(options.jobId))

if not options.dryRun:
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
