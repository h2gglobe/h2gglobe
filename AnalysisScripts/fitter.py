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

### ROOT.gSystem.Load("libRooFit.so")
### ROOT.gSystem.Load("libPhysics.so");
### ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("../libLoopAll.so");
#ROOT.gSystem.Load("../../../../../lib/slc5_amd64_gcc434/libHiggsAnalysisHiggsToGammaGamma.so")

ROOT.gROOT.SetBatch(1)

ROOT.gBenchmark.Start("Analysis");

config_file="datafiles.dat"
if options.inputDat:
    config_file = options.inputDat

ut = ROOT.LoopAll();
type_run = 0
if options.typeRun != -1:
  type_run = options.typeRun

options.preSearchPath.reverse()
seach_path=options.preSearchPath+options.searchPath.split(":")+options.postSearchPath
cfg = configProducer(ut,config_file,type_run,int(options.nJobs),int(options.jobId),seach_path)

ROOT.gROOT.cd()
if not options.dryRun:
  ut.LoopAndFillHistos();
  ROOT.gBenchmark.Show("Analysis");

  if type_run == 0:
    ut.WriteFits();  

    jsonname = "%s.json" % ( ut.histFileName.rsplit(".root")[0] )

    if "/" in str(ut.histFileName):
      path,name = str(ut.histFileName).rsplit("/",1)
      ut.histFileName=path+"/histograms_"+name
    else:
      ut.histFileName="histograms_"+ut.histFileName
    ut.WriteHist();  
    
    
    # for small files, Keep local
    print "Producing JSON file for lumi calculation.."

    runLines = dumpLumi([ut.histFileName])

    print "Processed lumi sections: "
    print "{" + ", ".join(runLines) + "}"

    json = open( "%s" % jsonname, "w+" )
    print >>json, "{" + ", ".join(runLines) + "}"
    json.close()

    if "/" in str(ut.histFileName):
      ut.histFileName="histograms_"+name
    ut.WriteCounters();  

else:
  cfg.print_members()

ROOT.gROOT.Reset()
