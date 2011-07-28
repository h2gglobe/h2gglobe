#!/usr/bin/env python

import ROOT
from python.configProducer import *
from python.runOptions import *
(options,args)=parser.parse_args()

ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("libRooFit.so");
ROOT.gSystem.Load("../libLoopAll.so");

ROOT.gBenchmark.Start("Analysis");

ut = ROOT.LoopAll();
config_file="inputfiles.dat"
if options.inputDat:
    config_file = options.inputDat
cfg = configProducer(ut,config_file,0,int(options.nJobs),int(options.jobId))
  
if not options.dryRun:
 ut.LoopAndFillHistos();
 ROOT.gBenchmark.Show("Analysis");

 ut.WriteHist();  
 ut.WriteCounters();  


