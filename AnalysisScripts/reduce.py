#!/usr/bin/env python

import ROOT
from python.configProducer import *
from python.runOptions import *
from sys import argv, exit
(options,args)=parser.parse_args()

ROOT.gROOT.SetBatch(1)
## ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.gSystem.Load("../libLoopAll.so")

ROOT.gBenchmark.Start("Reduction")

config_file="filestoreduce.dat"
if options.inputDat:
    config_file = options.inputDat

ut = ROOT.LoopAll();
cfg = configProducer(ut,config_file,1,int(options.nJobs),int(options.jobId),debug=options.verbose)

if not options.dryRun:
  ut.LoopAndFillHistos()
ROOT.gBenchmark.Show("Reduction")
