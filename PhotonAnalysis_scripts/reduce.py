#!/usr/bin/env python

import ROOT
from python.configProducer import *

from sys import argv, exit
config_file="filestoreduce.dat"
if len(argv) > 1:
    config_file = argv[1]

njobs=-1
jobId=0
if len(argv) > 2:
    if len(argv) == 3:
        print "Usage: reduce.py <config_file> [<njobs> <job_id>]" 
        exit(1)
    njobs = int(argv[2])
    jobId = int(argv[3])
    
#ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libPhysics.so")
ROOT.gSystem.Load("libCore.so")
ROOT.gSystem.Load("../libLoopAll.so")

ROOT.gBenchmark.Start("Reduction")

ut = ROOT.LoopAll();
cfg = configProducer(ut, config_file,1,njobs,jobId)
ut.LoopAndFillHistos()

ROOT.gBenchmark.Show("Reduction")
