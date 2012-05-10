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

ROOT.gROOT.SetBatch(1)

ROOT.gBenchmark.Start("Analysis");

config_file="filestocombine.dat"
if options.inputDat:
    config_file = options.inputDat

ut = ROOT.LoopAll();
cfg = configProducer(ut,config_file,0,int(options.nJobs),int(options.jobId))
#cfg = configProducer(ut,config_file,0,-1,0)

if not options.dryRun:
  ut.MergeContainers()
  ROOT.gBenchmark.Show("Analysis")
  ut.WriteFits()

ROOT.gROOT.Reset()
