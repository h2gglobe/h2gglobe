#!/usr/bin/env python

import ROOT
from sys import argv
from python.configProducer import *
from python.lumi import *
import sys
from python.runOptions import parser
parser.add_option("-f","--files",dest="files",type="string",action="store",default="")
parser.add_option("-o","--output",dest="output",type="string",action="store",default="")
(options,args)=parser.parse_args()
if (int(options.nJobs) > 0) and (int(options.jobId) >= int(options.nJobs)):
  sys.exit("Job id's must run from 0 -> %d when splitting into %d jobs"%(int(options.nJobs)-1,int(options.nJobs)))

ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
#ROOT.gSystem.Load("$CMSSW_BASE/lib/slc5_amd64_gcc472/libHiggsAnalysisCombinedLimit.so")
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
ROOT.gSystem.Load("../libLoopAll.so");

ROOT.gROOT.SetBatch(1)

ROOT.gBenchmark.Start("Analysis");

config_file="filestocombine.dat"
if options.inputDat:
    config_file = options.inputDat
fnames=[]
if options.files != "":
  fnames=[f.lstrip(" ").rstrip(" ") for f in options.files.split(',') if f != ""]

ut = ROOT.LoopAll();
cfg = configProducer(ut,config_file,0,int(options.nJobs),int(options.jobId),files=fnames,histfile=options.output,mountEos=options.mountEos,debug=options.verbose)
#cfg = configProducer(ut,config_file,0,-1,0)

if not options.dryRun:
  ut.MergeContainers()
  ROOT.gBenchmark.Show("Analysis")
  ut.WriteFits()

ROOT.gROOT.Reset()
