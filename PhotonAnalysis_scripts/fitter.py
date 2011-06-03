#!/usr/bin/env python

import ROOT
from python.configProducer import *

ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("../libLoopAll.so");

ROOT.gBenchmark.Start("Analysis");

config_file="datafiles.dat"
if len(argv) > 1:
    config_file = argv[1]

ut = ROOT.LoopAll();
cfg = configProducer(ut,config_file,0)
  
ut.LoopAndFillHistos();
ROOT.gBenchmark.Show("Analysis");

ut.WriteFits();  


