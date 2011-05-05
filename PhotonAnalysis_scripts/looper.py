#!/usr/bin/env python

import ROOT
from python.configProducer import *

ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("../libLoopAll.so");

ROOT.gBenchmark.Start("Analysis");

ut = ROOT.LoopAll();
cfg = configProducer(ut,"inputfiles.dat",0)
  
ut.LoopAndFillHistos();
ROOT.gBenchmark.Show("Analysis");

ut.WriteHist();  
ut.WriteCounters();  


