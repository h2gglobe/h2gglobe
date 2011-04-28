#!/usr/bin/env python

import ROOT
from python.configProducer import *

from sys import argv
config_file="filestoreduce.dat"
if len(argv) > 1:
    config_file = argv[1]

#ROOT.gSystem.Load("libRooFit.so")
ROOT.gSystem.Load("libPhysics.so");
ROOT.gSystem.Load("libCore.so");
ROOT.gSystem.Load("../libLoopAll.so");

ROOT.gBenchmark.Start("Reduction")

ut = ROOT.LoopAll();
cfg = configProducer(ut, config_file,1)
ut.LoopAndFillHistos()

ROOT.gBenchmark.Show("Reduction");
