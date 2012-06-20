#!/bin/env python

from sys import argv

target=argv[1]

if ".root" in target:
    target = target.replace(".root","")
elif not "/" in target:
    target = "%s/CMS-HGG" % target


import ROOT
ROOT.gSystem.SetIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include")

ROOT.gSystem.Load("../libLoopAll.so")
ROOT.gROOT.LoadMacro("InterpolateMass.C+")

ROOT.InterpolateMassRange(110,150,0.5,target)
