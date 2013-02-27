#!/bin/env python

import sys

fin = sys.argv[1]
syst = sys.argv[2]
sys.argv.append("-b")
import ROOT
ROOT.gROOT.ProcessLine(".L quadInterpolate.C+")
from ROOT import quadInterpolate

procs = [ "ggh", "vbf", "wzh", "tth" ]
ncat  = 4
ntot  = 9

nuis = {}

file = ROOT.TFile.Open(fin)

for cat in range(ncat):
    nuis[cat]={}
    for proc in procs:
        nominal = file.Get("th1f_sig_%s_mass_m125_cat%d" % (proc,cat) )
        up      = file.Get("th1f_sig_%s_mass_m125_cat%d_%sUp01_sigma" % (proc,cat,syst) )
        down    = file.Get("th1f_sig_%s_mass_m125_cat%d_%sDown01_sigma" % (proc,cat,syst) )
        if nominal.Integral() != 0:
            downE = quadInterpolate(1.,-3.,0.,3.,down.Integral(),nominal.Integral(),up.Integral())
            upE = quadInterpolate(-1.,-3.,0.,1.,down.Integral(),nominal.Integral(),up.Integral())
            nuis[cat][proc] = (downE,upE) 
        else:
            nuis[cat][proc] = (1.,1.)

print syst,  
for cat in range(ncat):
    for proc in procs:
        print "%1.3f/%1.3f " % nuis[cat][proc],
    print "-     ",

for cat in range(ncat,ntot):
    for proc in procs:
        print "%1.3f/%1.3f " % (1.,1.),
    print "-     ",

print

