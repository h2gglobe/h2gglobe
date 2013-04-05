#!/bin/env python

import sys

fin = sys.argv[1]
sys.argv.append("-b")
import ROOT

procs = [ "ggh", "vbf", "wzh", "tth" ]
ncat  = 4
ntot  = 9

nuis = {}

file = ROOT.TFile.Open(fin)

for cat in range(ncat):
    nuis[cat]={}
    for proc in procs:
        nominal = file.Get("th1f_sig_%s_mass_m125_cat%d" % (proc,cat) )
        up      = file.Get("th1f_sig_%s_mass_m125_cat%d_r9EffUp01_sigma" % (proc,cat) )
        down    = file.Get("th1f_sig_%s_mass_m125_cat%d_r9EffDown01_sigma" % (proc,cat) )
        if nominal.Integral() != 0:
            nuis[cat][proc] = (down.Integral()/nominal.Integral(), up.Integral()/nominal.Integral())
        else:
            nuis[cat][proc] = (1.,1.)

for cat in range(ncat):
    for proc in procs:
        print "%1.3f/%1.3f " % nuis[cat][proc],
    print "-     ",

for cat in range(ncat,ntot):
    for proc in procs:
        print "%1.3f/%1.3f " % (1.,1.),
    print "-     ",

print

