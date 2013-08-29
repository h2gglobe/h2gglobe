#!/usr/bin/env python

import sys

iname = sys.argv[1]
sys.argv.append("-b")

import ROOT 
## fin = ROOT.TFile.Open("moriond_parametric/CMS-HGG_sigfit_nonuis.root")
fin = ROOT.TFile.Open(iname)

wsig = fin.Get("wsig_8TeV")
vars = wsig.allVars()
it = vars.iterator()

for i in range( vars.getSize() ):
    var = it()
    name = var.GetName()
    if "nuis" in name or "globalscale" in name:
        print "here"
        var.setVal(0.)
        var.setConstant(True)
    print name, var.isConstant()

wsig.writeToFile(iname.replace(".root","_nonuis.root"))

