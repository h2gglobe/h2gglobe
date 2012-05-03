#!/bin/env python

print "Hello"

from ROOT import gROOT
gROOT.SetBatch(True)

from ROOT import TFile, gDirectory

from lip.Tools.myopt import options

from sys import argv
import os
import sys
import json
import numpy as np

print "We loaded all libraries"

# Options handled in python
o = options( [
    ("isdata",False,""),
    
    ("refLumi",1.,""),
    ("pufile","",""),

    ("puhisto","",""),
    ("npu","pu_n",""),
    ("evsel","",""),

    ("outfile","",""),
    
    ("samples_db","https://spreadsheets.google.com/pub?key=0AkBaWslL8vv3dGpOZXhzdnBlN3RCdy1NcnVWQWIzOEE&hl=en&output=csv&gid=0","CSV database with samples location"),
    ("xsec_db",   "https://spreadsheets.google.com/pub?key=0AkBaWslL8vv3dGpOZXhzdnBlN3RCdy1NcnVWQWIzOEE&hl=en&output=csv&gid=1", "CSV database with cross sections"),
    
    ("production", "Dec", "Production id in samples DB"),
    ("samples", "", "List of samples to analyse"),
    ("splitbins", False, ""),
    ("printonly", False, ""),

    ("maxEvent", -1, "")

    ]
             )

# parse options
o.parse(argv)

print "---------------------------------------------------------"
o.dump()

# make sure ROOT doesn't parse the command line options
sys.argv = []

# load ROOT
print "Loading ROOT libraries..."
from lip.Tools.rootutils import *
from lip.Tools.ntupleutils import *

# trees
print "Loading ntuples..."
datasets = []
if o.samples != "" :
    selection = []
    for s in o.samples.split(";"):
        (sample,bins) = s.split(":")
        bins = bins.split(",")
        selection.append((sample,bins))
    o.isdata = loadFromDb(datasets,o.production,o.refLumi,selection,o.splitbins,o.samples_db,o.xsec_db,o.maxEvent,o.printonly,useTChain=True)
elif o.infile != "":
    o.isdata = loadOneFile(datasets,o.infile,data)
else:
    print "Error: no input specified."
    sys.exit(1)

# unweight if required
for d in datasets:
    d.SetWeight(1.,"global")

# Load analysis libraries
fpu = ROOT.TFile.Open(o.pufile)
puOrig = fpu.Get(o.puhisto)

fout = ROOT.TFile.Open(o.outfile,"recreate")

for d in datasets:
    dsname = d.GetTitle()
    print "Computing weights for %s " %  dsname
#    for f in d.GetListOfFiles():
#        print  f.GetTitle()

#    fname0 = d.GetListOfFiles()[0]
#    f0 = TFile.Open(fname0.GetTitle())
    firstFileIdx=0
    files = d.GetListOfFiles()
    while firstFileIdx < len(files):
        fname0 = files[firstFileIdx]
        f0 = TFile.Open(fname0.GetTitle())
        if f0:
            break
        firstFileIdx+=1


    puf0 = f0.Get("pileup") 
    if puf0:
        useHistos=True
        if puf0.GetNbinsX()<puOrig.GetNbinsX():
            nbins = puf0.GetNbinsX()
        elif puf0.GetNbinsX()>puOrig.GetNbinsX():
            print "Unhandled case: Gen PU in the samples has more bins than the Target PU one.\n Quitting."
            sys.exit(-1);
        else:
            nbins = False
    else:
        useHistos=False
        
    f0.Close()
    
    fout.cd()
    fout.mkdir(dsname)
    fout.cd(dsname)
 
    if nbins:        
        pu = puOrig.Rebin(nbins,'',np.arange(-0.5,nbins+0.5,1.0))
    else:
        pu = puOrig
 
    samplepu = pu.Clone("generated_pu")
    samplepu.SetTitle(dsname)
    samplepu.Reset("ICE")
    samplepu.SetEntries(0)
    
    samplewei = pu.Clone("weights")
    samplewei.SetTitle(dsname)
    samplewei.Reset("ICE")
    samplewei.SetEntries(0)

    if useHistos:
        print "Hey! These files have the pileup histograms."
        for fname in d.GetListOfFiles():
            f = TFile.Open(fname.GetTitle())
            if f:
                print " But hey! I'll add also %s" % fname.GetTitle()
                samplepu.Add(f.Get('pileup'))
                f.Close()
            else:
                print " Error opening %s" % fname.GetTitle()
    else:
        print "Hey! These files have no pileup histo. Making it by looping over all events."
        draw_expr = "%s>>%s" % ( o.npu, samplepu.GetName() )
        print draw_expr
        d.Draw(draw_expr, o.evsel, "goff" )
    
#    pu.Print()
#    samplepu.Print()
    samplewei.Divide( pu, samplepu, 1./pu.Integral(), 1./samplepu.Integral() )

    max = 0.
    integral = samplepu.Integral()
    for b in range(1,samplepu.GetNbinsX()):
        content = samplepu.GetBinContent(b)
        if content / integral < 0.01:
            continue
        wei = samplewei.GetBinContent(b)
        if wei > max:
            max = wei
    ## samplewei.Scale( 1. / max )

    eff = samplewei.Clone("eff")
    eff.Multiply( samplepu )
    print "Unweighing efficiency: %1.3f; Weights integral %1.2f" % (eff.Integral() / samplepu.Integral(), samplewei.Integral())
    
    fout.cd(dsname)
    samplepu.Write()
    samplewei.Write()

fout.cd()
pu.Write("target_pu")
samplepu.Write()
samplewei.Write()

fout.Close()

# fix crash in TChain destructor when using TNetFiles 
for d in datasets : d.Reset()

print "All done!"
