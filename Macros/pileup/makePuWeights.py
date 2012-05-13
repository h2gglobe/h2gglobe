#!/bin/env python

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
from sys import argv
argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
argv.remove( '-b-' )

import os, re
from pprint import pprint

from sampleDict import sampleString
from fsUtils import canonicalize

from optparse import OptionParser
parser = OptionParser(usage="usage: %prog target_pileup.root mc_pileup.root [mc_pileup2.root ...] \nrun with --help to get list of options")
(options, args) = parser.parse_args()

if len(args) < 2:
    parser.print_usage()
    raise RunTimeError("Need at least 2 files.")

fNames = [canonicalize(x) for x in args]

dataPuFileName = fNames[0]
mcPuFileNames = fNames[1:]

from ROOT import TFile, gDirectory

def getPuHistos(tFile):
    fileObjs = list(tFile.GetListOfKeys())
    return [ tFile.Get(obj.GetName()) for obj in fileObjs if re.search('pileup', obj.GetName()) ]


dataPuFile = ROOT.TFile.Open(dataPuFileName)
dataPuHistosDict = dict( [(h.GetName(),h) for h in getPuHistos(dataPuFile)] )

fout = ROOT.TFile.Open("weights.root","recreate")
fout.cd()

for (name, hist) in dataPuHistosDict.iteritems():
    hist.Write('target_'+name)


def makeWeightHisto(data, mc):
    w = data.Clone(data.GetName()+'_weights')
    w.SetTitle(w.GetName())
    w.Reset("ICE")
    w.SetEntries(0)

    #TODO: check that the axis dimensions are compatible
#         if puf0.GetNbinsX()<puOrig.GetNbinsX():
#             nbins = puf0.GetNbinsX()
#         elif puf0.GetNbinsX()>puOrig.GetNbinsX():
#             print "Unhandled case: Gen PU in the samples has more bins than the Target PU one.\n Quitting."
#             sys.exit(-1);
#         else:
#             nbins = False
            
#         if nbins:
#import numpy as np
#             pu = puOrig.Rebin(nbins,'',np.arange(-0.5,nbins+0.5,1.0))
#         else:
#             pu = puOrig
                
    w.Divide( data, mc, 1./data.Integral(), 1./mc.Integral() )

#     eff = w.Clone('eff')
#     eff.Multiply( mc )
#     print "Unweighing efficiency: %1.3f; Weights integral %1.2f" % (eff.Integral() / mc.Integral(), w.Integral())

    return w
    
for mcPuFileName in mcPuFileNames:
    sampleName = sampleString(mcPuFileName)
    print sampleName
    
    mcPuFile = ROOT.TFile.Open(mcPuFileName)
    mcPuHistos = getPuHistos(mcPuFile)

    histToSave = list()
    for mcPuHisto in mcPuHistos:
        try:
            dataPuHisto = dataPuHistosDict[mcPuHisto.GetName()]
        except KeyError:
            continue

        histToSave += [ mcPuHisto, makeWeightHisto(dataPuHisto, mcPuHisto) ]

    if histToSave:
        fout.mkdir(sampleName)
        fout.cd(sampleName)
        for h in histToSave:
            h.Write()

    mcPuFile.Close()

fout.Close()


