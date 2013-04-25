#!/usr/bin/env python

import sys
import os
import fnmatch

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-d","--dir",dest="dir")
parser.add_option("-m","--mass",dest="mass",type="float")
(options,args) = parser.parse_args()

njobs=0
for root,dirs,files in os.walk(options.dir):
  #for filename in fnmatch.filter(files,'higgsCombineTest.HybridNew.mH125*.root'):
  for filename in fnmatch.filter(files,'sub*.sh'):
    if root==options.dir:
      njobs+=1

ndone=0
for job in range(njobs):
  if os.path.isfile('%s/sub%d.sh.done'%(options.dir,job)):
    ndone+=1

print ndone, '/', njobs, 'finished'

if ndone==njobs:
    
  import ROOT as r
  r.gROOT.ProcessLine('.L $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTree.cxx')
  from ROOT import hypoTestResultTree
  
  for job in range(njobs):
    filename = "higgsCombinejob%d.HybridNew.mH%3.1f.root"%(job,options.mass)
    print filename
    r.gROOT.ProcessLine('TFile::Open(\"'+options.dir+'/'+filename+'\")')

  r.gROOT.ProcessLine('hypoTestResultTree(\"qmu.root\",%3.1f,1.,\"x\")'%(options.mass))
