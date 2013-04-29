#!/usr/bin/env python

import sys
import os
import fnmatch

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-M","--method",dest="methods",default=[],action="append")
parser.add_option("-q","--qqbarPoints",dest="qqbarPoints",default=[],action="append")
parser.add_option("-d","--dir",dest="dir")
parser.add_option("-m","--mass",dest="mass",type="float")
(options,args) = parser.parse_args()

def doqqbar():
  
  njobs=0
  for root,dirs,files in os.walk(options.dir):
    for filename in fnmatch.filter(files,'sub*qqbar.sh'):
      if root==options.dir:
        njobs+=1

  ndone=0
  for job in range(njobs):
    if os.path.isfile('%s/sub%dqqbar.sh.done'%(options.dir,job)):
      ndone+=1

  print 'qqbar has', ndone, '/', njobs, 'jobs finished'

  if ndone==njobs:
    import ROOT as r
    r.gROOT.ProcessLine('.L $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTree.cxx')
    from ROOT import hypoTestResultTree
    
    for i,p in enumerate(options.qqbarPoints):
    
      for job in range(njobs):
        filename = "higgsCombinejob%dfqq%1.2f.HybridNew.mH%3.1f.root"%(job,p,options.mass)
        print filename
        r.gROOT.ProcessLine('TFile::Open(\"'+options.dir+'/'+filename+'\")')

      r.gROOT.ProcessLine('hypoTestResultTree(\"%s/qmu_qqbar%1.2f.root\",%3.1f,1.,\"x\")'%(options.dir,p,options.mass))
      
      r.gROOT.CloseFiles()

def doSeparation():
  
  njobs=0
  for root,dirs,files in os.walk(options.dir):
    for filename in fnmatch.filter(files,'sub*Sep.sh'):
      if root==options.dir:
        njobs+=1

  ndone=0
  for job in range(njobs):
    if os.path.isfile('%s/sub%dSep.sh.done'%(options.dir,job)):
      ndone+=1

  print 'Separation has', ndone, '/', njobs, 'jobs finished'

  if ndone==njobs:
      
    import ROOT as r
    r.gROOT.ProcessLine('.L $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTree.cxx')
    from ROOT import hypoTestResultTree
    
    for job in range(njobs):
      filename = "higgsCombinejob%dsep.HybridNew.mH%3.1f.root"%(job,options.mass)
      print filename
      r.gROOT.ProcessLine('TFile::Open(\"'+options.dir+'/'+filename+'\")')

    r.gROOT.ProcessLine('hypoTestResultTree(\"%s/qmu_sep.root\",%3.1f,1.,\"x\")'%(options.dir,options.mass))

if len(options.methods)==0: options.methods=['ChannelCompatibility','Separation','qqbar']
if len(options.qqbarPoints)==0: options.qqbarPoints=[0.,0.25,0.5,0.75,1.0]
if 'Separation' in options.methods: doSeparation()
if 'qqbar' in options.methods: doqqbar()
