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
parser.add_option("-l","--lowJob",type="int",default=0)
parser.add_option("-u","--highJob",type="int",default=100000)
parser.add_option("-s","--skipJobs",default='')
parser.add_option("--splitJobs",default=False,action="store_true")
(options,args) = parser.parse_args()

if options.skipJobs=='': options.skipJobs=[]
else: options.skipJobs = [int(x) for x in options.skipJobs.split(',')]

def doqqbar():
  
  njobs=0
  for root,dirs,files in os.walk(options.dir):
    for filename in fnmatch.filter(files,'sub*qqbar.sh'):
      if 'BF' in filename: continue
      if root==options.dir:
        njobs+=1
  
  done_jobs = []
  fail_jobs = []
  run_jobs = []
  for job in range(njobs):
    if os.path.isfile('%s/sub%dqqbar.sh.done'%(options.dir,job)):
      done_jobs.append(job)
    if os.path.isfile('%s/sub%dqqbar.sh.fail'%(options.dir,job)):
      fail_jobs.append(job)
    if os.path.isfile('%s/sub%dqqbar.sh.run'%(options.dir,job)):
      run_jobs.append(job)


  print 'qqbar has', len(done_jobs), '/', njobs, 'jobs finished'
  print 'qqbar has', len(fail_jobs), '/', njobs, 'jobs failed'
  print 'qqbar has', len(run_jobs), '/', njobs, 'jobs running'
  
  if len(done_jobs)!=njobs:
    yes = raw_input('Combine anyway? [y/n]')

  if len(done_jobs)==njobs or yes=='y' or yes=='Y' or yes=='yes' or yes=='Yes':
    import ROOT as r
    r.gROOT.ProcessLine('.L $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/plotting/hypoTestResultTree.cxx')
    from ROOT import hypoTestResultTree
    
    for i,point in enumerate(options.qqbarPoints):
      p = float(point)
      for job in range(njobs):
        if job in options.skipJobs: continue
        if job<options.lowJob or job>options.highJob: continue
        if options.mass==125 or options.mass==126: filename = "higgsCombinejob%dfqq%1.2f.HybridNew.mH%3.1f.root"%(job,p,options.mass)
        else: filename = "higgsCombinejob%dfqq%1.2f.HybridNew.mH%3.1f.0.root"%(job,p,options.mass)
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
      if job in options.skipJobs: continue
      if options.mass==125 or options.mass==126: filename = "higgsCombinejob%dsep.HybridNew.mH%3.1f.root"%(job,options.mass)
      else: filename = "higgsCombinejob%dsep.HybridNew.mH%3.1f.0.root"%(job,options.mass)
      print filename
      r.gROOT.ProcessLine('TFile::Open(\"'+options.dir+'/'+filename+'\")')

    r.gROOT.ProcessLine('hypoTestResultTree(\"%s/qmu_sep.root\",%3.1f,1.,\"x\")'%(options.dir,options.mass))
    
    r.gROOT.CloseFiles()

if len(options.methods)==0: options.methods=['ChannelCompatibility','Separation','qqbar']
if len(options.qqbarPoints)==0: options.qqbarPoints=[0.,0.25,0.5,0.75,1.0]
if 'Separation' in options.methods: doSeparation()
if 'qqbar' in options.methods: doqqbar()
