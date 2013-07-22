#!/usr/bin/env python

import os
import sys
import fnmatch

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-D","--dir")
parser.add_option("-o","--outfile")
parser.add_option("-d","--datfile")
parser.add_option("--expectSignal",type="float",default=0.)
parser.add_option("--runSpecificFiles",type="str")
parser.add_option("--eosWalk",type="int",help="Walk eos for files - provide njobs")
(options,args)=parser.parse_args()

import ROOT as r

r.gSystem.Load('lib/libBackgroundProfileFitting.so')

profiler = r.ProfileMultiplePdfs()

def makeHists(cat=0,meanB=50,meanL=-4.,meanH=4.,errB=50,errL=0.5,errH=1.5,pullB=50,pullL=-5.,pullH=5.):

  dir = options.dir
  list_of_files=[]

  if '/store/' in dir and options.eosWalk:
    dir = dir.strip('~')
    if '/eos/cms' not in dir:
      dir = '/eos/cms'+dir
    dir = 'root://eoscms/'+dir
  
  if '/store/' in options.outfile and options.eosWalk:
    options.outfile = options.outfile.strip('~')
    if '/eos/cms' not in options.outfile:
      options.outfile = '/eos/cms'+options.outfile
    options.outfile = 'root://eoscms/'+options.outfile
  
  if options.eosWalk:
    for job in range(options.eosWalk):
      list_of_files.append('BiasStudyOut_cat%d_job%d.root'%(cat,job))
  else:
    for root, dirs, files in os.walk(dir):
      for filename in fnmatch.filter(files,'*job*.root'):
        list_of_files.append(filename)

  outfile = r.TFile.Open('%s'%options.outfile,'RECREATE')

  truth_models = set()

  test_file = r.TFile.Open(dir+'/'+list_of_files[0])
  for key in test_file.GetListOfKeys():
    if 'truth' not in key.GetName(): continue
    truth= key.GetName().split('truth')[1].split('_cat')[0][1:]
    truth_models.add(truth)

  #print truth_models

  types = set()
  types.add('Fab')
  types.add('Paul')
  types.add('Chi2')
  types.add('AIC')

  util={}
  util['fab'] = 'Fab'
  util['paul'] = 'Paul'
  util['chi2'] = 'Chi2'
  util['aic'] = 'AIC'

  coverageValues=[0.1,0.5,1.,2.]

  histMap={}
  histErrMap={}
  histPullMap={}
  graphCovMap={}
  counterCovMap={}

  for type in types:
    histMap[type] = {}
    histErrMap[type] = {}
    histPullMap[type] = {}
    graphCovMap[type] = {}
    counterCovMap[type] = {}
    for mod in truth_models:
      histMap[type][mod] = r.TH1F('%s_mu%s'%(mod,type),'%s_mu%s'%(mod,type),meanB,meanL,meanH)
      histErrMap[type][mod] = r.TH1F('%s_mu%sErr'%(mod,type),'%s_mu%sErr'%(mod,type),errB,errL,errH)
      histPullMap[type][mod] = r.TH1F('%s_mu%sPull'%(mod,type),'%s_mu%sPull'%(mod,type),pullB,pullL,pullH)
      graphCovMap[type][mod] = []
      counterCovMap[type][mod] = []
      for c, cov in enumerate(coverageValues):
        graphCovMap[type][mod].append(r.TGraphErrors())
        graphCovMap[type][mod][c].SetName('%s_mu%sCov%3.1f'%(mod,type,cov))
        counterCovMap[type][mod].append([0,0])

  print dir
  for i, f in enumerate(list_of_files):
    print '\tJob', i+1,'/',len(list_of_files), '\r',
    sys.stdout.flush()
    file = r.TFile.Open(dir+'/'+f)
    for key in file.GetListOfKeys():
      graph = key.ReadObj()
      if 'Envelope' not in graph.GetName(): continue
      truth = graph.GetName().split('truth')[1].split('_cat')[0][1:]
      type = graph.GetName().split('Envelope')[0]
      mytype = util[type]
     
      muInfo = profiler.getMinAndErrorAsymmVec(graph,1.)
      muVal = muInfo.at(0)
      err_low = muInfo.at(1)
      err_high = muInfo.at(2)
      sym_err = (err_low+err_high)/2.
      pull = profiler.getPull(graph,options.expectSignal)

      #print truth, mytype, '%4.2f  %4.2f  %4.2f  %4.2f'%(muVal,err_low,err_high,sym_err)
      
      if muVal>=999. or sym_err>=999. or sym_err<0.001: continue

      histMap[mytype][truth].Fill(muVal)
      histErrMap[mytype][truth].Fill(sym_err)
      histPullMap[mytype][truth].Fill(pull)
      
      for c, cov in enumerate(coverageValues):
        counterCovMap[mytype][truth][c][1] += 1
        if r.TMath.Abs(pull)<cov:
          counterCovMap[mytype][truth][c][0] += 1
  
  print '\n',
  outfile.cd()

  for type, item in histMap.items():
    for truth, hist in item.items():
      hist.Write()
  for type, item in histErrMap.items():
    for truth, hist in item.items():
      hist.Write()
  for type, item in histPullMap.items():
    for truth, hist in item.items():
      hist.Write()
  for type, item in counterCovMap.items():
    for truth, covArray in item.items():
      for c, covCounter in enumerate(covArray):
        nPass = float(covCounter[0])
        nTotal = float(covCounter[1])
        covVal = nPass/nTotal
        covErr = r.TMath.Sqrt((covVal*(1.-covVal))/nTotal)
        graphCovMap[type][truth][c].SetPoint(0,options.expectSignal,covVal)
        graphCovMap[type][truth][c].SetPointError(0,0.5,covErr)
        graphCovMap[type][truth][c].Write()
  outfile.Close()

def countLines(filename):
  f = open(filename)
  count=0
  for line in f.readlines():
    if line.startswith('#') or line=='' or line =='\n': continue
    else: count+=1
  f.close()
  return count

# --main-- here
if not options.datfile:
  makeHists()
else:
  nfiles = countLines(options.datfile)
  f = open(options.datfile)
  sw = r.TStopwatch()
  list_of_specifics=[]
  if options.runSpecificFiles:
    for spec in options.runSpecificFiles.split(','):
      list_of_specifics.append(int(spec))
  i=0
  for line in f.readlines():
    if line.startswith('#') or line=='' or line =='\n': continue
    if options.runSpecificFiles and i not in list_of_specifics: 
      i+=1
      continue 
    info = line.split()
    cat = int(info[0])
    options.dir = info[1]
    options.outfile = info[2]
    options.expectSignal = float(info[3])
    sw.Reset()
    sw.Start()
    print 'Running file', i, 'of', nfiles
    makeHists(cat,int(info[4]),float(info[5]),float(info[6]),int(info[7]),float(info[8]),float(info[9]),int(info[10]),float(info[11]),float(info[12]))
    sw.Stop()
    print 'Took:', sw.Print()
    i+=1
  f.close()
    
