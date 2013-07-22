#!/usr/bin/env python

import os
import sys
import fnmatch

dir=sys.argv[1]

# check completed jobs
comp={}
fail={}
for root,dirs,files in os.walk(dir):
  for file in fnmatch.filter(files,'*.sh'):
    cat = int(file.split('cat')[-1].split('_')[0])
    job = int(file.split('job')[-1].split('.sh')[0])
    if os.path.isfile('%s/%s.done'%(dir,file)):
      if cat not in comp:
        comp[cat] = [job]
      else:
        comp[cat].append(job)
    if os.path.isfile('%s/%s.fail'%(dir,file)):
      if cat not in fail:
        fail[cat] = [job]
      else:
        fail[cat].append(job)

print 'Completed jobs:'
for key,item in comp.items():
  print '\t cat', key, ' - ', len(item)

print 'Failed jobs:'
for key,item in fail.items():
  print '\t cat', key, ' - ', len(item)

raw_input('Continue?')

import ROOT as r
# merge jobs
for cat, jobs in comp.items():
  
  ofile = r.TFile('%s/BiasStudyOut_cat%d.root'%(dir,cat),'RECREATE')
  chain = r.TChain("muTree_cat%d"%cat)
  for job in jobs:
    chain.Add('%s/BiasStudyOut_cat%d_job%d.root/muTree'%(dir,cat,job))

  chain.Write()
  ofile.Close()

