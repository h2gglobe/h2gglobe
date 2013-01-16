#!/usr/bin/env python

import os
import sys
import fnmatch

njobs=int(sys.argv[1])
dir=sys.argv[2]

completed_jobs=[]
failed_jobs=[]
haddline="hadd -f %s/qmu.root"%dir

for j in range(njobs):
  if os.path.exists('%s/job%d/sub%d.sh.done'%(dir,j,j)):
    print 'Job %d complete'%j
    completed_jobs.append(j)
  if os.path.exists('%s/job%d/sub%d.sh.log'%(dir,j,j)) and not os.path.exists('%s/job%d/sub%d.sh.done'%(dir,j,j)):
    failed_jobs.append(j)

var = raw_input('%d jobs failed. Resubmit? [y]/[n]\n'%len(failed_jobs))
if var=='y' or var=='Y' or var =='yes':
  for j in failed_jobs:
    os.system('bsub -q 1nh -o %s/%s/job%d/sub%d.sh.log %s/%s/job%d/sub%d.sh'%(os.getcwd(),dir,j,j,os.getcwd(),dir,j,j))

if len(completed_jobs)!=njobs:
  print 'Only %d of %d jobs complete'%(len(completed_jobs),njobs)
  var = raw_input('Do you want to hadd what\'s there? [y]/[n]\n')
  if var=='n' or var=='N' or var=='no': sys.exit()

for j in completed_jobs:
  for root, dirs, files in os.walk('%s/job%d'%(dir,j)):
    matched = fnmatch.filter(files,'LLOut.root')
    if len(matched)==0:
      sys.exit('No files LLOut.root found')
    elif len(matched)>1:
      sys.exit('Found mulitple matching files')
      print matched
    else:
      haddline += ' %s/job%d/LLOut.root'%(dir,j)

print haddline
var = raw_input('You definitely want to hadd? Will overwrite %s/qmu.root - [y]/[n]\n'%dir)
if var=='n' or var=='N' or var=='no': sys.exit()
os.system(haddline)
