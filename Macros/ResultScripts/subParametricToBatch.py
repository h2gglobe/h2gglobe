#!/usr/bin/env python

import os
import sys
import numpy
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p","--path",dest="path")
parser.add_option("-d","--datacard",dest="datacard")
parser.add_option("-m","--method",dest="method")
parser.add_option("-q","--queue",dest="queue")
parser.add_option("-e","--expectedOnly",dest="expectedOnly",default=False,action="store_true")
parser.add_option("-L","--mHlow",dest="mHlow",type="int",default=110)
parser.add_option("-H","--mHhigh",dest="mHhigh",type="int",default=150)
parser.add_option("","--runIC",dest="runIC",default=False,action="store_true")

(options,args)=parser.parse_args()

path = options.path 
datacard = options.datacard
method = options.method
queue = options.queue
expectedOnly = options.expectedOnly
ml = options.mHlow
mh = options.mHhigh

if method!='Asymptotic' and method!='AsymptoticNew' and method!='ProfileLikelihood' and method!='ExpProfileLikelihood':
  print method, 'is invalid'
  sys.exit()

if method=='ProfileLikelihood' and expectedOnly:
  method='ExpProfileLikelihood'
  expectedOnly=False

if method=='ExpProfileLikelihood':
  expectedOnly=False

if not os.path.isdir('%s/%s'%(path,method)):
  os.makedirs('%s/%s'%(path,method))

files=[]

for m in numpy.arange(ml,mh+1,1):
  f = open("%s/%s/%s%d.sh"%(path,method,method,m),"w")
  f.write('#!/bin/bash\n')
  f.write('cd %s\n'%path)
  f.write('eval `scramv1 runtime -sh`\n')
  line = '%s/runCombine.sh %s %s %s %d'%(path,path,method,datacard,m)
  if method=='ProfileLikelihood':
    line += ' --signif --pvalue'
  if method=='ExpProfileLikelihood':
    line += ' --signif --pvalue -t -1 --expectSignal=1'
  if expectedOnly:
    line += ' --run=expected'
  f.write(line+'\n')
  f.close()
  os.system('chmod +x %s'%f.name)
  files.append(f)

print 'Done'
print 'Submit?'
a = raw_input()

if a!='No' and a!='no' and a!='n' and a!='N':
  for f in files:
    os.system('rm %s.log'%(f.name))
    os.system('bsub -q %s -o %s.log %s'%(queue,f.name,f.name))  
    #print 'bsub -q %s -o %s.log %s'%(queue,f.name,f.name)  

