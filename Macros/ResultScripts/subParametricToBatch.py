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
parser.add_option("","--dryRun",dest="dryRun",default=False,action="store_true")

(options,args)=parser.parse_args()

path = options.path 
datacard = options.datacard
method = options.method
folder = options.method
queue = options.queue
expectedOnly = options.expectedOnly
ml = options.mHlow
mh = options.mHhigh

if method!='Asymptotic' and method!='AsymptoticNew' and method!='ProfileLikelihood' and method!='ExpProfileLikelihood':
  print method, 'is invalid'
  sys.exit()

if method=='ProfileLikelihood' and expectedOnly:
  method='ExpProfileLikelihood'
  folder='ProfileLikelihood'
  expectedOnly=False

if method=='ExpProfileLikelihood':
  folder='ProfileLikelihood'
  expectedOnly=False

if not os.path.isdir('%s/%s'%(path,folder)):
  os.makedirs('%s/%s'%(path,folder))

files=[]

for m in numpy.arange(ml,mh+1,1):
  f = open("%s/%s/%s%d.sh"%(path,folder,method,m),"w")
  f.write('#!/bin/bash\n')
  f.write('cd %s\n'%path)
  f.write('touch %s.run\n'%f.name)
  f.write('echo --------------------------------------------\n')
  f.write('echo   Running %s on %s at mass %d \n'%(method,datacard,m)) 
  f.write('echo --------------------------------------------\n')
  f.write('eval `scramv1 runtime -sh`\n')
  line = 'if ( combine %s/%s -M %s -m %d'%(path,datacard,method,m)
  #line = '%s/runCombine.sh %s %s %s %d'%(path,path,method,datacard,m)
  if method=='ProfileLikelihood':
    line += ' --signif --pvalue'
  if method=='ExpProfileLikelihood':
    line += ' --signif --pvalue -t -1 --expectSignal=1'
  if expectedOnly:
    line += ' --run=expected'
  f.write(line+' ) \n')
  f.write('\ttouch %s.done\n'%f.name)
  f.write('else\n')
  f.write('\t touch %s.fail\n'%f.name)
  f.write('fi\n')
  f.write('rm %s.run\n'%f.name)
  if method=='Asymptotic' or method=='ProfileLikelihood':
    f.write('mv %s/higgsCombineTest.%s.mH%d.root %s/%s/higgsCombineTest.%s.mH%d.root\n'%(path,method,m,path,folder,method,m))
  elif method=='ExpProfileLikelihood':
    f.write('mv %s/higgsCombineTest.%s.mH%d.root %s/%s/higgsCombineEXPECTED.%s.mH%d.root\n'%(path,folder,m,path,folder,folder,m))
  f.close()
  os.system('chmod +x %s'%f.name)
  files.append(f)

if options.dryRun: exit()

print 'Done'
print 'Submit?'
a = raw_input()

if a!='No' and a!='no' and a!='n' and a!='N':
  for f in files:
    if os.path.exists('%s.log'%f.name):
      os.system('rm %s.log'%(f.name))
    os.system('bsub -q %s -o %s.log %s'%(queue,f.name,f.name))  
    #print 'bsub -q %s -o %s.log %s'%(queue,f.name,f.name)  

