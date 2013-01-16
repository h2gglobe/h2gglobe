#!/usr/bin/env python
import sys
import os

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n","--njobs",dest="njobs",type="int",)
parser.add_option("-t","--toysperjob",dest="toyspj",type="int")
parser.add_option("-f","--file",dest="file",type="string")
parser.add_option("-q","--queue",dest="queue",type="string")
parser.add_option("-d","--directory",dest="directory",type="string",default=".")
parser.add_option("-b","--bdtcats",dest="bcats",type="int",default=4)
parser.add_option("-s","--spincats",dest="scats",type="int",default=2)
parser.add_option("--dryRun",dest="dryRun",action="store_true",default=False)
parser.add_option("--clean",dest="clean",action="store_true",default=False)
(options,args) = parser.parse_args()

for j in range(options.njobs):
  dir = '%s/job%d'%(options.directory,j)
  os.system('mkdir -p %s'%dir)
  os.system("rm -f %s/*.sh.fail"%dir)
  os.system("rm -f %s/*.sh.done"%dir)
  os.system("rm -f %s/*.sh.log"%dir)
  os.system("rm -f %s/*.sh.run"%dir)
  #os.system('cp %s/%s %s/job%d\n'%(os.getcwd(),options.file,os.getcwd(),j))
  os.system('cp %s/bin/diySeparation %s/%s'%(os.getcwd(),os.getcwd(),dir))
  f = open('%s/sub%d.sh'%(dir,j),'w')
  f.write('#!/bin/bash\n')
  f.write('cd %s/%s\n'%(os.getcwd(),dir))
  f.write('eval `scramv1 runtime -sh`\n')
  f.write('touch sub%d.sh.run\n'%j)
  f.write('if ( ./diySeparation %d %s/%s %d %d)\n'%(options.toyspj,os.getcwd(),options.file,options.bcats,options.scats))
  f.write('\tthen touch sub%d.sh.done\n'%j)
  f.write('\telse touch sub%d.sh.fail\n'%j)
  f.write('fi\n')
  f.write('rm -f sub%d.sh.run\n'%j)
  os.system('chmod +x %s'%f.name)
  if not options.dryRun: os.system('bsub -q %s -o %s/%s.log %s/%s'%(options.queue,os.getcwd(),f.name,os.getcwd(),f.name))
  #else: print 'bsub -q %s -o %s/%s.log %s/%s'%(options.queue,os.getcwd(),f.name,os.getcwd(),f.name)
