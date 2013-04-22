#!/usr/bin/env python
import sys
import os

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n","--njobs",dest="njobs",type="int",)
parser.add_option("-t","--toysperjob",dest="toyspj",type="int")
parser.add_option("-f","--datfile",dest="datfile",type="string")
parser.add_option("-q","--queue",dest="queue",type="string")
parser.add_option("-d","--directory",dest="directory",type="string",default=".")
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
  os.system('cp %s/bin/diySeparation %s/'%(os.getcwd(),dir))
  f = open('%s/sub%d.sh'%(dir,j),'w')
  f.write('#!/bin/bash\n')
  if dir[0] == '/':
    f.write('cd %s\n'%(dir))
  else:
    f.write('cd %s/%s\n'%(os.getcwd(),dir))
  f.write('eval `scramv1 runtime -sh`\n')
  f.write('touch sub%d.sh.run\n'%j)
  subline = ''
  if options.datfile[0] == '/':
    subline = './diySeparation %s %d'%(options.datfile,options.toyspj)
  else:
    subline = './diySeparation %s/%s %d'%(os.getcwd(),options.datfile,options.toyspj)
  f.write('if ( %s )\n'%subline)
  f.write('\tthen touch sub%d.sh.done\n'%j)
  f.write('\telse touch sub%d.sh.fail\n'%j)
  f.write('fi\n')
  f.write('rm -f sub%d.sh.run\n'%j)
  os.system('chmod +x %s'%f.name)

  filename=f.name
  if filename[0] != '/':
    filename = os.getcwd()+'/'+filename
  if not options.dryRun: os.system('bsub -q %s -o %s.log %s'%(options.queue,filename,filename))
  #else: print 'bsub -q %s -o %s/%s.log %s/%s'%(options.queue,os.getcwd(),f.name,os.getcwd(),f.name)
