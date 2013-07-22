#!/usr/bin/env python

import os

from optparse import OptionParser
parser=OptionParser()
parser.add_option("-D","--dir")
parser.add_option("-d","--datfile")
parser.add_option("--eosWalk",type="int")
parser.add_option("-q","--queue",default="8nm")
parser.add_option("--dryRun",default=False,action="store_true")
(options,args)=parser.parse_args()

os.system('mkdir -p %s'%options.dir)

def writeSpec(i):
  f = open('%s/sub%d.sh'%(os.path.abspath(options.dir),i),'w')
  f.write('#!/bin/bash\n')
  f.write('rm -f %s.done\n'%(f.name))
  f.write('rm -f %s.fail\n'%(f.name))
  f.write('rm -f %s.run\n'%(f.name))
  f.write('rm -f %s.log\n'%(f.name))
  f.write('cd %s\n'%os.getcwd())
  f.write('eval `scramv1 runtime -sh`\n')
  f.write('cd -\n')
  f.write('cp %s/scripts/make_hists_from_raw_files.py .\n'%os.getcwd())
  f.write('cp %s .\n'%os.path.abspath(options.datfile))
  f.write('mkdir lib\n')
  f.write('cp %s/lib/libBackgroundProfileFitting.so lib/\n'%os.getcwd())
  f.write('touch %s.run\n'%(f.name))
  subline = './make_hists_from_raw_files.py -d %s --runSpecificFiles=%d'%(os.path.basename(options.datfile),i)
  if options.eosWalk:
    subline += ' --eosWalk=%d'%options.eosWalk
  f.write('if ( %s ) then \n'%subline)
  f.write('\ttouch %s.done\n'%(f.name))
  f.write('\trm -f %s.run\n'%(f.name))
  f.write('else\n')
  f.write('\ttouch %s.fail\n'%(f.name))
  f.write('fi\n')
  f.close()
  os.system('chmod +x %s'%f.name)
  if not options.dryRun: os.system('bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name))
  else: print 'bsub -q %s -o %s.log %s'%(options.queue,f.name,f.name)

df = open(options.datfile)
i=0
for line in df.readlines():
  if line.startswith('#') or line=='' or line =='\n': continue
  writeSpec(i)
  i+=1
