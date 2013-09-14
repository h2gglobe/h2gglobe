#!/usr/bin/env python

import os
import sys
import fnmatch

dir = sys.argv[1]

njobs=0
for root,dirs,files in os.walk(dir):
  for filename in fnmatch.filter(files,'sub*qqbar.sh'):
    if 'BF' in filename: continue
    if root==dir:
      njobs+=1

print njobs

resub=0
for job in range(njobs):
  if not os.path.isfile('%s/sub%dqqbar.sh.done'%(dir,job)):
    script = os.path.abspath('%s/sub%dqqbar.sh'%(dir,job))
    print script
    os.system('bsub -q 1nd -o %s.log %s'%(script,script))
    resub+=1

print resub
