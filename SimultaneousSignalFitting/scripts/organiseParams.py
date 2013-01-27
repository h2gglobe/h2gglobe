#!/usr/bin/env python

import sys
import numpy
import os

datfile = open(sys.argv[1])
mh = int(sys.argv[2])

for line in datfile.readlines():
  if line.startswith('#') or line=='\n': continue
  if line.startswith('infile'): continue
  line_els  = line.split()
  if len(line_els)!=7: continue
  proc      = line_els[1]
  cat       = int(line_els[2])
  os.system('mv dat/initFit_%s_cat%d.dat dat/initFit_%s_cat%d_dump.dat'%(proc,cat,proc,cat))

  paramfile = open('dat/initFit_%s_cat%d_dump.dat'%(proc,cat))
  newoutfile = open('dat/initFit_%s_cat%d.dat'%(proc,cat),'w')
  mylines=[]
  for line in paramfile.readlines():
    if 'mh%d'%mh in line: mylines.append(line)

  for m in numpy.arange(110,151,5):
    for line in mylines:
      newoutfile.write(line.replace('mh125','mh%d'%m))

  newoutfile.close()
  paramfile.close()


