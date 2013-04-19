#!/usr/bin/env python

import os
import sys
import fnmatch

baseDirectory = sys.argv[1]
subDirs = [x[1] for x in os.walk(baseDirectory)]
subDirs = [baseDirectory+"/"+x for x in subDirs[0]]

home='$CMSSW_BASE/src/h2gglobe/SpinAnalysis'

for dir in subDirs:
  #jobs = [x[1] for x in os.walk(dir)]
  #jobs = [dir+"/"+x for x in jobs[0]]
  jobs = glob.glob(dir+"/job*")
  print "Checking "+dir+" with %d jobs."%len(jobs)

  processCommand = '%s/checkDIY.py %d %s'%(home,len(jobs),dir)
  #print processCommand
  os.system(processCommand)
