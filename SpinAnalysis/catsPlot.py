#!/usr/bin/env python

import os
import re
import sys
import glob
import fnmatch

baseDirectory = sys.argv[1]
subDirs = [x[1] for x in os.walk(baseDirectory)]
subDirs = [baseDirectory+"/"+x for x in subDirs[0]]

home='$CMSSW_BASE/src/h2gglobe/SpinAnalysis'
startDir = os.getenv("PWD")

print "Starting dir: "+startDir

haddLine = "hadd %s/%s/separation.root"%(startDir,baseDirectory)

for dir in subDirs:
  #jobs = [x[1] for x in os.walk(dir)]
  #jobs = [dir+"/"+x for x in jobs[0]]
  jobs = glob.glob(dir+"/job*")
  print "Checking "+dir+" with %d jobs."%len(jobs)
  #os.system("cd "+startDir+"/"+dir)
  os.chdir(startDir+"/"+dir)

  temp = re.findall(r'\d+',dir)
  print "  This one has "+temp[0]+" categories in cosTheta*."

  print "  File to process: "+dir+"/qmu.root"

  processCommand = '%s/bin/diyPlot %s'%(home,startDir+"/"+dir+"/qmu.root --nCats "+temp[0])
  #print processCommand
  os.system(processCommand)
  os.chdir(startDir)

  haddLine += " %s/%s/%sCats_separation.root"%(startDir, dir, temp[0])

#print haddLine
os.system(haddLine)
