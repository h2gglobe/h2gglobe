#!/usr/bin/env python

import os
import re
import sys
import glob
import fnmatch

baseDirectory = sys.argv[1]
subDirs = [x[1] for x in os.walk(baseDirectory)]
subDirs = [baseDirectory+"/"+x for x in subDirs[0]]

startDir = os.getenv("PWD")

home='$CMSSW_BASE/src/h2gglobe/SpinAnalysis'

for dir in subDirs:
  #jobs = [x[1] for x in os.walk(dir)]
  #jobs = [dir+"/"+x for x in jobs[0]]
  jobs = glob.glob(dir+"/job*")
  print "Checking "+dir+" with %d job* directories."%len(jobs)

  numCats = int(re.findall('\d+', dir.split("/")[-1])[0])

  running = []
  finished = []
  failed = []
  not_running = []

  for job in jobs:
    jobnum = int(re.findall('\d+', job.split("/")[-1])[0])
    if os.path.exists(dir+"/job%d"%jobnum):
      if os.path.exists(dir+"/job%d/sub%d.sh.done"%(jobnum, jobnum)):
        finished.append(jobnum)
      else:
        if os.path.exists(dir+"/job%d/sub%d.sh.fail"%(jobnum, jobnum)) or (os.path.exists(dir+"/job%d/sub%d.sh.log"%(jobnum, jobnum))):
          failed.append(jobnum)
        else:
          if os.path.exists(dir+"/job%d/sub%d.sh.run"%(jobnum, jobnum)):
            running.append(jobnum)
          else:
            not_running.append(jobnum)

  numRuns = len(running) + len(finished) + len(failed) + len(not_running)
  print "\tFound %d jobs, %d still running, %d finished, %d failed, %d in unknown status (probably not running or even ran)"%(numRuns, len(running), len(finished), len(failed), len(not_running))

  joblist = ""
  if len(failed) != 0:
    var = raw_input('\tDo you want to resubmit the failed jobs? [y]/[n]\n\t  ')
    if var=='y' or var=='Y' or var =='yes':
      for j in failed:
        joblist+="%d,"%j

  if len(not_running) != 0:
    var = raw_input('\tDo you want to submit the jobs with unknown status? [y]/[n]\n\t  ')
    if var=='y' or var=='Y' or var =='yes':
      for j in not_running:
        joblist+="%d,"%j

  if joblist != "":
    joblist = joblist[0:-1]
    print "\tJobs to be submitted: ", joblist
    var = raw_input("\t  To which queue do you want to submit the jobs? [default: 1nh]\n\t    ")
    if var == "":
      var = "1nh"
    os.chdir(dir)
    jobListName = "%dCategories[%s]"%(numCats,joblist)
    logFile = "%s/job%%I/sub%%I.sh.log"%(os.getcwd())
    processCommand = 'bsub -J "%s" -q %s -o "%s" %d_Cats.sh'%(jobListName,var,logFile,numCats)
    #print processCommand
    os.system(processCommand)
    os.chdir(startDir)

  if len(finished) != 0:
    print "\tFinished jobs: ",finished
    var = raw_input("\t  Do you want to hadd the finished jobs? [y]/[n]\n\t    ")
    if var == "":
      var = "y"
    if var=='y' or var=='Y' or var =='yes':
      haddline="hadd -f %s/qmu.root"%dir
      for jobnum in finished:
        if not os.path.exists(dir+"/job%d/LLOut.root"%jobnum):
          print "\t  Unable to find LLOut.root for job %d"%jobnum
        else:
          haddline += " %s/job%d/LLOut.root"%(dir, jobnum)
      #print haddline
      os.system(haddline)


  print ""
  #processCommand = '%s/checkDIY.py %d %s'%(home,len(jobs),dir)
  #print processCommand
  #os.system(processCommand)

"""
foo = "3"
bar = 1
bar += int(foo)
print bar"""
