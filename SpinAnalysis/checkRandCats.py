#!/usr/bin/env python

import os
import re
import sys
import glob
import json
import fnmatch

baseDirectory = sys.argv[1]
subDirs = [x[1] for x in os.walk(baseDirectory)]
subDirs = [baseDirectory+"/"+x for x in subDirs[0]]

startDir = os.getenv("PWD")

catInfo = {}
if os.path.exists("%s/catInfo.json"%baseDirectory):
  fp = open("%s/catInfo.json"%baseDirectory)
  catInfo = json.load(fp)

for cat in catInfo:
  print "   Checking %s"%catInfo[cat][0]
  fp = open('%s/jobs.json'%catInfo[cat][1])
  jobs = json.load(fp)

  processedJobs = []
  if os.path.exists('%s/processedJobs.json'%catInfo[cat][1]):
    fp = open('%s/processedJobs.json'%catInfo[cat][1])
    processedJobs = json.load(fp)

  submittedJobs = []
  if os.path.exists('%s/submittedJobs.json'%catInfo[cat][1]):
    fp = open('%s/submittedJobs.json'%catInfo[cat][1])
    submittedJobs = json.load(fp)

  haddedJobs = []
  if os.path.exists('%s/haddedJobs.json'%catInfo[cat][1]):
    fp = open('%s/haddedJobs.json'%catInfo[cat][1])
    haddedJobs = json.load(fp)

  print "\tOut of %d job arrays, %d have been processed and %d have been submitted."%(len(jobs), len(processedJobs), len(submittedJobs))

  for job in jobs:
    if job in processedJobs and job in submittedJobs:
      print "\tJob Array: %s"%jobs[job][0]
      jobdir = jobs[job][1]

      numRuns = jobs[job][4]
      running = []
      finished = []
      failed = []
      not_running = []

      for index in range(1, numRuns+1):
        if os.path.exists("%s/job%d"%(jobdir, index)):
          if os.path.exists("%s/job%d/sub%d.sh.done"%(jobdir, index, index)):
            finished.append(index)
          else:
            if os.path.exists("%s/job%d/sub%d.sh.fail"%(jobdir, index, index)) or (os.path.exists("%s/job%d/sub%d.sh.log"%(jobdir, index, index))):
              failed.append(index)
            else:
              if os.path.exists("%s/job%d/sub%d.sh.run"%(jobdir, index, index)):
                running.append(index)
              else:
                not_running.append(index)
        else:
          print "Couldn't find job %d"%index

      print "\t  %d jobs: %d still running, %d finished, %d failed, %d in unknown status (probably not running or even ran)"%(numRuns, len(running), len(finished), len(failed), len(not_running))

      joblist = ""
      if len(failed) != 0:
        var = raw_input('\t    Do you want to resubmit the failed jobs? [y]/[n]\n\t      ')
        if var=='y' or var=='Y' or var =='yes':
          for j in failed:
            joblist+="%d,"%j

      if len(not_running) != 0:
        var = raw_input('\t    Do you want to submit the jobs with unknown status? [y]/[n]\n\t      ')
        if var=='y' or var=='Y' or var =='yes':
          for j in not_running:
            joblist+="%d,"%j

      if joblist != "":
        joblist = joblist[0:-1]
        print "\t    Jobs to be submitted: ", joblist

        os.chdir(jobdir)
        jobListName = "%s[%s]"%(jobs[job][0],joblist)
        logFile = "%s/job%%I/sub%%I.sh.log"%(jobdir)
        processCommand = 'bsub -J "%s" -q %s -o "%s" %s.sh'%(jobListName,jobs[job][6],logFile,jobs[job][0])
        #print processCommand
        os.system(processCommand)
        os.chdir(startDir)

      if len(finished) != 0 and job not in haddedJobs:
        print "\t    Finished jobs: ",finished
        var = raw_input("\t      Do you want to hadd the finished jobs? [y]/[n]\n\t        ")
        if var == "":
          var = "y"
        if var=='y' or var=='Y' or var =='yes':
          haddline="hadd -f %s/qmu.root"%jobdir
          for jobnum in finished:
            if not os.path.exists("%s/job%d/LLOut.root"%(jobdir,jobnum)):
              print "\t  Unable to find LLOut.root for job %d"%jobnum
            else:
              haddline += " %s/job%d/LLOut.root"%(jobdir, jobnum)
          #print haddline
          os.system(haddline)
          if len(finished) == numRuns:
            haddedJobs.append(job)

  with open('%s/haddedJobs.json'%catInfo[cat][1], 'wb') as fp:
    json.dump(haddedJobs, fp)

  print ""
