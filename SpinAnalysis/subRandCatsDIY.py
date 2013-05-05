#!/usr/bin/env python
import sys
import json
import random

from subCatsCommon import *

if not options.catToys:
  print "No catToys defined"
  exit()

if not options.minDist:
  print "No minDist defined: Using default (0.01)"
  options.minDist = 0.01

if not options.minCats:
  print "No minCats defined: Using default (4)"
  options.minCats = 4
else:
  if options.minCats < 2:
    print "minCats option is too low, setting it to 2"
    options.minCats = 2

if not options.maxCats:
  print "No maxCats defined: Using default (7)"
  options.maxCats = 7

config = readDatfile(options.datfile)

print " "

catInfo = {}
if not options.new:
  if os.path.exists("%s/catInfo.json"%options.directory):
    fp = open("%s/catInfo.json"%options.directory)
    catInfo = json.load(fp)

print "Building directory structure and config files"
for j in range(options.minCats,options.maxCats+1):
  catName = "%dCategories"%j
  catDir = "%s/%s/%s"%(startDir,options.directory,catName)
  print "   %s"%catName

  jobs = {}
  if options.new:
    command = "rm -R %s"%catDir
    os.system(command)
    command = "mkdir -p %s"%catDir
    os.system(command)
  else:
    if os.path.exists('%s/jobs.json'%catDir):
      fp = open('%s/jobs.json'%catDir)
      jobs = json.load(fp)

  previousJobs = len(jobs)
  for job in range(previousJobs, previousJobs+options.catToys):
    jobname = "Boundaries%d"%job
    jobdir = "%s/%s"%(catDir,jobname)
    print "\tJob: %s"%jobname

    boundaries = []
    rejected = True
    while rejected:
      rejected = False
      boundaries = sorted([random.random() for i in range(j-1)])
      print "\t  Generated boundaries: [%s]"%', '.join(map(str, boundaries))

      if boundaries[0] < options.minDist or 1-boundaries[j-2] < options.minDist:
        rejected = True
      for i in range(j-2):
        if boundaries[i+1]-boundaries[i] < options.minDist:
          rejected = True
      if rejected:
        print "\t   Rejected boundaries for not respecting minimum distance between boundaries"

      if not rejected:
        for testJob in jobs:
          #distances between boundaries
          dists = [boundaries[i]-jobs[testJob][2][i] for i in range(j-1)]
          dists = [i/options.minDist for i in dists]
          dists = [i*i for i in dists]
          distance = sum(dists)/(j-1)
          if distance < 1.5:
            rejected = True
            print "\t   Rejected boundaries for being too similar to a previous boundary"
            break

    bounds = "1 "
    for i in boundaries:
      bounds += "%f "%(1-i)
    bounds += "0"

    os.system('mkdir -p %s'%jobdir)

    datfile='%s/%dCats_%s'%(jobdir,j,options.datfile)
    f = open('%s'%(datfile),'w')
    copyDatfile(config, f, ["wsfile", "nSpinCats", "catBoundaries"])
    f.write("wsfile=%s/Workspace.root\n"%(jobdir))
    f.write("nSpinCats=%d\n"%j)
    f.write("catBoundaries=%s\n"%bounds)
    f.close()

    prepareJobArray(jobname, jobdir, datfile, options.njobs, options.toyspj)

    jobs[job] = (jobname, jobdir, boundaries, datfile, options.njobs, options.toyspj, options.queue)

  with open('%s/jobs.json'%catDir, 'wb') as fp:
    json.dump(jobs, fp)

  catInfo[str(j)] = (catName, catDir)
  print " "

with open("%s/catInfo.json"%options.directory, 'wb') as fp:
  json.dump(catInfo, fp)


print " "
print "Building workspaces"
for cat in catInfo:
  print "   %s"%catInfo[cat][0]
  fp = open('%s/jobs.json'%catInfo[cat][1])
  jobs = json.load(fp)

  processedJobs = []
  if os.path.exists('%s/processedJobs.json'%catInfo[cat][1]):
    fp = open('%s/processedJobs.json'%catInfo[cat][1])
    processedJobs = json.load(fp)

  for job in jobs:
    if job not in processedJobs:
      print "\tJob: %s"%jobs[job][0]

      os.chdir(jobs[job][1])
      processCommand = "%s/bin/diyBuildWorkspace %s &> buildWorkspace.log"%(startDir, jobs[job][3])
      if options.dryRun: print processCommand
      else:
        os.system(processCommand)
        processedJobs.append(job)
      os.chdir(startDir)

  with open('%s/processedJobs.json'%catInfo[cat][1], 'wb') as fp:
    json.dump(processedJobs, fp)

  print " "


print " "
print "Submitting jobs"
for cat in catInfo:
  print "   %s"%catInfo[cat][0]
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

  for job in jobs:
    if job in processedJobs and job not in submittedJobs:
      print "\tJob: %s"%jobs[job][0]
      submitJobArray(jobs[job][0], jobs[job][1], jobs[job][4], jobs[job][6])
      if not options.dryRun:
        submittedJobs.append(job)

  with open('%s/submittedJobs.json'%catInfo[cat][1], 'wb') as fp:
    json.dump(submittedJobs, fp)
