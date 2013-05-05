#!/usr/bin/env python
import sys
import json

from subCatsCommon import *

if not options.jsonFile:
  print "No JSON file defined"
  exit()

if not options.minCats:
  print "No minCats defined: Using lowest binning in json"

if not options.maxCats:
  print "No maxCats defined: Using highest binning in json"

config = readDatfile(options.datfile)

json_file = open(options.jsonFile)
data1 = json.load(json_file)

print "Building directory structure and config files"
for key, entry in data1.items():
  key = int(key)
  if options.minCats:
    if key < options.minCats:
      continue
  if options.maxCats:
    if key > options.maxCats:
      continue
  jobname = "%dCategories"%(key)
  jobdir = "%s/%s/%s"%(startDir,options.directory,jobname)
  print "\tJob: %s"%jobname

  boundaries=""
  bounds = entry["boundaries"]
  bounds = [min(max(x,0),1) for x in bounds]
  bounds = sorted(set(bounds))
  for i in bounds:
    boundaries += "%f "%(1-i)
  print "\t  %d boundaries: %s"%(len(bounds),boundaries)

  os.system('mkdir -p %s'%jobdir)

  datfile='%s/%dCats_%s'%(jobdir,key,options.datfile)
  f = open('%s'%(datfile),'w')
  copyDatfile(config, f, ["wsfile", "nSpinCats", "catBoundaries"])
  f.write("wsfile=%s/Workspace.root\n"%(jobdir))
  f.write("nSpinCats=%d\n"%key)
  f.write("catBoundaries=%s\n"%boundaries)
  f.close()

  prepareJobArray(jobname, jobdir, datfile, options.njobs, options.toyspj)

print "Building workspaces"
for key, entry in data1.items():
  key = int(key)
  if options.minCats:
    if key < options.minCats:
      continue
  if options.maxCats:
    if key > options.maxCats:
      continue
  jobname = "%dCategories"%(key)
  jobdir = "%s/%s/%s"%(startDir,options.directory,jobname)
  datfile='%s/%dCats_%s'%(jobdir,key,options.datfile)
  print "\tJob: %s"%jobname

  os.chdir(jobdir)
  processCommand = "%s/bin/diyBuildWorkspace %s &> buildWorkspace.log"%(startDir, datfile)
  if options.dryRun: print processCommand
  else: os.system(processCommand)
  os.chdir(startDir)

print "Submitting jobs"
for key, entry in data1.items():
  key = int(key)
  if options.minCats:
    if key < options.minCats:
      continue
  if options.maxCats:
    if key > options.maxCats:
      continue
  jobname = "%dCategories"%(key)
  jobdir = "%s/%s/%s"%(startDir,options.directory,jobname)
  datfile='%s/%dCats_%s'%(jobdir,key,options.datfile)
  print "\tJob: %s"%jobname

  submitJobArray(jobname, jobdir, options.njobs, options.queue)
