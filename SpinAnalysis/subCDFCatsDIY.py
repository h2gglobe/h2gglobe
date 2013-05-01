#!/usr/bin/env python
import sys
import pickle

from subCatsCommon import *

if not options.pickleFile:
  print "No pickle file defined"
  exit()

if not options.minCats:
  print "No minCats defined: Using lowest binning in pickle"

if not options.maxCats:
  print "No maxCats defined: Using highest binning in pickle"

config = readDatfile(options.datfile)

pkl_file = open(options.pickleFile, 'rb')
data1 = pickle.load(pkl_file)

print "Building directory structure and config files"
for key, entry in data1.items():
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
  for i in entry["quantiles"]:
    boundaries += "%f "%(1-i)
  print "\t  Boundaries: "+boundaries

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
