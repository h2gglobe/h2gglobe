#!/usr/bin/env python
import os
import StringIO
import ConfigParser

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-n","--njobs",dest="njobs",type="int",)
parser.add_option("-t","--toysperjob",dest="toyspj",type="int")
parser.add_option("-f","--datfile",dest="datfile",type="string")
parser.add_option("-q","--queue",dest="queue",type="string")
parser.add_option("--minCats",dest="minCats",type="int")
parser.add_option("--minDist",dest="minDist",type="float")
parser.add_option("--maxCats",dest="maxCats",type="int")
parser.add_option("--catToys",dest="catToys",type="int")
parser.add_option("-d","--directory",dest="directory",type="string",default="./testDir")
parser.add_option("--dryRun",dest="dryRun",action="store_true",default=False)
parser.add_option("--clean",dest="clean",action="store_true",default=False)
parser.add_option("--new",dest="new",action="store_true",default=False)
parser.add_option("-p","--pickle",dest="pickleFile",type="string")
parser.add_option("-j","--json",dest="jsonFile",type="string")
(options,args) = parser.parse_args()

if not options.datfile:
  print "No datfile defined"
  exit()

if not options.directory:
  print "No directory defined"
  exit()

if not options.njobs:
  print "No number of jobs defined: Using default (50)"
  options.njobs = 50

if not options.toyspj:
  print "No toys per job defined: Using default (40)"
  options.toyspj = 40

if not options.queue:
  print "No queue defined: Using default queue (1nh)"
  options.queue = "1nh"

home='$CMSSW_BASE/src/h2gglobe/SpinAnalysis'
startDir = os.getenv("PWD")
#startDir = os.environ['PWD']
#startDir = os.getcwd()

def readDatfile(file):
  ini_str = '[root]\n' + open(file, 'r').read()
  ini_fp = StringIO.StringIO(ini_str)
  config = ConfigParser.RawConfigParser()
  config.optionxform=str
  config.readfp(ini_fp)
  return config

def copyDatfile(config, file, exclude=[]):
  for option in config.options('root'):
    if option not in exclude:
      file.write("%s=%s\n"%(option,config.get('root', option)))

def prepareJobArray(jobname, dir, datfile, njobs, ntoys):
  os.system('cp %s/bin/diySeparation %s/'%(os.getcwd(),dir))

  #Prepare the job directories
  for i in range(1,njobs+1):
    jobdir = "%s/job%d"%(dir,i)

    os.system("mkdir -p %s"%jobdir)

    os.system("rm -f %s/*.sh.fail"%jobdir)
    os.system("rm -f %s/*.sh.done"%jobdir)
    os.system("rm -f %s/*.sh.log"%jobdir)
    os.system("rm -f %s/*.sh.run"%jobdir)

  #Make the script to run
  f = open('%s/%s.sh'%(dir,jobname),'w')
  f.write('#!/bin/bash\n')
  f.write('cd %s/job$LSB_JOBINDEX\n'%(dir))
  f.write('touch sub$LSB_JOBINDEX.sh.run\n')
  f.write('eval `scramv1 runtime -sh`\n')
  subline = '../diySeparation %s %d'%(datfile,ntoys)
  f.write('if ( %s )\n'%subline)
  f.write('\tthen touch sub$LSB_JOBINDEX.sh.done\n')
  f.write('\telse touch sub$LSB_JOBINDEX.sh.fail\n')
  f.write('fi\n')
  f.write('rm -f sub$LSB_JOBINDEX.sh.run\n')
  os.system('chmod +x %s'%f.name)

def submitJobArray(jobname, dir, njobs, queue):
  os.chdir(dir)
  jobListName = "%s[1-%d]"%(jobname,njobs)
  logFile = "%s/job%%I/sub%%I.sh.log"%(os.getcwd())
  processCommand = 'bsub -J "%s" -q %s -o "%s" %s.sh'%(jobListName,queue,logFile,jobname)
  if not options.dryRun: os.system(processCommand)
  else: print processCommand
  os.chdir(startDir)

def prepareAndSubmitJobArray(jobname, dir, datfile, njobs, ntoys, queue):
  prepareJobArray(jobname, dir, datfile, njobs, ntoys)
  submitJobArray(jobname, dir, njobs, queue)
