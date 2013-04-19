#!/usr/bin/env python
import sys
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
parser.add_option("--maxCats",dest="maxCats",type="int")
parser.add_option("-d","--directory",dest="directory",type="string",default=".")
parser.add_option("--dryRun",dest="dryRun",action="store_true",default=False)
parser.add_option("--clean",dest="clean",action="store_true",default=False)
(options,args) = parser.parse_args()

home='$CMSSW_BASE/src/h2gglobe/SpinAnalysis'
startDir = os.getenv("PWD")
#os.environ['PWD']

ini_str = '[root]\n' + open(options.datfile, 'r').read()
ini_fp = StringIO.StringIO(ini_str)
config = ConfigParser.RawConfigParser()
config.readfp(ini_fp)

for j in range(options.minCats,options.maxCats+1):
  #print j
  dir=startDir+'/%s/%dCategories'%(options.directory,j)
  os.system('mkdir -p %s'%dir)

  datfile='%s/%dCats_%s'%(dir,j,options.datfile)
  f = open('%s'%(datfile),'w')
  f.write("treefile=%s\n"%config.get('root', 'treefile'))
  f.write("wsfile=%s/Workspace.root\n"%(dir))
  f.write("isMassFac=%d\n"%config.getint('root', 'isMassFac'))
  f.write("globePDFs=%d\n"%config.getint('root', 'globePDFs'))
  f.write("fullSMproc=%d\n"%config.getint('root', 'fullSMproc'))
  f.write("useSMpowheg=%d\n"%config.getint('root', 'useSMpowheg'))
  f.write("useSpin2LP=%d\n"%config.getint('root', 'useSpin2LP'))
  f.write("nBDTCats=%d\n"%config.getint('root', 'nBDTCats'))
  f.write("nSpinCats=%d\n"%j)
  f.close()

  os.chdir(dir)
  processCommand = "%s/bin/diyBuildWorkspace %s"%(home, datfile)
  #print processCommand
  os.system(processCommand)
  os.chdir(startDir)

  processCommand = '%s/subDIY.py -n %d -t %d -f %s -q %s -d %s'%(home, options.njobs, options.toyspj, datfile, options.queue, dir)
  if options.dryRun:
    processCommand += " --dryRun"
  if options.clean:
    processCommand += " --clean"
  #print processCommand
  os.system(processCommand)

