#!/bin/env python

from optparse import OptionParser
from pprint import pprint
import os.path

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    
    return None

parser = OptionParser(usage="usage: %prog [options] some_data.json\nrun with --help to get list of options")
parser.add_option("--skipLumi", dest="skipLumi", action="store_true", default=False, help="Skip lumi calculation [default: %default].")
parser.add_option("--skipPileup", dest="skipPu", action="store_true", default=False, help="Skip PU calculation [default: %default].")
parser.add_option("--minBiasXsec", dest="mbXsec",  default=68000,  type="int", help="Minimum bias cross section to use in calculating the PU [default: %default].")
parser.add_option("--maxPileupBin", dest="maxPuBin",  default=100,  type="int", help="Largest PU value [default: %default].")
parser.add_option("--numPileupBins", dest="nPuBin",  default=100,  type="int", help="Number of PU bins [default: %default].")
parser.add_option("--puJson", dest="puJson",  default='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_JSON_DCSONLY_190389-193575.txt',  type="string", help="Pileup JSON file [default: %default].")

(options, args) = parser.parse_args()

if options.skipPu and options.skipLumi:
    raise RuntimeError("Skip PU and Lumi calculation? Nothing to do.")


if len(args) == 0:
    parser.print_usage()
    import sys
    sys.exit(1)
options.inputJson = args[0]
options.inputJsonName = os.path.basename(options.inputJson)

jobs = list()

#pu
if not options.skipPu:
    calcModes = [ 'true', 'observed' ]
    cmd = "pileupCalc.py\
    --calcMode %(calcMode)s   --minBiasXsec %(mbXsec)s\
    -i %(inputJson)s   --inputLumiJSON %(puJson)s\
    --maxPileupBin %(maxPuBin)d   --numPileupBins %(nPuBin)d\
    %(inputJsonName)s.%(mbXsec)s.%(calcMode)s.pileup.root &> %(inputJsonName)s.%(mbXsec)s.%(calcMode)s.pileup.root.log"
    
    for options.calcMode in calcModes:
        print "Queuing %(calcMode)s PU." % vars(options)
        jobs.append(cmd % vars(options))


#lumi
if not options.skipLumi:
    calculators = [ 'lumiCalc2', 'pixelLumiCalc' ]
    cmd = "%(calculator)s.py -i %(inputJson)s --nowarning overview &> %(inputJsonName)s.%(calculator)s.lumi"

    for options.calculator in calculators:
        program = "%(calculator)s.py" % vars(options)
        if not which(program):
            raise RuntimeError("%s not found. Check https://twiki.cern.ch/twiki/bin/view/CMS/LumiCalc for how to install lumi stuff." % program)
        print "Queuing %(calculator)s lumi." % vars(options)
        jobs.append(cmd % vars(options))

#function to be used in parallel
def runCmd(cmd):
    if not cmd: return 0
    from subprocess import call
    print cmd
    return call(cmd, shell=True)

from multiprocessing import Pool
pool = Pool()

# run all jobs in parallel
ret = pool.map(runCmd, jobs)
if reduce(lambda x, y: x+y, ret):
        raise RuntimeError, "Non-zero return code in parallel loop. Check logs."

