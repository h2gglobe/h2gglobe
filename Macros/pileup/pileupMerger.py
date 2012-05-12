#!/bin/env python

from optparse import OptionParser
from pprint import pprint
import os.path as p
from subprocess import check_call as call

parser = OptionParser(usage="usage: %prog [options] EOS_source_directory\nrun with --help to get list of options")
parser.add_option("--putBack", dest="putBack", action="store_true", default=False, help="Put back merged file in source directory [default: %default].")

(options, args) = parser.parse_args()

if len(args) == 0:
    parser.print_usage()
    import sys
    sys.exit(1)
options.inDir = args[0]
options.inDirName = p.basename(options.inDir)

cmd = """cmsLs %(inDir)s | awk '{ print $5 }' | xargs cmsPfn | sed 's/\?.*$//' > %(inDirName)s.files.txt
hadd -T %(inDirName)s.pileup.root @%(inDirName)s.files.txt &> %(inDirName)s.pileup.root.log""" 

if not options.inDirName:
    raise RuntimeError("Empty end directory name (which defined the sample name). Check path.")

print cmd % vars(options)
call(cmd % vars(options), shell=True)

if options.putBack:
    cmd = """xrdcp %(inDirName)s.pileup.root `cmsPfn %(inDir)s | sed 's/\?.*$//'`/%(inDirName)s.pileup.root"""
    print "Copying back root"
    call(cmd % vars(options), shell=True)
    cmd = """xrdcp %(inDirName)s.pileup.root.log `cmsPfn %(inDir)s | sed 's/\?.*$//'`/%(inDirName)s.pileup.root.log"""
    print "Copying back log"
    call(cmd % vars(options), shell=True)

