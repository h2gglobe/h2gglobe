#!/usr/bin/env python

import sys, glob, os, commands
from mk_fitter import Conf 

taskdir = sys.argv[1]
docombine=None
if( len(sys.argv) > 2 ):
    docombine=sys.argv[2]

jobs = glob.glob( "%s/sub*.sh" % taskdir )
status = {}

for j in jobs:
    status[j] = ""
    for stat in "run","fail","done":
        if os.path.isfile("%s.%s" % (j,stat) ):
            status[j] = stat

groups = { "" : [], "run" : [], "fail" : [], "done" : [] } 
for i,s in status.iteritems():
    groups[s].append(i)

for g,jo in groups.iteritems():
    print
    print "%d jobs in status %s" % (len(jo), g)
    print "     ",
    for j in jo: print "%s" % j,
    print


if len(groups["done"]) == len(jobs):
    print "All jobs completed"
    filestocmb=glob.glob("%s/filestocombine_*.dat" % taskdir)
    if len(filestocmb) == 1:
        cfg = Conf() 
        cfile = open( filestocmb[0], "r" )

        for line in cfile.read().split("\n"):
            if "histfile" in line:
		cfg.read_histfile(line)
		line = line.replace(cfg.histdir,"./")
                break


        if not os.path.isfile( "%s.%s" % ( os.path.join(cfg.histdir,cfg.histfile[0]), cfg.histfile[1] ) ):
            ## sys.argv = [ "combiner.py", "-i", filestocmb[0] ]
            ## import combiner
            os.system("python combiner.py -i %s" % filestocmb[0] )
            os.system("%s %s | tee -a %s/do_combine.log" % ( docombine, taskdir, taskdir ) )

    sys.exit(0)
            
            
sys.exit(-1)
