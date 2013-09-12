#!/usr/bin/env python

import sys, glob, os, commands

### def usage():
###     print "Usage: check_fitter.py <task_folder> [post_processing_script]"
###     print ""
###     print " The script checks the status of all fitter jobs in <task_folder>"
###     print "  Upon completion of all jobs, all workspaces are combined and the script returrns with 0."
###     print "  Otherwise the script returns with -1."
###     print ""
###     print " A continuos monitoring of the folder can be performed with the following command line:"
###     print "   while( ! ./check_fitter.py <task_folder> [post_processing_script] ); do sleep 360; done"
###     print ""
###     print " Optionally, a post-processing script can be called after the workspaces have been merged."
###     print " The script is rung with <taks_folder> as argument."
###     print "  See do_intepolation.sh for an example."
###     print ""
    

if len(sys.argv) == 1:
    usage()
    sys.exit(-1)


taskdir = sys.argv[1]

samples = glob.glob( "%s/*.dat" % taskdir )

for sample in samples:
    status = {}
    jobs = glob.glob( "%s_[0-9]*.[a-z]*" % sample )
    for j in jobs:
        toks = j.split("_")
        id,stat = toks[-1].split(".")
        if stat in ["sub","run","fail","done"]:
            status[id] = stat


    groups = { "sub" : [], "run" : [], "fail" : [], "done" : [] } 
    for i,s in status.iteritems():
        groups[s].append(i)

    print sample
    for g,jo in groups.iteritems():
        if len(jo) > 0:
            print "%d jobs %s" % (len(jo), g),
            print "[",
            for j in jo: print j,
            print "]"
    print
    
