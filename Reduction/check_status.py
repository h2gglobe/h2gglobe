#!/usr/bin/env python

import sys, glob, os, commands

if len(sys.argv) == 1:
    usage()
    sys.exit(-1)


taskdir = sys.argv[1]

samples = glob.glob( "%s/*.dat" % taskdir )

for sample in samples:
    status = {}
    tot = 0
    jobs = glob.glob( "%s_[0-9]*.[a-z]*" % sample )
    for j in jobs:
        toks = j.split("_")
        id,stat = toks[-1].split(".")
        if stat in ["sub","run","fail","done"]:
            status[id] = stat

    groups = { "sub" : [], "run" : [], "fail" : [], "done" : [] } 
    for i,s in status.iteritems():        
        groups[s].append(i)
        tot += 1

    print sample
    name = os.path.basename(sample).rsplit(".",1)[0]
    for g,jo in groups.iteritems():
        if len(jo) > 0:
            print "%d jobs %s" % (len(jo), g),
            print "[",
            for j in jo:
                print j,
            print "]"
            if g == "fail":
                print "find %s " % name,
                print "-name %s_%s.root" % (name,j[0]) ,
                for j in jo[1:]:
                    print "-or -name %s_%s.root" % (name,j) ,
                print " -exec rm -v \\{\\} \\; ;"
                print "./submit_reduction.sh %s %s %d" % ( taskdir, name, tot ),
                for j in jo:
                    print j,
                print
    print
    
