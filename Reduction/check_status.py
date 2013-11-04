#!/usr/bin/env python

import sys, glob, os, commands

if len(sys.argv) == 1:
    usage()
    sys.exit(-1)

for taskdir in sys.argv[1:]:
    samples = glob.glob( "%s/*.dat" % taskdir )

    for sample in samples:
        status = {}
        tot = 0
        if os.path.exists("%s.njobs" % sample):
            nj = open("%s.njobs" % sample)
            tot = int(nj.read())
            nj.close()
        jobs = glob.glob( "%s_[0-9]*.[a-z]*" % sample )
        
        for j in jobs:
            toks = j.split("_")
            id,stat = toks[-1].split(".")
            id = int(id)
            if stat in ["sub","run","fail","done"]:
                status[id] = stat

        groups = { "sub" : [], "run" : [], "fail" : [], "done" : [] } 
        for i,s in status.iteritems():
            groups[s].append(i)


        print "%s %d %1.0f%% " % ( sample, tot, float(len(groups["done"]))/float(tot)*100. )
        name = os.path.basename(sample).rsplit(".",1)[0]
        for g,jo in groups.iteritems():
            jo = sorted(jo)
            if len(jo) > 0:
                print "%d jobs %s" % (len(jo), g),
                print "[",
                for j in jo:
                    print j,
                print "]"
                if g == "fail":
                    print "find %s " % name,
                    print "-name %s_%d.root" % (name,jo[0]) ,
                    for j in jo[1:]:
                        print "-or -name %s_%s.root" % (name,j) ,
                    print " -exec rm -v \\{\\} \\; ;"
                    print "./submit_reduction.sh %s %s %d" % ( taskdir, name, tot ),
                    for j in jo:
                        print j,
                    print                
        print
    
