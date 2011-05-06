#!/usr/bin/env python

from sys import argv
from os import mkdir, system

def mk_outdir( name ):
    if name.startswith("rfio:") or name.startswith("/castor"):
        system("rfmkdir %s" % name.lstrip("rfio:") )
    else:
        try:
            mkdir( outdir )
        except Exception, e:
            print e


datasets="datasets.txt"
indir=argv.pop(1)
outdir=argv.pop(1)

if len(argv) > 1:
    datasets=argv.pop(1)

datasetsdir=datasets.rsplit(".",1)[0]
try:
    mkdir( datasetsdir )
except Exception, e:
    print e

ds=open(datasets)

mk_outdir( outdir )
for d in ds.read().split("\n"):
    s = d.lstrip(" ").lstrip(" ")
    if s == "" or s.startswith("#"):
        continue

    sl = [ t for t in d.split(" ") if t != "" ]
    dname, dtype = sl
    datname="%s/%s.dat" % (datasetsdir, dname)
    print "Making configuration for %s (type %s)" % ( datname, dtype )
    f = open(datname ,"w+")
    print >>f, """output=%s/%s/%s.root
        
CaDir=%s/%s typ=%s

analyzer PhotonAnalysis photonanalysis.dat

inputBranches reduction_input.dat
outputBranches reduction_output.dat
        """ % ( outdir, dname, dname, indir, dname, dtype)
    f.close()

    mk_outdir( "%s/%s" % ( outdir, dname )  )

ds.close()

