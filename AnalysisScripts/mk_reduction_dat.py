#!/usr/bin/env python

from sys import argv
from os import mkdir, system, path

def mk_outdir( name ):
    if name.startswith("rfio:") or name.startswith("/castor"):
        system("rfmkdir %s" % name.lstrip("rfio:") )
    else:
        try:
            mkdir( outdir )
        except Exception, e:
            pass


datasets="datasets.txt"
indir=argv.pop(1)
outdir=argv.pop(1)

if indir != "-":
    indir = "%s/" % indir
else:
    indir = ""

if len(argv) > 1:
    datasets=argv.pop(1)

datasetsdir=datasets.rsplit(".",1)[0]
try:
    mkdir( datasetsdir )
except Exception, e:
    pass

ds=open(datasets)

mk_outdir( outdir )
for d in ds.read().split("\n"):
    s = d.replace("\t","").lstrip(" ").lstrip(" ")
    if s == "" or s.startswith("#"):
        continue

    props = ""
    sl = [ t for t in s.split(" ") if t != "" ]
    dname = sl.pop(0)
    dtype = sl.pop(0)
    if ":" in dname:
        iname,dname = dname.rsplit(":",1)
    else:
        iname = dname
        dname = dname.replace("_AODSIM","").replace("_AOD","")
    analyzer = "analyzer PhotonAnalysis photonanalysis.dat"
    getanalyzer = False
    for s in sl:
        if "append" in s:
            dname += s.rsplit("=",1)[1]
        elif "analyzer" in s:
            analyzer = "%s" % s
            getanalyzer = True
        elif getanalyzer:
            analyzer += " %s" % s
        else:
            print s
            props += " %s" % s
    basname = path.basename(dname)
    datname="%s/%s.dat" % (datasetsdir, basname )
    print "Making configuration for %s (type %s)" % ( datname, dtype )
    f = open(datname ,"w+")
    print >>f, """output=%s/%s/%s.root
    
CaDir=%s%s typ=%s %s

%s

inputBranches reduction_input.dat
outputBranches reduction_output.dat
        """ % ( outdir, basname, basname, indir, iname, dtype, props, analyzer)
    f.close()

    mk_outdir( "%s/%s" % ( outdir, dname )  )

ds.close()

