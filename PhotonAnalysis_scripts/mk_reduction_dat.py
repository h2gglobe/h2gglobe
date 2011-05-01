#!/usr/bin/env python

from sys import argv

datasets="datasets.txt"
datasetsdir="datasets"
indir=argv.pop(1)
outdir=argv.pop(1)

ds=open(datasets)

for d in ds.read().split("\n"):
    dname = d.strip(" ").lstrip(" ")
    if dname!= "":
        datname="%s/%s.dat" % (datasetsdir, dname)
        print "Making %s" % datname
        f = open(datname ,"w+")
        print >>f, """output=%s/%s.root
        
CaDir=%s/%s typ=1

vtxAlgoParams vertex_selection.dat

inputBranches reduction_input.dat
outputBranches reduction_output.dat
        """ % ( outdir, dname, indir, dname)
        f.close()

ds.close()

