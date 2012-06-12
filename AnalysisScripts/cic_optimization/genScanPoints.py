#!/bin/env python

import itertools

#### enum phoCiCIDLevel { phoNOCUTS=0, phoLOOSE, phoMEDIUM, phoTIGHT, phoSUPERTIGHT, phoHYPERTIGHT1, phoHYPERTIGHT2,
####                      phoHYPERTIGHT3, phoHYPERTIHT4, phoNCUTLEVELS };
toScan = [
    ### [4,5,6,7],
    ### [4,5,6,7],
    ### [4,5,6,7],
    ### [4,5,6,7]
    ]

simple = [4,5,6,7]

combinations = [ (s,s,s,s) for s in simple  ]
if len(toScan) != 0:
    combinations+=list(itertools.product(*toScan))

print "Number of combinations: %d " % ( len(combinations) )

print combinations

for ic in range(len(combinations)):
    f = open("point%d.dat" % ( ic ), "w+" )
    f.write( "cicCutLevels=%d,%d,%d,%d\n" % combinations[ic] )
    f.close()

