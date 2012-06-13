#!/bin/env python

import itertools

#### enum phoCiCIDLevel { phoNOCUTS=0, phoLOOSE, phoMEDIUM, phoTIGHT, phoSUPERTIGHT, phoHYPERTIGHT1, phoHYPERTIGHT2,
####                      phoHYPERTIGHT3, phoHYPERTIHT4, phoNCUTLEVELS };
toScan = [
    [4,5,6],
    [4,5,6],
    [4,5,6],
    [4,5,6]
    ]

simple = []
## simple = [4,5,6]

allcombinations = [ (s,s,s,s) for s in simple  ]
if len(toScan) != 0:
    allcombinations+=list(itertools.product(*toScan))

combinations = [ c for c in allcombinations if c[0]<=c[2] and c[1]<=c[3] and c[0]<=c[1] and c[2]<=c[3] ]

print "Number of combinations: %d " % ( len(combinations) )

print combinations

for ic in range(len(combinations)):
    f = open("cic_optimization/point%d.dat" % ( ic ), "w+" )
    f.write( "cicCutLevels=%d,%d,%d,%d\n" % combinations[ic] )
    f.close()

