#!/bin/env python

import sys

iname=sys.argv[1]

r9map = { "gold" : ("Gold",(0.94,999)),
          "bad"  : ("Bad",(-999,0.94)),
          }
etamap = { "0_1"   : "LowEta",
           "1_1.5" : "HighEta",
           "1.5_2" : "LowEta",
           "2_3"   : "HighEta",
    }

input=open(iname,'r')
lastcat = None
for line in input.read().split("\n"):
    toks = [ l for l in line.split("\t") if l != "" ]
    if len(toks) == 6:
        cat,first,last,corr,err = toks[0], int(toks[2]), int(toks[3]), float(toks[4]), float(toks[5])
        corr = 1. - corr
        det, eta, r9 = cat.split('-')
        etaRng = eta.split("_")[1:3]
        if etaRng[1] == '2.5':
            etaRng[1] = '3'
        elif etaRng[1] == '1.4442':
            etaRng[1] = '1.5'
        elif etaRng[0] == '1.566':
            etaRng[0] = '1.5'
        
        r9Lab, r9Rng = r9map[r9]
        typ = 0
        etaLab = etamap[ "%s_%s" % ( etaRng[0], etaRng[1] ) ]
        cat = "%s%s%s" % ( det, etaLab, r9Lab )
        if cat != lastcat:
            print
        lastcat = str(cat)
        print "%s %d %s %s %1.2f %1.2f %d %d %1.3g %1.3g" % (cat,  typ, etaRng[0], etaRng[1], r9Rng[0], r9Rng[1], first, last,  corr, err )
        ## print det, etaRng, r9Lab, r9Rng, first,last,corr,err
        

