#!/usr/bin/env python

import sys
import math

def getEffAcc(fname):
    fin=open(fname)
    effAcc = {}
    
    for line in fin.read().split("\n"):
        if "125.0" in line:
            toks = [ t for t in line.split(" ") if t != "" ]
            cat = int(toks[0])
            ggh = float(toks[2])
            qqh = float(toks[3])
            effAcc[cat] = [ggh,qqh]
    fin.close()
    return effAcc


allTunes = {}

for file in sys.argv[1:]:
    tune = file.split("/")[-2].split("_")[-1]
    allTunes[tune] = getEffAcc(file)

tunes = set([ t.replace("UEOFF","") for t in allTunes.keys() ])

ratios = {}

for tune in tunes:
    nominal = allTunes[tune]
    ueoff = allTunes["%sUEOFF" % tune ]
    
    totggh = sum( map( lambda x: x[0], nominal.itervalues() ) )
    totqqh = sum( map( lambda x: x[1], nominal.itervalues() ) )

    totgghUEOFF = sum( map( lambda x: x[0], ueoff.itervalues() ) )
    totqqhUEOFF = sum( map( lambda x: x[1], ueoff.itervalues() ) )

    print tune, "(qqh,qqhUEOFF,ggh,gghUEOFF)"
    print "exA", totqqh, totqqhUEOFF, totggh, totgghUEOFF

    gghcat4 = nominal[4][0]
    qqhcat4 = nominal[4][1]

    gghcat5 = nominal[5][0]
    qqhcat5 = nominal[5][1]

    qqhVBF = qqhcat4 + qqhcat5
    gghVBF = gghcat4 + gghcat5

    gghcat4UEOFF = ueoff[4][0]
    qqhcat4UEOFF = ueoff[4][1]

    gghcat5UEOFF = ueoff[5][0]
    qqhcat5UEOFF = ueoff[5][1]

    qqhVBFUEOFF = qqhcat4UEOFF + qqhcat5UEOFF
    gghVBFUEOFF = gghcat4UEOFF + gghcat5UEOFF

    print "cat4", qqhcat4, qqhcat4UEOFF, gghcat4, gghcat5UEOFF
    print "cat5", qqhcat5, qqhcat5UEOFF, gghcat5, gghcat4UEOFF

    print "ratio", qqhVBFUEOFF / qqhVBF, gghVBFUEOFF / gghVBF
    print "doubleRatio", (qqhVBFUEOFF/totqqhUEOFF) / (qqhVBF/totqqh), (gghVBFUEOFF/totgghUEOFF) / (gghVBF/totggh)
    print "migration",   (qqhcat4UEOFF/qqhVBFUEOFF) / (qqhcat4/qqhVBF), (gghcat4UEOFF / gghVBFUEOFF) / (gghcat4/gghVBF)
    
    qqRatioVBF=qqhVBFUEOFF / qqhVBF
    ggRatioVBF=gghVBFUEOFF / gghVBF
    
    qqDoubleRatioVBF=(qqhVBFUEOFF/totqqhUEOFF) / (qqhVBF/totqqh)
    ggDoubleRatioVBF=(gghVBFUEOFF/totgghUEOFF) / (gghVBF/totggh)

    qqMigration=(qqhcat4UEOFF / qqhVBFUEOFF) / (qqhcat4/qqhVBF)
    ggMigration=(gghcat4UEOFF / gghVBFUEOFF) / (gghcat4/gghVBF)

    qqDoubleRatio4=(qqhcat4UEOFF/totqqhUEOFF) / (qqhcat4/totqqh)
    qqDoubleRatio5=(qqhcat5UEOFF/totqqhUEOFF) / (qqhcat5/totqqh)

    ggDoubleRatio4=(gghcat4UEOFF/totgghUEOFF) / (gghcat4/totggh)
    ggDoubleRatio5=(gghcat5UEOFF/totgghUEOFF) / (gghcat5/totggh)

    ratios[tune] = [ qqRatioVBF, qqDoubleRatioVBF, qqMigration, qqDoubleRatio4, qqDoubleRatio5,
                     ggRatioVBF, ggDoubleRatioVBF, ggMigration, ggDoubleRatio4, ggDoubleRatio5, ]

    print
    
print "uncertainties (ratio, doubleRatio, migration, doubleRatio4, doubleRatio5)"
print "qqh", 
for i in range(10):
    print ("%1.1f%%" % (100.*max( map( lambda x : math.fabs(1.-x[i]), ratios.itervalues() ))) ) ,
    if i == 4:
        print "\nggh",

print "\n"

print "qqh"
for tune in tunes:
    print tune, 
    for i in range(5):
        print ("%1.1f%%" % (100*math.fabs(1.-ratios[tune][i]) ) ),
    print
print "\n"
    
print "ggh"
for tune in tunes:
    print tune, 
    for i in range(5):
        print ("%1.1f%%" % (100*math.fabs(1.-ratios[tune][i+5]) ) ),
    print
print "\n"
