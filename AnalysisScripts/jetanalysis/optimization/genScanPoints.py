#!/bin/env python

import itertools
from optparse import OptionParser, make_option
import sys
import numpy

def main(options,args):
    lastConfig=set([1.])

    npoints,min,max=[float(f) for f in options.scanRange.split(",") if f!="" ]
    step=(max-min)/npoints
            
    if options.input != "":
        bstring=""
        if "," in options.input:
            bstring=options.input
        else:
            inp=open(options.input)
            for l in inp.read().split("\n"):
                if options.vectorName in l:
                    bstring=l.lstrip("").split("=")[1]
        lastConfig=set([float(b) for b in bstring.split(",") if b!=""])

    points=numpy.arange(min,max,step)
    ip=0
    for p in points:
        newb = round(p,3)

        thisPoint=lastConfig.copy()
        thisPoint.add(newb)
        if len(thisPoint) < 2:
            continue
        ip+=1

        bkgPolySrt="bkgPolOrderByCat=4,5,5,5"
        
        f = open("%s/point%d.dat" % ( options.output,ip ), "w+" ) 
        f.write( "%s=" % (options.vectorName) )
        srted=sorted(thisPoint)
        f.write("%1.4g" % srted[0])
        for b in srted[1:]:
            f.write(",%1.4g" % b)
            bkgPolySrt+=",2"
        f.write("\n")
        f.write("%s\n" % bkgPolySrt)
        f.close()


if __name__ == "__main__":
    parser = OptionParser(option_list=[
        make_option("-i", "--input",
                    action="store", type="string", dest="input",
                    default="1.,",
                    help="default=%default", metavar=""
                    ),
        make_option("-o", "--output",
                    action="store", type="string", dest="output",
                    default="jetanalysis/optimization",
                    help="default=%default", metavar=""
                    ),
        make_option("-r", "--scanRange",
                    action="store", type="string", dest="scanRange",
                    default="90,0.85,1.",
                    help="default=%default", metavar=""
                    ),
        make_option("-s", "--step",
                    action="store", type="int", dest="step",
                    default=0,
                    help="default=%default", metavar=""
                    ),
        make_option("-n", "--vectorName",
                    action="store", type="string", dest="vectorName",
                    default="mvaVbfCatBoundaries",
                    help="default=%default", metavar=""
                    ),
        
        ])

    (options, args) = parser.parse_args()

    sys.argv.append("-b")
    main( options, args ) 
