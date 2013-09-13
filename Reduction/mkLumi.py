#!/bin/env python

from mkPileup2D import Parallel, cpu_count

from optparse import OptionParser, make_option
import sys, os
import commands
import array

from python.makeFilelist import makeCaFiles, black_list
from python.lumi import dumpLumi

def run(fname,parallel):

    remote = os.environ.get("storeremote")
    
    fin = open(fname)
    vals = {}
    for l in fin.read().split("\n"):
        if "output=" in l :
            for t in l.split(" "):
                if "=" in t:
                    nam,val = t.split("=",1)
                    vals[nam] = val
    dir = "%s/%s" % ( remote, os.path.dirname(vals["output"]).replace("./datastore",""))
    output = "%s.json" % os.path.dirname(vals["output"])
    
    print "listing ", dir
    files = [ f[0] for f in makeCaFiles(dir) if f[0] not in black_list ]
    print "filling", output
    if parallel.njobs == 1:
        makeLumi([output]+files)
    else:
        cmd = "./mkLumi.py --do"
        for a in [output]+files:
            cmd += " %s " % a
        sta,out= commands.getstatusoutput(cmd)
        print "filled", output, out
        return sta,out
    return None

def runAll(args,options=None):

    import ROOT
    ROOT.gROOT.ProcessLine(".L ../Macros/pileup/pileup.C+g")
    
    parallel = Parallel(cpu_count(),len(args))
    ## parallel = Parallel(1)
    
    for fname in args:
        parallel.run( run, (fname,parallel) )
    
def makeLumi(args):

    output = args[0]
    outdir = os.path.dirname(output)
    if not os.path.exists(outdir):
        commands.getstatusoutput("mkdir -p %s" % outdir)

    filelist = args[1:]

    runLines = dumpLumi(filelist)

    out = open(output,"w+")
    out.write( "{" + ", ".join(runLines) + "}\n" )
    out.close()
    
if __name__ == "__main__":

    parser = OptionParser(option_list=[
        make_option("--do", 
                    action="store_true", dest="do",
                    default=False,
                    help="Input file name", metavar=""
                    ),
        ]
                          )

    (options, args) = parser.parse_args()
    ## print options, args
    sys.argv.append("-b")

    if options.do:
        makeLumi(args)
    else:
        runAll(args)

    
