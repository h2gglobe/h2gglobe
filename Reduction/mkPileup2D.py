#!/bin/env python

from optparse import OptionParser, make_option
import sys, os
import commands
import array

from python.makeFilelist import makeCaFiles, black_list


from Queue import Queue

from threading import Thread, Semaphore
from multiprocessing import cpu_count

class Wrap:
    def __init__(self, func, args, queue):
        self.queue = queue
        self.func = func
        self.args = args
        
    def __call__(self):
        ret = self.func( *self.args )
        self.queue.put( ret  )

    
class Parallel:
    def __init__(self,ncpu,njobs):
        self.running = Queue(ncpu)
        self.returned = Queue()
        self.njobs = njobs

    def run(self,cmd,args):
        wrap = Wrap( self, (cmd,args), self.returned )
        thread = Thread(None,wrap)
        thread.start()
        
    def __call__(self,cmd,args):
        if type(cmd) == str:
            print cmd
            for a in args:
                cmd += " %s " % a
            args = (cmd,)
            cmd = commands.getstatusoutput
        self.running.put((cmd,args))
        ret = cmd( *args ) 
        self.running.get()
        self.running.task_done()
        return ret
    
def run(fname,parallel):
    ## print "run", fname
    
    fin = open(fname)
    vals = {}
    for l in fin.read().split("\n"):
        if "output=" in l or "CaDir=" in l:
            for t in l.split(" "):
                if "=" in t:
                    nam,val = t.split("=",1)
                    vals[nam] = val
    dir = vals["CaDir"]
    output = "%s.pileup.root" % os.path.dirname(vals["output"])
    site = vals.get("site","cern.ch")

    print "listing ", dir
    files = [ f[0] for f in makeCaFiles(dir,site=site) if f[0] not in black_list ]
    print "filling", output
    if parallel.njobs == 1:
        makeHistos([output]+files)        
    else:
        cmd = "./mkPileup2D.py --do"
        for a in [output]+files:
            cmd += " %s " % a
        sta,out= commands.getstatusoutput(cmd)
        print "filled", output, out
        return sta,out
    return None
    ## parallel.run( "./mkPileup2D.py --do", [output]+files )


def runAll(args,options=None):

    import ROOT
    ROOT.gROOT.ProcessLine(".L ../Macros/pileup/pileup.C+g")
    
    parallel = Parallel(cpu_count(),len(args))
    ## parallel = Parallel(1)
    
    for fname in args:
        parallel.run( run, (fname,parallel) )
        
def makeHistos(args):
    print "filling ", args[0]
    import ROOT
    ROOT.gSystem.Load("../Macros/pileup/pileup_C.so")
    fillHists = ROOT.fillHists

    runBinning = [0, 197495, 203767, 210000]

    output = args[0]
    outdir = os.path.dirname(output)
    if not os.path.exists(outdir):
        commands.getstatusoutput("mkdir -p %s" % outdir)

    tout = ROOT.TFile.Open(output,"recreate")
    tout.cd()
    h1 = ROOT.TH1D("pileup","pileup",100,0,100)
    h2 = ROOT.TH2D("pu_2D", "pu_2D", 100, 0, 100, 3, array.array('d',runBinning))
    
    chain = ROOT.TChain("event")
    for f in args[1:]:
        chain.AddFile(f)
        
    fillHists(chain,h1,h2)
    
    h1.Write()
    h2.Write()
    tout.Close()
    
    
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
        makeHistos(args)
    else:
        runAll(args)

    
