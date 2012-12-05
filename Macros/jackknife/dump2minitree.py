#!/usr/bin/env python

from optparse import OptionParser, make_option
import sys
import json 
from pprint import pprint 
from math import fabs
import operator
import random
import ROOT

minrun=999999
maxrun=0

lastrun=999999

isMVA=False


def mkBranch(tree, nt, name, type):
    tree.Branch(name,ROOT.AddressOf(nt,name),"%s/%s" % (name,type) )

def dumpTree(file,lst):

    file.cd()
    nt = ROOT.Entry()
    tree = ROOT.TTree("hgg_mini_tree","hgg_mini_tree")
    mkBranch(tree,nt,"run","I")
    mkBranch(tree,nt,"lumis","I")
    mkBranch(tree,nt,"event","I")
    mkBranch(tree,nt,"category","I")
    if isMVA: mkBranch(tree,nt,"diphotonMVA","F")
    mkBranch(tree,nt,"CMS_hgg_mass","F")

    for idx,vars in lst.iteritems():
        if len(vars.keys()) == 0:
            continue
        nt.run, nt.lumis, nt.event = idx
        nt.category = vars["evcat"]
        if isMVA: nt.diphotonMVA = vars["diphoBDT"]
        nt.CMS_hgg_mass = vars["mgg"]
        tree.Fill()
    tree.Write()

def getlist(input):
    lst = {}

    global minrun, maxrun, isMVA
    

    for line in input.split("\n"):
        vars = {}
        for i in line.replace("="," ").replace("\t"," ").split(" "):
            try:
                if i != "":
                    j = i.split(":")
                    if j[0] == "run" or j[0] == "lumi" or j[0] == "event":
                        globals()[j[0]] = int(j[1])
                    else:
                        try:
                            vars[j[0]] = float(j[1])
                        except Exception, e:
                            ### print e
                            pass
                    if "diphoBDT" in vars:
                        isMVA = True
                    else:
                        isMVA = False
                        
            except Exception, e:
                print e
        
        try:
            if run > lastrun:
                continue
            
            if run<minrun:
                minrun=run
            if run>maxrun:
                maxrun=run

            if "diphoBDT" in vars and vars["diphoBDT"] < -0.05:
                  continue
            lst[  (run, lumi, event) ] = vars
        except Exception, e:
            pass
        

    return lst

def main(options,args):

    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetPalette(1)
    ROOT.gStyle.SetOptStat(111111)
    ROOT.gStyle.SetOptFit(1)
    
    ROOT.gROOT.ProcessLine( \
       "struct Entry{          \
           int run;            \
           int lumis;          \
           int event;          \
           int category;       \
           float CMS_hgg_mass; \
           float diphotonMVA;  \
          };"
       )
    
    fn1 = args.pop(0)
    file1 = open(fn1)

    fn2 = ""
    file2 = None
    if len(args)>0:
        fn2 = args.pop(0)
        file2 = open(fn2)
    
    print "reading list1"
    list1 = getlist( file1.read() )
    if options.makeTrees:
        print "isMVA ",isMVA
        print "making mini tree"
        fout1 = ROOT.TFile.Open(fn1.replace("txt","root"),"recreate")
        dumpTree(fout1,list1)
        fout1.Close()


    if file2:
        print "reading list2"
        list2 = getlist( file2.read() )
        if options.makeTrees:
            print "isMVA ",isMVA
            print "making mini tree"
            fout2 = ROOT.TFile.Open(fn2.replace("txt","root"),"recreate")
            dumpTree(fout2,list2)
            fout2.Close()
    else:
        list2 = {}
        
    print "getting event lists"
    events1 = list1.keys()
    events2 = list2.keys()
    
    print "computing union, intersection and differences"
    common = set(events1).intersection(  set(events2) )
    all    = list(set(events1).union(  set(events2) ))
    
    only1 = list(set(events1) -  set(events2))
    only2 = list(set(events2) -  set(events1))
    
    only1.sort()
    only2.sort()
    
    for id in options.d:
        ## all.sort()
        random.shuffle(all)

        d = int(id)
        n = len(all)
        g = n / d
        print "computing partitions: n=%d d=%d g=%d" % (n,d,g)
        if n % d != 0:
            g += 1
        parts = [ [] for k in range(g) ]
        for j in range(n):
            parts[ j % g ].append( all[j] )

        open("%s/partitions_%d.json" % (options.jsondir,d),"w+").write(json.dumps(parts))



if __name__ == "__main__":
    parser = OptionParser(option_list=[
        make_option("-d", "--delete",
                    action="append", dest="d",
                    default=[],
                    help="default: %default", metavar=""
                    ),
        make_option("-t", "--makeTrees",
                    action="store_true", dest="makeTrees",
                    default=False,
                    help="default: %default", metavar=""
                    ),
        make_option("-o", "--jsondir",
                    action="store", type="string", dest="jsondir",
                    default="./",
                    help="default : %default", metavar=""
                    ),
        ])
    
    (options, args) = parser.parse_args()

    print options, args

    sys.argv.append("-b")
    main( options, args ) 
