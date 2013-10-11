#!/bin/env python

import sys, types, os
import numpy
from math import sqrt, log
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '.3f')

from plotVbfVar import *

from optparse import OptionParser, make_option
from  pprint import pprint

objs = []

# -----------------------------------------------------------------------------------------------------------
def main(options,args):

    ## setTDRStyle()
    ROOT.gStyle.SetOptStat(0)
        
    fin = ROOT.TFile.Open(options.infile)

    samples = {
        "vbf" :  {"style" : [(setcolors,ROOT.kBlue) ,("SetFillStyle",0)]},
        "ggh" :  {"style" : [(setcolors,ROOT.kRed)  ,("SetFillStyle",0)]},
        "data" : {"style" : [(setcolors,ROOT.kBlack),("SetFillStyle",0)]}
        }

    catdef = open(options.catdef)
    summary = json.loads(catdef.read())

    cats = summary[options.ncat]["boundaries"]
    print cats
    
    vars,bins = options.vars.split("::")
    variables = [ v for v in vars.split(":") if v != "" ]

    ncat = int(options.ncat)
    boxes = []

    ymins = []
    xmins = []
    for icat in range(ncat+1):
        ymins.append( cats[icat] )
        xmins.append( cats[icat+ncat+1] )

    ymins.sort()
    xmins.sort()
    for icat in range(ncat+1):
        ymin = ymins[icat]
        xmin = xmins[icat]
        print xmin,ymin

        ## box = ROOT.TBox(xmin,ymin,1,1)
        box = ROOT.TBox(1-xmin,1-ymin+1,1-1,1-1)
        applyModifs( box, [("SetLineColor",ROOT.kMagenta+4), ("SetFillStyle",0), ("SetLineWidth",3) ] )
        boxes.append( box  )

    objs.append(boxes)
    
    for s in samples.keys():
        samples[s]["tree"] = fin.Get(s)
        
    for name, cont in samples.iteritems():
        tree = cont["tree"]
        tree.Draw("%s>>h_%s(%s)" % ( vars,name,bins ), "weight", "goff" )
        cont["style"].append(("SetTitle",name))
        cont["style"].append((xtitle,variables[1]))
        cont["style"].append((ytitle,variables[0]))
        
        hist = ROOT.gDirectory.Get("h_%s" % name)
        hist.SetDirectory(0)
        cont["hist"] = hist
        
        applyModifs(hist,cont["style"])

        canv = ROOT.TCanvas("cats_%s_%s" % (options.ncat,name),"cats_%s_%s" % (options.ncat,name) )
        canv.SetLogz()
        canv.SetLogy()
        cont["canv"] = canv

        canv.cd()
        hist.DrawNormalized("colz")
        for box in boxes:
            box.Draw("same")

        for fmt in "C","png","pdf":
            canv.SaveAs("%s/%s.%s" % ( options.outdir, canv.GetName(), fmt) )
        
    objs.append(samples)
        
                     
if __name__ == "__main__":

    parser = OptionParser(option_list=[
            make_option("-i", "--infile",
                        action="store", type="string", dest="infile",
                        default="",
                        help="input file",
                        ),
            make_option("-o", "--outdir",
                        action="store", type="string", dest="outdir",
                        default="",
                        help="categories definition file",
                        ),
            make_option("-d", "--cat-def",
                        action="store", type="string", dest="catdef",
                        default="",
                        help="categories definition file",
                        ),
            make_option("-n", "--ncat",
                        action="store", type="string", dest="ncat",
                        default="",
                        help="number of categories",
                        ),
            make_option("-v", "--variables",
                        action="store", type="string", dest="vars",
                        default="1-vbfMva:1-diphoMva::50,0,2,50,0,2",
                        help="variables"
                        ),
            make_option("-s", "--selections",
                        action="store", type="string", dest="sel",
                        default="",
                        help="variables"
                        ),
            
            ]
                          )

    (options, args) = parser.parse_args()
    ## sys.argv.append("-b")
    
    pprint(options.__dict__)

    import ROOT
    
    main(options,args)
        
