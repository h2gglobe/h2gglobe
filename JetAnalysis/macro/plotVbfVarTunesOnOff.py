#!/bin/env python

import sys
import os
import ROOT

from plotVbfVar import *

if __name__ == "__main__":

    #
    # open input files
    #
    nominal = [ "TuneZ2Star" ]
    tunes = [ "TuneZ2Star", "TuneUEOFF",  "TuneProQ20", "TuneProPT0",  "TuneP0", "TuneD6T" ]
    nameTemplate = "/tmp/musella/histograms_CMS-HGG_mva_%s.root"

    styles = [
        [(setcolors,ROOT.kBlack),("SetFillStyle",3004),("SetLineWidth",2)],
        [(setcolors,ROOT.kRed),("SetFillStyle",0),("SetLineWidth",2)],
        [(setcolors,ROOT.kBlue),("SetFillStyle",0),("SetLineWidth",2)],
        [(setcolors,ROOT.kGreen),("SetFillStyle",0),("SetLineWidth",2)],
        [(setcolors,ROOT.kMagenta),("SetFillStyle",0),("SetLineWidth",2)],
        [(setcolors,ROOT.kOrange),("SetFillStyle",0),("SetLineWidth",2)],
        ]
    
    files = [ (t,ROOT.TFile(nameTemplate % t),styles[i]) for i,t in enumerate(tunes) ]

    #
    # prepare output folder
    #
    outdir = sys.argv[1]
    try:
        os.mkdir(outdir)
    except:
        pass
    os.chdir(outdir)


    plots = { "cut_VBFLeadJPt_sequential": [],
              "cut_VBFSubJPt_sequential": [],
              "cut_VBF_dEta_sequential": [],
              "cut_VBF_Zep_sequential": [],
              "cut_VBF_Mjj_sequential": [],
              "cut_VBF_dPhi_sequential": [],
              "cut_VBF_DiPhoPtOverM_sequential": [],
              "cut_VBF_LeadPhoPtOverM_sequential": [],
              "cut_VBF_SubPhoPtOverM_sequential": [],
              
              "cut_VBFLeadJPt_nminus1": [],
              "cut_VBFSubJPt_nminus1": [],
              "cut_VBF_dEta_nminus1": [],
              "cut_VBF_Zep_nminus1": [],
              "cut_VBF_Mjj_nminus1": [],
              "cut_VBF_dPhi_nminus1": [],
              "cut_VBF_DiPhoPtOverM_nminus1": [],
              "cut_VBF_LeadPhoPtOverM_nminus1": [],
              "cut_VBF_SubPhoPtOverM_nminus1": [],

              "vbf_mva": [("Rebin",10),(xrange,(-1.,1.))],
              "nvtx"   : [(xrange,(0.,50.))]
             }

    legPos = (0.137,0.525,0.474,0.888)
    objs = []

    print files

    for cat in [0]:
        for plot,plotmodifs in plots.iteritems():
            ggh = []
            qqh = []
            for tune, fin, style in files:
                ggh += [ readProc(fin,subproc={"ggh_m125_8TeV":[]},name="%s_%s" % ( plot, tune ), style=style,title=tune,plot=plot,plotmodifs=plotmodifs,category=cat) ]
                qqh += [ readProc(fin,subproc={"vbf_m125_8TeV":[]},name="%s_%s" % ( plot, tune ), style=style,title=tune,plot=plot,plotmodifs=plotmodifs,category=cat) ]
    
            gghStk = makeStack( "ggh", ggh )
            qqhStk = makeStack( "qqh", qqh )

            ggHcanv,ggHleg = makeCanvAndLeg("ggh_%s_cat%d" % ( plot, cat), legPos )
            ggHcanv.cd()
            drawStack( gghStk, "Draw", "hist nostack")
            
            qqHcanv,qqHleg = makeCanvAndLeg("qqh_%s_cat%d" % ( plot, cat), legPos )
            qqHcanv.cd()
            drawStack( qqhStk, "Draw", "hist nostack")
            
            objs += [ gghStk, qqhStk, ggHcanv, ggHleg, qqHcanv, qqHleg ]
            
            for canv in ggHcanv,qqHcanv:
                for fmt in "png", "C", "pdf":
                    canv.SaveAs( "%s.%s" % ( canv.GetName(), fmt ) )
                    
