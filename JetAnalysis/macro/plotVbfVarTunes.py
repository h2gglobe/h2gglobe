#!/bin/env python

import sys
import os
import ROOT
from array import array

from plotVbfVar import *

def normalize(h):
    ## h.Scale( 1./h.Integral() )
    pass

def overflow(h,den=False):
    return
    nbins = h.GetNbinsX()
    if den:
        h.SetBinContent( nbins, (h.GetBinContent(nbins+1)*h.GetBinWidth(nbins+1)*nbins +h.GetBinContent(nbins)*h.GetBinWidth(nbins)) / (h.GetBinWidth(nbins+1)+h.GetBinWidth(nbins+1)*nbins) )
        
    else:
        h.SetBinContent( nbins, h.GetBinContent(nbins+1)+h.GetBinContent(nbins) )
    h.SetBinContent( nbins+1, 0. )

def rebin(h,bins):
    nb =  int(len(bins.tolist())-1)
    print nb, bins
    rb = h.Rebin( nb,"pippo",bins )
    rb.Print()
    return rb

def density(h,presel=None):
    for ibin in range(h.GetNbinsX()+1):
        lowEdge = h.GetBinLowEdge(ibin)
        width = h.GetBinWidth(ibin)
        if presel and presel>lowEdge and presel<lowEdge+width:
            width -= (presel-lowEdge)
        h.SetBinContent( ibin, h.GetBinContent(ibin) / width )

if __name__ == "__main__":

    setTDRStyle()
    
    #
    # open input files
    #
    tunes = [ "TuneZ2Star", "TuneProQ20", "TuneProPT0",  "TuneP0", "TuneD6T" ]
    ## nameTemplate = "/tmp/musella/histograms_CMS-HGG_mva_%s.root"
    nameTemplate = "/afs/cern.ch/user/m/musella/Analysis/CMGTools/CMSSW_5_2_3_patch2/src/h2gglobe_HEAD/Reduction/AnalysisScripts/yr3_systematics_v1/cutbtag_%s/histograms_CMS-HGG.root"
    ## nameTemplate = "/afs/cern.ch/user/m/musella/Analysis/CMGTools/CMSSW_5_2_3_patch2/src/h2gglobe_HEAD/Reduction/AnalysisScripts/yr3_systematics_v1/mavtag_%s/histograms_CMS-HGG.root"
    doUEOFF = True
    ### doUEOFF = False
    ## doRatio = True
    doRatio = False
    
    styles = [
        [(setcolors,ROOT.kBlack),("SetFillStyle",0),("SetLineWidth",2)],
        [(setcolors,ROOT.kRed),("SetFillStyle",0),("SetLineWidth",2)],
        [(setcolors,ROOT.kBlue),("SetFillStyle",0),("SetLineWidth",2)],
        [(setcolors,ROOT.kGreen+2),("SetFillStyle",0),("SetLineWidth",2)],
        [(setcolors,ROOT.kOrange),("SetFillStyle",0),("SetLineWidth",2)],
        [(setcolors,ROOT.kMagenta),("SetFillStyle",0),("SetLineWidth",2)],
        ]

    ## styles = [
    ##     [(setcolors,ROOT.kRed+3),("SetFillStyle",0),("SetLineWidth",2)],
    ##     [(setcolors,ROOT.kGreen+3),("SetFillStyle",0),("SetLineWidth",2)],
    ##     [(setcolors,ROOT.kCyan+3),("SetFillStyle",0),("SetLineWidth",2)],
    ##     [(setcolors,ROOT.kOrange-1),("SetFillStyle",0),("SetLineWidth",2)],
    ##     [(setcolors,ROOT.kMagenta),("SetFillStyle",0),("SetLineWidth",2)],
    ##     [(setcolors,ROOT.kOrange),("SetFillStyle",0),("SetLineWidth",2)],
    ##     ]

    files = [ (t,ROOT.TFile(nameTemplate % t),styles[i]) for i,t in enumerate(tunes) ]
    if doUEOFF:
        filesUEOFF = [ (t,ROOT.TFile(nameTemplate % (t+"UEOFF")),styles[i]) for i,t in enumerate(tunes) ]

    #
    # prepare output folder
    #
    outdir = sys.argv[1]
    try:
        os.mkdir(outdir)
    except:
        pass
    os.chdir(outdir)


    plots = { "cut_VBFLeadJPt_sequential": [normalize,("Rebin",4),(density,20.),(overflow,True),
                                            (xtitle,"p_{Tj1} (GeV)"),(ytitle,"dN / d%(xtitle)s")],
              "cut_VBFSubJPt_sequential":  [normalize,("Rebin",4),(density,20.),(overflow,True),
                                            (xtitle,"p_{Tj2} (GeV)"),(ytitle,"dN / d%(xtitle)s")],
              "cut_VBF_dEta_sequential":   [normalize,
                                            (rebin,array('d',[0.,2.]+[3.+i for i in range(5)])),
                                            density,(overflow,True),
                                            (xtitle,"#Delta#eta_{jj}"),(ytitle,"dN / d%(xtitle)s ")],
              "cut_VBF_Zep_sequential":    [normalize,("Rebin",2),density,(overflow,True),
                                            (xtitle,"#eta_{#gamma#gamma} - (#eta_{j1}+#eta_{j2})/2"),(ytitle,"dN / d%(xtitle)s ")],
              "cut_VBF_Mjj_sequential":    [normalize,density,(overflow,True),("Rebin",4),
                                            (xtitle,"M_{jj} (GeV)"),(ytitle,"dN / d%(xtitle)s")],
              "cut_VBF_dPhi_sequential":   [normalize,
                                            ## ("Rebin",2),
                                            ## (rebin,array('d',[0.,2.6]+[2.8+0.2*i for i in range(3)])),
                                            density,(overflow,True),
                                            (xtitle,"#Delta\phi_{jj#gamma#gamma}"),(ytitle,"dN / d%(xtitle)s")],
              ## "cut_VBF_DiPhoPtOverM_sequential": [normalize,density,(overflow,True),(ytitle,"dN / d%(xtitle)s")],
              ## "cut_VBF_LeadPhoPtOverM_sequential": [normalize,density,(overflow,True),(ytitle,"dN / d%(xtitle)s")],
              ## "cut_VBF_SubPhoPtOverM_sequential": [normalize,density,(overflow,True),(ytitle,"dN / d%(xtitle)s")],
              
              "cut_VBFLeadJPt_nminus1": [normalize,("Rebin",4),(density,20.),(overflow,True),
                                         (xtitle,"p_{Tj1} (GeV)"),(ytitle,"dN / d%(xtitle)s")],                 
              "cut_VBFSubJPt_nminus1":  [normalize,("Rebin",4),(density,20.),(overflow,True),
                                         (xtitle,"p_{Tj2} (GeV)"),(ytitle,"dN / d%(xtitle)s")],                  
              "cut_VBF_dEta_nminus1":   [normalize,
                                         (rebin,array('d',[0.,2.]+[3.+i for i in range(5)])),
                                         density,(overflow,True),
                                         (xtitle,"#Delta#eta_{jj}"),(ytitle,"dN / d%(xtitle)s ")],                                 
              "cut_VBF_Zep_nminus1":    [normalize,density,(overflow,True),("Rebin",2),
                                         (xtitle,"#eta_{#gamma#gamma} - (#eta_{j1}+#eta_{j2})/2"),(ytitle,"dN / d%(xtitle)s")], 
              "cut_VBF_Mjj_nminus1":    [normalize,("Rebin",4),density,(overflow,True),
                                         (xtitle,"M_{jj} (GeV)"),(ytitle,"dN / d%(xtitle)s")],                       
              "cut_VBF_dPhi_nminus1":   [normalize,
                                         ## ("Rebin",2),
                                         ## (rebin,array('d',[0.,2.6]+[2.8+0.2*i for i in range(3)])),
                                         density,(overflow,True),
                                         (xtitle,"#Delta\phi_{jj#gamma#gamma}"),(ytitle,"dN / d%(xtitle)s")],
              
              ## "cut_VBF_DiPhoPtOverM_nminus1": [normalize],
              ## "cut_VBF_LeadPhoPtOverM_nminus1": [normalize],
              ## "cut_VBF_SubPhoPtOverM_nminus1": [normalize],

              "vbf_mva": [("Rebin",10),(xrange,(-1.,1.))],
              "nvtx"   : [(xrange,(0.,50.))]
             }

    legPos = (0.137,0.525,0.474,0.888)
    objs = []

    print files
    if doUEOFF:
        print filesUEOFF

    for cat in [0]:
        for plot,plotmodifs in plots.iteritems():
            ggh = []
            qqh = []
            gghUEOFF = []
            qqhUEOFF = []
            for tune, fin, style in files:
                ggh += [ readProc(fin,subproc={"ggh_m125_8TeV":[]},name="%s_%s" % ( plot, tune ), style=style,title=tune,plot=plot,plotmodifs=plotmodifs,category=cat) ]
                qqh += [ readProc(fin,subproc={"vbf_m125_8TeV":[]},name="%s_%s" % ( plot, tune ), style=style,title=tune,plot=plot,plotmodifs=plotmodifs,category=cat) ]
            if doUEOFF:
                for tune, fin, style in filesUEOFF:
                    gghUEOFF += [ readProc(fin,subproc={"ggh_m125_8TeV":[]},name="%s_%sUEOFF" % ( plot, tune ), style=style
                                           ## +[("SetLineStyle",ROOT.kDashed)]
                                           ## +[("SetFillStyle",3001)]
                                           ,title=tune,plot=plot,plotmodifs=plotmodifs,category=cat) ]
                    qqhUEOFF += [ readProc(fin,subproc={"vbf_m125_8TeV":[]},name="%s_%sUEOFF" % ( plot, tune ), style=style
                                           ## +[("SetLineStyle",ROOT.kDashed)]
                                           ## +[("SetFillStyle",3001)]
                                           ,title=tune,plot=plot,plotmodifs=plotmodifs,category=cat) ]

            ### gghStk = makeStack( "ggh", ggh )
            ### qqhStk = makeStack( "qqh", qqh )

            
            if doUEOFF:
                ## gghStk = makeStack( "ggh", ggh+gghUEOFF )
                ## qqhStk = makeStack( "qqh", qqh+qqhUEOFF )
                if doRatio:
                    qqhLegend = makeLegend((0.34,0.73,0.64,0.94))
                    gghLegend = makeLegend((0.34,0.73,0.64,0.94))
                    
                    gghShifted = []
                    qqhShifted = []
                    for itune in range(len(ggh)):
                        gghShifted.append( ratioHisto(gghUEOFF[itune],ggh[itune],"MPI off / MPI on") )
                        qqhShifted.append( ratioHisto(qqhUEOFF[itune],qqh[itune],"MPI off / MPI on") )

                    gghStk = makeStack( "ggh", gghShifted )
                    qqhStk = makeStack( "qqh", qqhShifted )
                    
                    for h in qqhShifted:
                        qqhLegend.AddEntry( h, "", "l" )
                    for h in gghShifted:
                        gghLegend.AddEntry( h, "", "l" )

                else:
                    qqhLegend = makeLegend((0.34,0.87,0.64,0.94))
                    gghLegend = makeLegend((0.34,0.87,0.64,0.94))

                    gghShifted = [ggh[0]]
                    qqhShifted = [qqh[0]]
                    for itune in range(len(ggh)):
                        gghShifted.append( shiftHisto(ggh[0],ggh[itune],gghUEOFF[itune]) )
                        qqhShifted.append( shiftHisto(qqh[0],qqh[itune],qqhUEOFF[itune]) )
                
                    gghStk = makeEnvelope( "ggh", gghShifted, styles[1], styles[1] )
                    qqhStk = makeEnvelope( "qqh", qqhShifted, styles[2], styles[2] )
                
                    qqhLegend.AddEntry( qqhStk.GetHists()[2], "POWHEG+PYTHIA", "l" )
                    qqhLegend.AddEntry( qqhStk.GetHists()[0], "MPI uncertainty", "l" )
                    gghLegend.AddEntry( gghStk.GetHists()[2], "POWHEG+PYTHIA", "l" )
                    gghLegend.AddEntry( gghStk.GetHists()[0], "MPI uncertainty", "l" )
                    
            else:
                gghStk = makeEnvelope( "ggh", ggh, styles[1], styles[1] )
                qqhStk = makeEnvelope( "qqh", qqh, styles[1], styles[1] )

            
            ggHcanv,ggHleg = makeCanvAndLeg("ggh_%s_cat%d" % ( plot, cat), legPos )
            ggHcanv.cd()
            gghFrame = gghStk.GetHists()[0].Clone()
            gghFrame.Reset("ICE")
            gghFrame.SetEntries(0)
            gghFrame.Draw()
            if doRatio:
                drawStack( gghStk, "Draw", "same hist e1 nostack")
            else:
                drawStack( gghStk, "Draw", "same hist nostack")
            gghLegend.Draw("same")
            ## stackTitles( gghStk )
            gghFrame.GetYaxis().SetRangeUser(0.8*gghStk.GetMinimum("nostack"),1.4*gghStk.GetMaximum("nostack"))
            ggHcanv.RedrawAxis()
            
            qqHcanv,qqHleg = makeCanvAndLeg("qqh_%s_cat%d" % ( plot, cat), legPos )
            qqHcanv.cd()
            qqhFrame = qqhStk.GetHists()[0].Clone()
            qqhFrame.Reset("ICE")
            qqhFrame.SetEntries(0)
            qqhFrame.Draw()
            if doRatio:
                drawStack( qqhStk, "Draw", "same hist e1 nostack")
            else:
                drawStack( qqhStk, "Draw", "same hist nostack")
            qqhLegend.Draw("same")
            stackTitles( qqhStk )
            qqhFrame.GetYaxis().SetRangeUser(0.8*gghStk.GetMinimum("nostack"),1.4*qqhStk.GetMaximum("nostack"))
            qqHcanv.RedrawAxis()
            
            objs += [ gghStk, qqhStk, ggHcanv, ggHleg, qqHcanv, qqHleg ]
            
            for canv in ggHcanv,qqHcanv:
                for fmt in "png", "root", "pdf":
                    canv.SaveAs( "%s.%s" % ( str(canv.GetName()).replace("_cat0","").replace("_cut",""), fmt ) )
                    
