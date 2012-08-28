#!/bin/env python

import sys
import ROOT
import numpy
from math import sqrt
import json

def readPlot(fin, name, cat=0, which=[""], samples=["diphojet_8TeV","ggh_m125_8TeV","vbf_m125_8TeV",]):

    ret = []
    for p in which:
        hists = []
        for s in samples:
            nam = "%s%s_cat%d_%s" % ( name, p, cat, s )
            h = fin.Get(nam)
            ## print nam, h
            hists.append(h)
        ret.append(hists)
            
    return ret

def optmizeCats(func, ws, args):
    
    grS = ROOT.TGraph()
    grS.SetName("zVsNcat")
    summary = {}
    for iter in range(1,7):
        boundaries,z = func(iter,ws,*args)
        ncat = len(boundaries)-1
        if not ncat in summary or summary[ncat]["z"] < z: 
            summary[ncat] =  { "z" : z, "boundaries" : boundaries, "ncat": ncat }
    for ncat,val in summary.iteritems():
        grS.SetPoint( grS.GetN(), ncat, val["z"] )

    canv1 = ROOT.TCanvas("canv1","canv1")
    canv1.cd()
    grS.SetMarkerStyle(ROOT.kFullCircle)
    grS.Draw("AP")

    ws.rooImport(canv1)
    ws.rooImport(grS)
    
    raw_input("aa")

    sout = open("vbfmva_opt.json","w+")
    sout.write( json.dumps(summary) )

    ws.writeToFile("vbfmva_opt.root")

def optimizeNcat2D(ncat, ws, sigPx, sigPy, sigIntegral, bkgIntegral1, bkgIntegral2, xmin1, xmin2, cutoff1, cutoff2,
                   sigNorm, bkg1Norm, bkg2Norm, relNorm=0.1, syst=0.5):

    nbound = ncat+1
    quantilesX = numpy.zeros(nbound)
    quantilesY = numpy.zeros(nbound)
    prbs      = numpy.arange(0.,1.,1./ncat)
    prbs[0]   = 0.01 
    sigPx.GetQuantiles(ncat,quantilesX,prbs)
    sigPy.GetQuantiles(ncat,quantilesY,prbs)
    quantilesX[ncat] = 1.
    quantilesY[ncat] = 1.
    print quantilesX, quantilesY

    fomDen = ROOT.TF3("fomden","x[1]+[1]*x[2] + 2*[0]*[1]*sqrt(x[1]+[1]*x[2])*x[2] + [0]*[0]*[1]*[1]*x[2]*x[2]",0.,1.e+3,0.,1.e+3,0.,1.e+3)
    fomDen.SetParameter(0,syst)
    fomDen.SetParameter(1,relNorm)

    pws = ROOT.PieceWise2DSignif(nbound,sigIntegral,bkgIntegral1,bkgIntegral2,fomDen,cutoff1,cutoff2,sigNorm,bkg1Norm,bkg2Norm)
        
    minimizer = ROOT.TMinuitMinimizer()
    minimizer.SetFunction(pws)
    minimizer.SetFixedVariable(0,"fixedBoundX",1.)
    minimizer.SetFixedVariable(nbound,"fixedBoundY",1.)
    for ivar in range(1,nbound):
        minimizer.SetLimitedVariable(ivar,       "deltaBoundX%d" % ivar, quantilesX[nbound-ivar]-quantilesX[nbound-ivar-1],1.e-3,cutoff1,1.)
        minimizer.SetLimitedVariable(ivar+nbound,"deltaBoundY%d" % ivar, quantilesY[nbound-ivar]-quantilesY[nbound-ivar-1],1.e-3,cutoff2,1.)

    minimizer.PrintResults()
    minimizer.Minimize()
    minimizer.PrintResults()
    maxval = sqrt(-minimizer.MinValue())
    xvar = minimizer.X()
    boundariesX = [ xvar[0] ]
    boundariesY = [ xvar[nbound] ]
    for ivar in range(1,nbound):
        boundariesX.append( boundariesX[-1] - xvar[ivar] )
        boundariesY.append( boundariesY[-1] - xvar[ivar+nbound] )
        
    print "---------------------------------------------"
    print "input  ncat: ", ncat
    print "output ncat: ", len(boundariesX)-1, len(boundariesY)-1
    ## print "max: %1.3g %1.3g" % ( maxval, maxval2 )
    print "max: %1.3g" % ( maxval )
    print "boundariesX: ",
    for b in boundariesX:
        print "%1.3g" % b,
    print
    print "boundariesY: ",
    for b in boundariesY:
        print "%1.3g" % b,
    print
    print 
    
    return (boundariesX,boundariesY),maxval



def optimize2D(fin):
    xmin1 = 0.
    xmin2 = 0.
    cutoff1 = 0.015
    cutoff2 = 0.015

    try:
        ws = fin.Get("ws")
        
        sigDs = ws.data("sigDs")
        bkg1Ds = ws.data("bkg1Ds")
        bkg2Ds = ws.data("bkg2Ds")

        sigPdf  = ws.pdf("sigPdf")
        bkg1Pdf = ws.pdf("bkg1Pdf")
        bkg2Pdf = ws.pdf("bkg2Pdf")

        sigHisto = sigDs.createHistogram(mva1,mva2)
        sigHistoX = sigHisto.ProjectionX()
        sigHistoY = sigHisto.ProjectionY()
        
    except:
        ws = ROOT.RooWorkspace("ws","ws")
        ws.rooImport = getattr(ws,'import')
        mva1 = ws.factory("mva1[1.,%f.,1.]" % xmin1 )
        mva2 = ws.factory("mva2[1.,%f.,1.]" % xmin2 )
        varSet = ROOT.RooArgSet(mva1,mva2)
        varList = ROOT.RooArgList(mva1,mva2)

        sig = fin.Get("sig")
        bkg1 = fin.Get("bkg1")
        bkg2 = fin.Get("bkg2")
        
        sigDs = ROOT.RooDataSet("sigDs","sigDs",sig,varSet)
        bkg1Ds = ROOT.RooDataSet("bkg1Ds","bkg1Ds",bkg1,varSet)
        bkg2Ds = ROOT.RooDataSet("bkg2Ds","bkg2Ds",bkg2,varSet)
        
        sigHisto = sigDs.createHistogram(mva1,mva2)
        sigHistoX = sigHisto.ProjectionX()
        sigHistoY = sigHisto.ProjectionY()

        print "Fitting RooKeyPdfs: this takes a while (and they cannot be persisted.....)"
        sigPdf = ROOT.RooNDKeysPdf("sigPdf","sigPdf",varList,sigDs,"amvv")
        print "sigPdf done"
        bkg1Pdf = ROOT.RooNDKeysPdf("bkg1Pdf","bkg1Pdf",varList,bkg1Ds,"amvv")
        print "bkg1Pdf done"
        bkg2Pdf = ROOT.RooNDKeysPdf("bkg2Pdf","bkg2Pdf",varList,bkg2Ds,"amvv")
        print "bkg2Pdf done"
        
        ws.rooImport(mva1)
        ws.rooImport(mva2)
        ws.rooImport(sig)
        ws.rooImport(bkg1)
        ws.rooImport(bkg2)
        ws.rooImport(sigPdf)
        ws.rooImport(bkg1Pdf)
        ws.rooImport(bkg2Pdf)

        ##ws.writeToFile("vbfmva_input.root")

    sigTF2 = sigPdf.asTF(varList)
    bkg1TF2 = bkg1Pdf.asTF(varList)
    bkg2TF2 = bkg2Pdf.asTF(varList)

    canv2 = ROOT.TCanvas("canv2","canv2")
    canv2.cd()
    sigTF2.SetLineColor(ROOT.kRed)
    bkg1TF2.SetLineColor(ROOT.kMagenta)
    sigTF2.Draw("")
    bkg1TF2.Draw("SAME")
    bkg2TF2.Draw("SAME")


    raw_input("aa")
    optmizeCats( optimizeNcat2D, ws, (sigHistoX,sigHistoY,sigTF2,bkg1TF2,bkg2TF2,xmin1,xmin2,cutoff1,cutoff2,
                                      sigDs.sumEntries(), bkg1Pdf.sumEntries(), bkg2Pdf.sumEntries() ) )

    return ws

def optimizeNcat(ncat, ws, sig, sigIntegral, bkgIntegral, xmin, cutoff):

    nbound = ncat+1
    quantiles = numpy.zeros(nbound)
    prbs      = numpy.arange(0.,1.,1./ncat)
    prbs[0]   = 0.01 
    sig.GetQuantiles(ncat,quantiles,prbs)
    quantiles[ncat] = 1.

    pws = ROOT.PieceWiseSignif(nbound,sigIntegral,bkgIntegral,cutoff)
    
    minimizer = ROOT.TMinuitMinimizer()
    minimizer.SetFunction(pws)
    minimizer.SetFixedVariable(0,"fixedBound",1.)
    minimizer.PrintResults()
    for ivar in range(1,nbound):
        minimizer.SetLimitedVariable(ivar,"deltaBound%d" % ivar, quantiles[nbound-ivar]-quantiles[nbound-ivar-1],1.e-3,cutoff,1.)

    minimizer.Minimize()
    minimizer.PrintResults()
    maxval = sqrt(-minimizer.MinValue())
    xvar = minimizer.X()
    boundaries = [ xvar[0] ]
    for ivar in range(1,nbound):
        boundaries.append( boundaries[-1] - xvar[ivar] )
        
    print "---------------------------------------------"
    print "input  ncat: ", ncat
    print "output ncat: ", len(boundaries)-1
    ## print "max: %1.3g %1.3g" % ( maxval, maxval2 )
    print "max: %1.3g" % ( maxval )
    print "boundaries: ",
    for b in boundaries:
        print "%1.3g" % b,
    print
    print 
    
    return boundaries,maxval

def optimize1D(fin):
    xmin = -1.
    cutoff = 0.015
    
    ## histos = readPlot(fin,"vbf_mva",0,[""],["vbf_m125_8TeV","diphojet_8TeV"])[0] ## ,"ggh_m125_8TeV"
    histos = readPlot(fin,"vbf_mva",0,[""],["vbf_m125_8TeV","Data"])[0] ## ,"ggh_m125_8TeV"
    ## histos = readPlot(fin,"vbf_mva",0,[""],["vbf_m125_8TeV","ggh_m125_8TeV"])[0]
    sig = histos[0]
    bkg = histos[1].Clone()
    ## bkg.Scale(2./80.)
    bkg.Scale(1./80.)
    for b in histos[2:]:
        bkg.Add(b)

    ws = ROOT.RooWorkspace("ws","ws")
    ws.rooImport = getattr(ws,'import')
    mva = ws.factory("mva[0.,%f.,1.]" % -1. )
    xvar = ROOT.RooArgList(mva)
    
    sigDset = ROOT.RooDataHist("sigDset","sigDset",xvar,sig)
    bkgDset = ROOT.RooDataHist("bkgDset","bkgDset",xvar,bkg)

    ws.rooImport(sigDset)
    ws.rooImport(bkgDset)
    
    hsig = sigDset.createHistogram("hsig",mva)
    hbkg = bkgDset.createHistogram("hbkg",mva)
    hbkg.Smooth(30)

    sigHistoToTF1 = ROOT.HistoToTF1("sigHToTF1",hsig)
    sigTF1        = ROOT.TF1("sigTF1",sigHistoToTF1,-1.,1,0,"HistoToTF1")

    bkgHistoToTF1 = ROOT.HistoToTF1("bkgHToTF1",hbkg)
    bkgTF1        = ROOT.TF1("bkgTF1",bkgHistoToTF1,-1.,1,0,"HistoToTF1")

    print sigTF1.Eval(1.)
    print bkgTF1.Eval(1.)

    ## ws.rooImport(canv)
    ws.rooImport(sigTF1)
    ws.rooImport(bkgTF1)

    ## raw_input("aa")

    optmizeCats( optimizeNcat, ws, (hsig,sigTF1,bkgTF1,xmin,cutoff) ) 

    canv2 = ROOT.TCanvas("canv2","canv2")
    canv2.cd()
    sigTF1.SetLineColor(ROOT.kRed)
    sigTF1.Draw("")
    bkgTF1.Draw("SAME")
    ws.rooImport(canv2)
    
    #### grS = ROOT.TGraph()
    #### grS.SetName("zVsNcat")
    #### summary = {}
    #### for iter in range(1,9):
    ####     boundaries,z = optimizeNcat(iter,ws,mva,hsig,sigTF1,bkgTF1,xmin,cutoff)
    ####     z*=sqrt(6.7)
    ####     ncat = len(boundaries)-1
    ####     if not ncat in summary or summary[ncat]["z"] < z: 
    ####         summary[ncat] =  { "z" : z, "boundaries" : boundaries, "ncat": ncat }
    #### for ncat,val in summary.iteritems():
    ####     grS.SetPoint( grS.GetN(), ncat, val["z"] )
    #### 
    #### canv1 = ROOT.TCanvas("canv1","canv1")
    #### canv1.cd()
    #### grS.SetMarkerStyle(ROOT.kFullCircle)
    #### grS.Draw("AP")
    #### 
    #### canv2 = ROOT.TCanvas("canv2","canv2")
    #### canv2.cd()
    #### sigTF1.SetLineColor(ROOT.kRed)
    #### sigTF1.Draw("")
    #### bkgTF1.Draw("SAME")
    #### 
    #### ws.rooImport(canv1)
    #### ws.rooImport(canv2)
    #### ws.rooImport(grS)
    #### 
    #### raw_input("aa")
    #### 
    #### sout = open("vbfmva_opt.json","w+")
    #### sout.write( json.dumps(summary) )
    #### 
    #### ws.writeToFile("vbfmva_opt.root")

    return ws
        
def main(fname):

    fin = ROOT.TFile.Open(fname)
    ROOT.gSystem.SetIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include")
    ROOT.gROOT.LoadMacro("NaiveBoundaryOptimization.C+")

    ## ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    
    ## ws = optimize1D(fin)
    ws = optimize2D(fin)

    
    return ws
    
if __name__ == "__main__":
    
    
    ws=main(sys.argv[1])
