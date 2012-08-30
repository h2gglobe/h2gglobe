#!/bin/env python

import sys
import ROOT
import numpy
from math import sqrt
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '.3f')

# -----------------------------------------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------------------------------------
def optmizeCats(func, ws, args):
    
    grS = ROOT.TGraph()
    grS.SetName("zVsNcat")
    grS.SetTitle(";n_{cat};f.o.m [A.U.]")
    summary = {}
    for iter in range(1,9):
    ## for iter in [2]:
        boundaries,z = func(iter,ws,*args)
        ncat = iter ## len(boundaries)-1
        if not ncat in summary or summary[ncat]["z"] < z: 
            summary[ncat] =  { "z" : z, "boundaries" : boundaries, "ncat": ncat }
    for ncat,val in summary.iteritems():
        grS.SetPoint( grS.GetN(), ncat, val["z"] )

    canv1 = ROOT.TCanvas("canv1","canv1")
    canv1.cd()
    grS.SetMarkerStyle(ROOT.kFullCircle)
    grS.Draw("AP")

    canv1.SaveAs("vbfmva_opt_fom.png")
    canv1.SaveAs("vbfmva_opt_fom.C")
    
    ws.rooImport(canv1)
    ws.rooImport(grS)
    
    raw_input("Done")

    sout = open("vbfmva_opt.json","w+")
    sout.write( json.dumps(summary,sort_keys=True, indent=4) )

    ws.writeToFile("vbfmva_opt.root")

# -----------------------------------------------------------------------------------------------------------
def optimizeNcat2D(ncat, ws, isCdf, sigPx, sigPy, sigIntegral, bkgIntegral1, bkgIntegral2, cutoffX, cutoffY,
                   sigNorm, bkg1Norm, bkg2Norm, relNorm=0.5, syst=0.5, cnst=0.):

    nbound = ncat+2 ## book one fake category, since the last boundary is not floated by Minuit for some reason
    quantilesX = numpy.zeros(nbound)
    quantilesY = numpy.zeros(nbound)
    prbs      = numpy.arange(0.,1.,1./ncat)
    prbs[0]   = 0.01 
    sigPx.GetQuantiles(ncat,quantilesX,prbs)
    sigPy.GetQuantiles(ncat,quantilesY,prbs)
    quantilesX[ncat] = 1.
    quantilesY[ncat] = 1.

    fomDen = ROOT.TF3("fomden","x[1]+[1]*x[2] + 2*[0]*[1]*sqrt(x[1]+[1]*x[2])*x[2] + [0]*[0]*[1]*[1]*x[2]*x[2] + [2]*[2]*x[0]*x[0]",0.,1.e+3,0.,1.e+3,0.,1.e+3)
    fomDen.SetParameter(0,syst)
    fomDen.SetParameter(1,relNorm)
    fomDen.SetParameter(2,cnst)

    pws = ROOT.PieceWise2DSignif(nbound,isCdf,sigIntegral,bkgIntegral1,bkgIntegral2,fomDen,cutoffX,cutoffY,sigNorm,bkg1Norm,bkg2Norm)
        
    minimizer = ROOT.TMinuitMinimizer()
    minimizer.SetFunction(pws)
    minimizer.SetFixedVariable(0,"fixedBoundX",1.)
    minimizer.SetFixedVariable(nbound,"fixedBoundY",1.)
    for ivar in range(1,nbound):
        minimizer.SetLimitedVariable(ivar,       "deltaBoundX%d" % ivar, quantilesX[nbound-ivar]-quantilesX[nbound-ivar-1],1.e-3,cutoffX,1.)
        minimizer.SetLimitedVariable(ivar+nbound,"deltaBoundY%d" % ivar, quantilesY[nbound-ivar]-quantilesY[nbound-ivar-1],1.e-3,cutoffY,1.)

    minimizer.Minimize()
    minimizer.PrintResults()
    maxval = sqrt(-minimizer.MinValue())
    xvar = minimizer.X()
    boundariesX = [ xvar[0] ]
    boundariesY = [ xvar[nbound] ]
    efficiencies = []
    for ivar in range(1,nbound-1):
        boundariesX.append( boundariesX[-1] - xvar[ivar] )
        boundariesY.append( boundariesY[-1] - xvar[ivar+nbound] )
        efficiencies.append( sigIntegral.Eval(boundariesX[-2], boundariesY[-2]) - sigIntegral.Eval(boundariesX[-1], boundariesY[-1]) )
        
    
    print "---------------------------------------------"
    print "input  ncat: ", ncat
    print "output ncat: ", len(boundariesX)-1, len(boundariesY)-1
    print "max: %1.3g" % ( maxval )
    print "boundariesX: ",
    for b in boundariesX:
        print "%1.3g" % b,
    print
    print "boundariesY: ",
    for b in boundariesY:
        print "%1.3g" % b,
    print
    print "efficiencies: ",
    for b in efficiencies:
        print "%1.3g" % b,
    print
    print

    
    return (boundariesX,boundariesY,efficiencies),maxval


# -----------------------------------------------------------------------------------------------------------
def optimize2D(fin):
    minX = 0.
    minY = 0.
    cutoffX = 0.015
    cutoffY = 0.01

    ws = ROOT.RooWorkspace("ws","ws")
    ws.rooImport = getattr(ws,'import')
    mva0 = ws.factory("mva0[1.,%f.,1.001]" % minX )
    mva2 = ws.factory("mva2[1.,%f.,1.001]" % minY )
    varSet = ROOT.RooArgSet(mva0,mva2)
    varList = ROOT.RooArgList(mva0,mva2)

    mva0.setBins(1001)
    mva2.setBins(1001)

    sig = fin.Get("sig")
    bkg1 = fin.Get("bkg1")
    bkg2 = fin.Get("bkg2")
    
    sigDs = ROOT.RooDataSet("sigDs","sigDs",sig,varSet)
    bkg1Ds = ROOT.RooDataSet("bkg1Ds","bkg1Ds",bkg1,varSet)
    bkg2Ds = ROOT.RooDataSet("bkg2Ds","bkg2Ds",bkg2,varSet)

    sigHisto = sigDs.createHistogram(mva0,mva2)
    bkg1Histo = bkg1Ds.createHistogram(mva0,mva2)
    bkg2Histo = bkg2Ds.createHistogram(mva0,mva2)
    sigHistoX = sigHisto.ProjectionX()
    sigHistoY = sigHisto.ProjectionY()

    print
    print "-----------------------------------------------------------------------------"
    print

    sigHistoIntegral  = ROOT.integrate2D(sigHisto)
    bkg1HistoIntegral = ROOT.integrate2D(bkg1Histo)
    bkg2HistoIntegral = ROOT.integrate2D(bkg2Histo)

    ws.rooImport(mva0)
    ws.rooImport(mva2)
    ws.rooImport(sig)
    ws.rooImport(bkg1)
    ws.rooImport(bkg2)

    ### sigHistoToTF2  = ROOT.SimpleHistoToTF2("sigHistoToTF2",sigHisto)
    ### sigTF2         = ROOT.TF2("sigTF2",sigHistoToTF2,minX,1,minY,1,0,"SimpleHistoToTF2")
    ### bkg1HistoToTF2 = ROOT.SimpleHistoToTF2("bkg1HistoToTF2",bkg1Histo)
    ### bkg1TF2        = ROOT.TF2("bkg1TF2",bkg1HistoToTF2,minX,1,minY,1,0,"SimpleHistoToTF2")
    ### bkg2HistoToTF2 = ROOT.SimpleHistoToTF2("bkg2HistoToTF2",bkg2Histo)
    ### bkg2TF2        = ROOT.TF2("bkg2TF2",bkg2HistoToTF2,minX,1,minY,1,0,"SimpleHistoToTF2")
    ### isCdf = False

    sigHistoToTF2  = ROOT.SimpleHistoToTF2("sigHistoToTF2",sigHistoIntegral)
    sigTF2         = ROOT.TF2("sigTF2",sigHistoToTF2,minX,1,minY,1,0,"SimpleHistoToTF2")
    bkg1HistoToTF2 = ROOT.SimpleHistoToTF2("bkg1HistoToTF2",bkg1HistoIntegral)
    bkg1TF2        = ROOT.TF2("bkg1TF2",bkg1HistoToTF2,minX,1,minY,1,0,"SimpleHistoToTF2")
    bkg2HistoToTF2 = ROOT.SimpleHistoToTF2("bkg2HistoToTF2",bkg2HistoIntegral)
    bkg2TF2        = ROOT.TF2("bkg2TF2",bkg2HistoToTF2,minX,1,minY,1,0,"SimpleHistoToTF2")
    isCdf = True

    canv2 = ROOT.TCanvas("canv2","canv2")
    canv2.cd()
    ## sigHisto.Draw("colz")
    sigTF2.SetLineColor(ROOT.kBlue)
    bkg1TF2.SetLineColor(ROOT.kRed)
    bkg2TF2.SetLineColor(ROOT.kMagenta)
    sigTF2.Draw("")
    bkg1TF2.Draw("SAME")
    bkg2TF2.Draw("SAME")

    raw_input("Go?")
    optmizeCats( optimizeNcat2D, ws, (isCdf,sigHistoX,sigHistoY,sigTF2,bkg1TF2,bkg2TF2,cutoffX,cutoffY,
                                      sigDs.sumEntries(), bkg1Ds.sumEntries(), bkg2Ds.sumEntries() ) )

    return ws

# -----------------------------------------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------------------------------------
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

    return ws

# -----------------------------------------------------------------------------------------------------------
def main(fname):

    fin = ROOT.TFile.Open(fname)
    ROOT.gSystem.SetIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include")
    ROOT.gROOT.LoadMacro("NaiveBoundaryOptimization.C+")

    ### ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    
    ## ws = optimize1D(fin)
    ws = optimize2D(fin)

    
    return ws
    
if __name__ == "__main__":
    
    
    ws=main(sys.argv[1])
