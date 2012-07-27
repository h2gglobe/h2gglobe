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



def optimizeNcat(ncat, ws, mva, sig, sigIntegral, bkgIntegral, xmin, cutoff):

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
    for ivar in range(1,nbound):
        minimizer.SetLimitedVariable(ivar,"deltaBound%d" % ivar, quantiles[nbound-ivar]-quantiles[nbound-ivar-1],1.e-3,cutoff,1.)

    minimizer.Minimize()
    minimizer.PrintResults()
    maxval = sqrt(-minimizer.MinValue())
    xvar = minimizer.X()
    boundaries = [ xvar[0] ]
    for ivar in range(1,nbound):
        boundaries.append( boundaries[-1] - xvar[ivar] )
        
    ### for ibound in range(nbound):
    ###     minimizer.SetLimitedVariable(ibound,"bound%d" % ibound,quantiles[ibound],1.e-3,xmin,0.999)
    ### minimizer.SetFixedVariable(nbound-1,"fixedBound",1.)

    #### minimizer.Minimize()
    #### maxval = sqrt(-minimizer.MinValue())
    #### xvar = minimizer.X()
    ### boundaries = sorted(set([ round(xvar[i]*400.)/400. for i in range(nbound) ]))

    #### pwsMin = ROOT.PieceWiseSignif(len(boundaries),sigIntegral,bkgIntegral,cutoff)
    #### minimizer.SetFunction(pwsMin)
    #### for ibound in range(len(boundaries)):
    ####     minimizer.SetLimitedVariable(ibound,"bound%d" % ibound,boundaries[ibound],1.e-3,xmin,0.999)
    #### minimizer.SetFixedVariable(len(boundaries)-1,"fixedBound",1.)
    #### 
    #### minimizer.Minimize()
    #### maxval2 = sqrt(-minimizer.MinValue())
    #### xvar = minimizer.X()
    #### boundaries2 = sorted(set([ round(xvar[i]*300.)/300. for i in range(nbound) ]))
    #### 
    #### if abs(1. - maxval2/maxval) > 5.e-2 and len(boundaries2) < len(boundaries):
    ####     boundaries = boundaries2
    ####     maxval = maxval2

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
        
def main(fname):

    fin = ROOT.TFile.Open(fname)
    ROOT.gSystem.SetIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include")
    ROOT.gROOT.LoadMacro("NaiveBoundaryOptimization.C+")

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

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

    ws = ROOT.RooWorkspace()
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

    ### sigHistoIntegral = ROOT.HistoIntegral("sigHistoIntegral",  hsig )
    ### sigIntegral = ROOT.TF2("sigIntegral", sigHistoIntegral, xmin, 1., xmin, 1., 1, "HistoIntegral" )
    ### sigIntegral.SetParameter(0,cutoff)
    ### 
    ### bkgHistoIntegral = ROOT.HistoIntegral("bkgHistoIntegral",  hbkg )
    ### bkgIntegral = ROOT.TF2("bkgIntegral", bkgHistoIntegral, xmin, 1., xmin, 1., 1, "HistoIntegral" )
    ### bkgIntegral.SetParameter(0,cutoff)

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

    grS = ROOT.TGraph()
    grS.SetName("zVsNcat")
    summary = {}
    for iter in range(1,9):
        boundaries,z = optimizeNcat(iter,ws,mva,hsig,sigTF1,bkgTF1,xmin,cutoff)
        z*=sqrt(6.7)
        ncat = len(boundaries)-1
        if not ncat in summary or summary[ncat]["z"] < z: 
            summary[ncat] =  { "z" : z, "boundaries" : boundaries, "ncat": ncat }
    for ncat,val in summary.iteritems():
        grS.SetPoint( grS.GetN(), ncat, val["z"] )

    canv1 = ROOT.TCanvas("canv1","canv1")
    canv1.cd()
    grS.SetMarkerStyle(ROOT.kFullCircle)
    grS.Draw("AP")

    canv2 = ROOT.TCanvas("canv2","canv2")
    canv2.cd()
    sigTF1.SetLineColor(ROOT.kRed)
    sigTF1.Draw("")
    bkgTF1.Draw("SAME")

    ws.rooImport(canv1)
    ws.rooImport(canv2)
    ws.rooImport(grS)
    
    raw_input("aa")

    sout = open("vbfmva_opt.json","w+")
    sout.write( json.dumps(summary) )

    ws.writeToFile("vbfmva_opt.root")
    
    return ws
    
if __name__ == "__main__":
    
    
    ws=main(sys.argv[1])
