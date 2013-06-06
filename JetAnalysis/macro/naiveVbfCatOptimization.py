#!/bin/env python

import sys
import ROOT
import numpy
from math import sqrt
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '.3f')


# -----------------------------------------------------------------------------------------------------------
def getQuantilesGraphs(histo,probs,twosided=False):
    ## histo.Print("all")
    graphs = [ ROOT.TGraphErrors(histo.GetNbinsX()) for p in probs ]
    if twosided:
        qtiles = []
        for p in probs:
            t = 0.5 - p*0.5
            qtiles.append( t )
            qtiles.append( 1-t )
    else:
        qtiles=probs
        
    nq = len(qtiles)
    graph = ROOT.TGraph(nq+2)
    for iq in range(nq):
        graph.SetPoint(iq,qtiles[iq],0.)
    for g in graphs:        
        g.GetXaxis().SetTitle(histo.GetXaxis().GetTitle())
        g.SetMarkerStyle(histo.GetMarkerStyle())
    graph.SetPoint(nq,0.25,0.)
    graph.SetPoint(nq+1,0.75,0.)
        
    for ix in range(1,histo.GetNbinsX()+1):
        proj = histo.ProjectionY("qtiles",ix,ix)
        
        ## proj.Print("all")
        ## graph.Print("all")
        proj.GetQuantiles(nq+2,graph.GetY(),graph.GetX())
        ntot = proj.Integral()
        if ntot == 0: continue
        h = 1.2*( graph.GetY()[nq+1] - graph.GetY()[nq] ) * pow(ntot,-0.2)
        
        if twosided:
            for ig in range(nq/2):                
                quant1 = graph.GetY()[ig]
                quant2 = graph.GetY()[ig+1]
                quant = (quant2 - quant1)*0.5                
                quant1mh = proj.FindBin( quant1 - h*0.5 )
                quant1ph = proj.FindBin( quant1 + h*0.5 )
                quant2mh = proj.FindBin( quant2 - h*0.5 )
                quant2ph = proj.FindBin( quant2 + h*0.5 )
                
                nint = proj.Integral( quant1mh, quant1ph ) + proj.Integral( quant2mh, quant2ph )
                if nint > 0 and ntot > 0:
                    fq = nint / (2.*h*ntot)
                    err = 1./(2.*sqrt(ntot)*fq)
                else:
                    err = 0.

                graphs[ig/2].SetPoint(ix-1,histo.GetXaxis().GetBinCenter(ix),quant)
                graphs[ig/2].SetPointError(ix-1,histo.GetXaxis().GetBinWidth(ix)*0.5,err)
                
        else:
            for ig in range(nq):
                quant = graph.GetY()[ig]
                quantmh = proj.FindBin( quant - h )
                quantph = proj.FindBin( quant + h )
                nint = proj.Integral( quantmh, quantph )
                
                if nint > 0 and ntot > 0:
                    fq = nint / (2.*h*ntot)
                    err = 1./(2.*sqrt(ntot)*fq)
                else:
                    err = 0.
                
                graphs[ig].SetPoint(ix-1,histo.GetXaxis().GetBinCenter(ix),quant)
                graphs[ig].SetPointError(ix-1,histo.GetXaxis().GetBinWidth(ix)*0.5,err)
                
    return graphs

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
    for iter in range(1,6):
        boundaries,z = func(iter,ws,*args)
        ncat = iter
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
                   sigNorm, bkg1Norm, bkg2Norm, relNorm=0.1, syst=0.5, cnst=0.06):

    nbound = ncat+2 ## book one fake category, since the last boundary is not floated by Minuit for some reason
    quantilesX = numpy.zeros(nbound)
    quantilesY = numpy.zeros(nbound)
    prbs      = numpy.arange(0.,1.,1./ncat)
    prbs[0]   = 0.01 
    sigPx.GetQuantiles(ncat,quantilesX,prbs)
    sigPy.GetQuantiles(ncat,quantilesY,prbs)
    quantilesX[ncat] = 1.
    quantilesY[ncat] = 1.
    
    fomDen = ROOT.TF3("fomden","x[1]+[1]*x[2]+[2]*x[0] + 2*[0]*[1]*sqrt(x[1]+[1]*x[2]+[2]*x[0])*x[2] + [0]*[0]*[1]*[1]*x[2]*x[2]",0.,1.e+3,0.,1.e+3,0.,1.e+3) ## denominator of z^2 (s,b1,b2 = x[0:3])
    fomDen.SetParameter(0,syst)
    fomDen.SetParameter(1,relNorm)
    fomDen.SetParameter(2,cnst)

    pws = ROOT.PieceWise2DSignif(nbound,isCdf,sigIntegral,bkgIntegral1,bkgIntegral2,fomDen,cutoffX,cutoffY,sigNorm,bkg1Norm,bkg2Norm)
        
    minimizer = ROOT.TMinuitMinimizer()
    minimizer.SetFunction(pws)
    minimizer.SetFixedVariable(0,"fixedBoundX",1.)
    minimizer.SetFixedVariable(nbound,"fixedBoundY",1.)
    for ivar in range(1,nbound):
        minimizer.SetLimitedVariable(ivar,       "deltaBoundX%d" % ivar, quantilesX[nbound-ivar]-quantilesX[nbound-ivar-1],1.e-3,cutoffX,2.)
        minimizer.SetLimitedVariable(ivar+nbound,"deltaBoundY%d" % ivar, quantilesY[nbound-ivar]-quantilesY[nbound-ivar-1],1.e-3,cutoffY,2.)

    minimizer.Minimize()
    minimizer.PrintResults()
    maxval = sqrt(-minimizer.MinValue())
    xvar = minimizer.X()
    boundariesX = [ xvar[0] ]
    boundariesY = [ xvar[nbound] ]
    efficienciesSig = []
    efficienciesBkg1 = []
    efficienciesBkg2 = []
    for ivar in range(1,nbound-1):
        boundariesX.append( boundariesX[-1] - xvar[ivar] )
        boundariesY.append( boundariesY[-1] - xvar[ivar+nbound] )
        efficienciesSig.append( sigIntegral.Eval(boundariesX[-1], boundariesY[-1]) - sigIntegral.Eval(boundariesX[-2], boundariesY[-2]) )
        efficienciesBkg1.append( bkgIntegral1.Eval(boundariesX[-1], boundariesY[-1]) - bkgIntegral1.Eval(boundariesX[-2], boundariesY[-2]) )
        efficienciesBkg2.append( bkgIntegral2.Eval(boundariesX[-1], boundariesY[-1]) - bkgIntegral2.Eval(boundariesX[-2], boundariesY[-2]) )
        
    
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
    print "efficienciesSig: ",
    for b in efficienciesSig:
        print "%1.3g" % b,
    print
    print "efficienciesBkg1: ",
    for b in efficienciesBkg1:
        print "%1.3g" % b,
    print
    print "efficienciesBkg2: ",
    for b in efficienciesBkg2:
        print "%1.3g" % b,
    print
    print

    
    return (boundariesX,boundariesY,efficienciesSig,efficienciesBkg1,efficienciesBkg2),maxval


# -----------------------------------------------------------------------------------------------------------
def optimize2D(fin):
    minX = -1
    minY = -1
    cutoffX = 0.015
    cutoffY = 0.015

    ws = ROOT.RooWorkspace("ws","ws")
    ws.rooImport = getattr(ws,'import')
    mva0 = ws.factory("mva0[1.,%f.,1.001]" % minX )
    mva2 = ws.factory("mva1[1.,%f.,1.001]" % minY )
    varSet = ROOT.RooArgSet(mva0,mva2)
    varList = ROOT.RooArgList(mva0,mva2)

    mva0.setBins(1001)
    mva2.setBins(1001)

    sig = fin.Get("sig")
    bkg1 = fin.Get("bkg0")
    bkg2 = fin.Get("bkg1")
    
    sigDs  = ROOT.RooDataSet("sigDs","sigDs",sig,varSet)
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

    canv3 = ROOT.TCanvas("canv3","canv3",1200,800)
    canv3.Divide(3,2)
    canv3.cd(1)
    sigHistoIntegral.Draw("colz")
    canv3.cd(3+1)
    sigHisto.DrawNormalized("colz")
    canv3.cd(2)
    bkg1HistoIntegral.Draw("colz")
    canv3.cd(3+2)
    bkg1Histo.DrawNormalized("colz")
    canv3.cd(3)
    bkg2HistoIntegral.Draw("colz")
    canv3.cd(3+3)
    bkg2Histo.DrawNormalized("colz")

    canv3.SaveAs("vbfmva_opt_input.png")

    raw_input("Go?")
    optmizeCats( optimizeNcat2D, ws, (isCdf,sigHistoX,sigHistoY,sigTF2,bkg1TF2,bkg2TF2,cutoffX,cutoffY,
                                      sigDs.sumEntries(), bkg1Ds.sumEntries(), bkg2Ds.sumEntries() ) )

    return ws


# -----------------------------------------------------------------------------------------------------------
def optimize1DMassWidth(fin):
    minX = -0.5
    cutoffX = 0.05
    ## cut = "vbfMVA>-2."
    cut = ""
    
    ws = ROOT.RooWorkspace("ws","ws")
    ws.rooImport = getattr(ws,'import')
    diphotonMVA = ws.factory("diphotonMVA[1.,%f.,1.001]" % minX )
    vbfMVA = ws.factory("vbfMVA[1.,-2.,1.001]" )
    mass = ws.factory("CMS_hgg_mass[125.,110.,150.]")
    weight = ws.factory("weight[1.,0.,10.]")
    varSet = ROOT.RooArgSet(diphotonMVA,vbfMVA,mass,weight)
    varList = ROOT.RooArgList(diphotonMVA,vbfMVA,mass,weight)

    diphotonMVA.setBins(100)
    mass.setBins(320)

    sig = [ "ggh_m124_8TeV" , "vbf_m124_8TeV", "wzh_m124_8TeV", "tth_m124_8TeV" ]
    bkg = [ "gjet_40_8TeV_pf", "diphojet_8TeV" ]
    
    sig0 = fin.Get(sig.pop(0))
    bkg0 = fin.Get(bkg.pop(0))

    sigHisto = ROOT.TH2D("sigHisto","sigHisto",100,minX,1.,320,110,150)
    bkgHisto = ROOT.TH2D("bkgHisto","bkgHisto",100,minX,1.,320,110,150)
    tmpHisto = ROOT.TH2D("tmpHisto","tmpHisto",100,minX,1.,320,110,150)

    tweight = "weight"
    if cut != "":
        tweight += " * ( %s )" % cut
    sig0.Draw("CMS_hgg_mass:diphotonMVA>>sigHisto",tweight)
    bkg0.Draw("CMS_hgg_mass:diphotonMVA>>bkgHisto",tweight)
    
    sigDs = ROOT.RooDataSet("sigDs","sigDs",sig0,varSet,cut,"weight")
    bkgDs = ROOT.RooDataSet("bkgDs","bkgDs",bkg0,varSet,cut,"weight")

    for sample in sig:
        sampleTree = fin.Get(sample)
        sampleDs = ROOT.RooDataSet(sample,sample,sampleTree,varSet,cut,"weight")
        sigDs.append(sampleDs)
        sampleTree.Draw("CMS_hgg_mass:diphotonMVA>>tmpHisto",tweight)
        sigHisto.Add(tmpHisto)
        
    for sample in bkg:
        sampleTree = fin.Get(sample)
        sampleDs = ROOT.RooDataSet(sample,sample,sampleTree,varSet,cut,"weight")
        bkgDs.append(sampleDs)
        sampleTree.Draw("CMS_hgg_mass:diphotonMVA>>tmpHisto",tweight)
        bkgHisto.Add(tmpHisto)
        
    sigHistoX = sigHisto.ProjectionX()
    bkgHistoX = bkgHisto.ProjectionX()
    bkgHistoX.SetLineColor(ROOT.kRed)
    sigHistoX.Scale(1./sigHistoX.Integral())
    bkgHistoX.Scale(1./bkgHistoX.Integral())

    print sigDs.sumEntries(), sigHisto.Integral()
    print bkgDs.sumEntries(), bkgHisto.Integral()


    c = ROOT.TCanvas()
    c.cd()
    sigHistoX.DrawNormalized("")
    bkgHistoX.DrawNormalized("same")
    raw_input("Go?")
        
    sigCdfX = ROOT.integrate1D(sigHistoX)
    sigW = getQuantilesGraphs(sigHisto, [0.683], True)[0]

    sigHisto.RebinX(sigHisto.GetNbinsX())
    sigW0 = getQuantilesGraphs(sigHisto, [0.683], True)[0]
    sigW0.Print("all")
    sigW0 = sigW0.GetY()[0]
    
    f = ROOT.TCanvas()
    f.cd()
    sigW.DrawClone("AP")

    for ibin in range(sigW.GetN()):
        x,y = sigW.GetX()[ibin], sigW.GetY()[ibin]
        if y == 0.:
            continue
        xerr, yerr = sigW.GetEX()[ibin], sigW.GetEY()[ibin]
        yerr /= y
        y *= y * sigHistoX.GetBinContent(ibin) / ( sigHistoX.Integral() * sigW0 * sigW0 )
        yerr *= 2.*y
        sigW.SetPoint(ibin, x, y)
        sigW.SetPointError(ibin, xerr, yerr)

    d = ROOT.TCanvas()
    d.cd()
    sigW.Draw("AP")
    e = ROOT.TCanvas()
    e.cd()
    sigCdfX.Draw("")
    raw_input("Go?")

    sigHistoToTF1 = ROOT.HistoToTF1("sigHToTF1",sigHistoX)
    sigTF1        = ROOT.TF1("sigTF1",sigHistoToTF1,-1.,1,0,"HistoToTF1")

    bkgHistoToTF1 = ROOT.HistoToTF1("bkgHToTF1",bkgHistoX)
    bkgTF1        = ROOT.TF1("bkgTF1",bkgHistoToTF1,-1.,1,0,"HistoToTF1")

    sigWFactor = ROOT.SignalWidthFactor("sigWFactor",sigW,sigCdfX)
    
    optmizeCats( optimizeNcat, ws, (sigHistoX,sigTF1,bkgTF1,minX,cutoffX,sigWFactor) ) 

    canv2 = ROOT.TCanvas("canv2","canv2")
    canv2.cd()
    sigTF1.SetLineColor(ROOT.kRed)
    sigTF1.Draw("")
    bkgTF1.Draw("SAME")
    ws.rooImport(canv2)

    return ws

# -----------------------------------------------------------------------------------------------------------
def optimizeNcat(ncat, ws, sig, sigIntegral, bkgIntegral, xmin, cutoff, sigW=None):

    nbound = ncat+1
    quantiles = numpy.zeros(nbound)
    prbs      = numpy.arange(0.,1.,1./ncat)
    prbs[0]   = 0.01 
    sig.GetQuantiles(ncat,quantiles,prbs)
    quantiles[ncat] = 1.

    if sigW:
        pws = ROOT.PieceWiseSignif(nbound,sigIntegral,bkgIntegral,cutoff,sigW)
    else:
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
    
    histos = readPlot(fin,"vbf_mva",0,[""],["vbf_m125_8TeV","diphojet_8TeV"])[0] ## ,"ggh_m125_8TeV"
    ## histos = readPlot(fin,"vbf_mva",0,[""],["vbf_m125_8TeV","Data"])[0] ## ,"ggh_m125_8TeV"
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

    ROOT.gStyle.SetOptStat(0)

    ### ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    
    ## ws = optimize1D(fin)
    ## ws = optimize2D(fin)
    ws = optimize1DMassWidth(fin)

    
    return ws

    
if __name__ == "__main__":
    ws=main(sys.argv[1])
