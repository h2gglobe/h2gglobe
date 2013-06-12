#!/bin/env python

import sys
import ROOT
import numpy
from math import sqrt, log
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '.3f')


# -----------------------------------------------------------------------------------------------------------
def getQuantilesGraphs(histo,probs,twosided=False):
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
def getSecondOrder(h):
    proj = h.ProjectionX()
    prof = h.ProfileX('s')
    effs = getQuantilesGraphs( h, [0.683], True)[0]
    hsum  = proj.Clone()
    hsum2 = proj.Clone()
    for ibin in range(0,hsum.GetNbinsX()):
        norm = proj.GetBinContent(ibin)
        mean = prof.GetBinContent(ibin)
        rms = effs.GetY()[ibin]
        hsum.SetBinContent(ibin, norm*mean)
        hsum2.SetBinContent(ibin, norm*(mean*mean-rms*rms))

    return proj, hsum, hsum2

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
def optmizeCats(optimizer,ws,rng,args,readBack=False,reduce=False):
    
    summary = {}
    if readBack:
        try:
            sin = open("cat_opt.json","r")
            summary = json.loads(sin.read())
            sin.close()
        except:
            summary = {}
        
    for iter in rng:
        optimizer.optimizeNCat(iter,*args)
        boundaries = numpy.array([0. for i in range(iter+1) ])
        z = optimizer.getBoundaries(iter, boundaries)
        if not iter in summary or summary[iter]["fom"] < z: 
            summary[iter] =  { "fom" : z, "boundaries" : list(boundaries), "ncat": iter }

    for ncat,val in summary.iteritems():
        printBoundaries(val["boundaries"],val["fom"])

    if reduce:
        maxncat = 0
        for ncat,val in summary.iteritems():
            boundaries = numpy.array([float(b) for b in val["boundaries"]])
            optimizer.reduce(int(ncat)+1, boundaries, numpy.array([0. for i in range(len(boundaries))]) )
            maxncat = max(int(ncat),maxncat)
        rng = range(1,maxncat+1)
        summary = {}
    
    for iter in rng:
        boundaries = numpy.array([0. for i in range(iter+1) ])
        z = optimizer.getBoundaries(iter, boundaries)
        printBoundaries(boundaries,z)
        ncat = iter
        if not ncat in summary or summary[ncat]["fom"] < z: 
            summary[ncat] =  { "fom" : z, "boundaries" : list(boundaries), "ncat": ncat }
                   
    sout = open("cat_opt.json","w+")
    sout.write( json.dumps(summary,sort_keys=True, indent=4) )

    return summary

# -----------------------------------------------------------------------------------------------------------
def printBoundaries(boundaries, maxval):
    
    print "---------------------------------------------"
    print "ncat: ", len(boundaries)-1
    print "max: %1.5f" % ( maxval )
    print "boundaries: ",
    for b in boundaries:
        print "%1.3g" % b,
    print
    print 


# -----------------------------------------------------------------------------------------------------------
def optimize1D(fin):
    minX = 0.
    maxX = 1.
    nbinsX = int((maxX-minX)/0.01)
    cutoff = 0.015
    ndim = 1

    minNormX = 0.
    maxNormX = 1.
        
    minMass = 100.
    maxMass = 180.
    nMassBins = 160
    ## nMassBins = 320
        
    ### histos = readPlot(fin,"vbf_mva",0,[""],["vbf_m125_8TeV","ggh_m125_8TeV"])[0] ## ,"ggh_m125_8TeV"
    ## sig = histos[0]
    ## bkg = histos[1].Clone()
    ## ### ## bkg.Scale(1./80.)
    ## for b in histos[2:]:
    ##     bkg.Add(b)

    ### ###########################################################################################################
    ### Toy model:
    ###  - S & B gaussian in discriminant dimension
    ###  - B exponential in mass
    ###  - S gaussian in mass
    nexpSig = 300.
    soB = 1.e-1
    sigToy = ROOT.TF2("sigToy","[0]*exp(-0.5*((x-[1])/[2])**2)*exp(-0.5*((y-[3])/([4]+[5]*(x-[1])**2))**2)",minX,maxX,minMass,maxMass)
    ### sigToy.SetParameters(1.,0.7*(minX+maxX),0.1*(minX-maxX),125.,0.016*125.)
    sigToy.SetParameters(1.,0.7,0.1,125.,0.016*125.,+0.01*125.)
    sigToy.SetParameter(0, 1./sigToy.Integral(minNormX, maxNormX,minMass,maxMass) )

    bkgToy = ROOT.TF2("bkgToy","[0]*exp(-0.5*((x-[1])/[2])**2)*exp(-y/[3])",minX,maxX,minMass,maxMass)
    ### bkgToy.SetParameters(1.,0.3*(minX+maxX),0.5*(minX-maxX),(maxMass-minMass)/log(10.))
    bkgToy.SetParameters(1.,0.3,0.5,(maxMass-minMass)/log(10.))
    bkgToy.SetParameter(0, 1./bkgToy.Integral(minNormX, maxNormX,minMass,maxMass) )
    
    ### bkgToy = ROOT.TF2("bkgToy","[0]*exp(-0.5*((x-[1])/[2])**2)*pow(y,-[3])",minX,maxX,minMass,maxMass)
    ### ### bkgToy.SetParameters(1.,0.3*(minX+maxX),0.5*(minX-maxX),4.)
    ### bkgToy.SetParameters(1.,0.3,0.5,4.)
    ### bkgToy.SetParameter(0, 1./bkgToy.Integral(minNormX, maxNormX,minMass,maxMass) )

    soB *= bkgToy.Eval(0.7,125.)/sigToy.Eval(0.7,125.)
    
    ### Fill 2D histogram from toy model
    histos = [ ROOT.TH2F("hsig","hsig",nbinsX+1,minX-0.5*(maxX-minX)/nbinsX,maxX+0.5*(maxX-minX)/nbinsX,nMassBins,minMass,maxMass),
               ROOT.TH2F("hbkg","hbkg",nbinsX+1,minX-0.5*(maxX-minX)/nbinsX,maxX+0.5*(maxX-minX)/nbinsX,nMassBins,minMass,maxMass) ]
    histos[0].FillRandom( "sigToy", 100000 )
    histos[0].Scale( nexpSig/histos[0].Integral() )
    histos[1].FillRandom( "bkgToy", 100000 )
    histos[1].Scale( nexpSig/(soB*histos[1].Integral()) )

    ### Extract N, Sum X and Sum X^2 (X=mass) in bins of the discriminant
    sig = histos[0]
    bkg = histos[1].Clone()
    for b in histos[2:]:
        bkgN.Add(b)
    bkgN,bkgX,bkgX2 = getSecondOrder(bkg)
    sigN,sigX,sigX2 = getSecondOrder(sig)

    ### Set up workspace
    ws = ROOT.RooWorkspace("ws","ws")
    ws.rooImport = getattr(ws,'import')
    
    mva = ws.factory("mva[%f,%f.,%f]" % (minX,minX,maxX) )
    mass = ws.factory("mass[%f,%f.,%f]" % (125.,minMass,maxMass) )
    mass.setBins(nMassBins)
    xvar = ROOT.RooArgList(mva,mass)
    
    mu = ws.factory("mu[%f,%f.,%f]" % (1.,0.,10.) )
    
    sigDset = ROOT.RooDataHist("sigDset","sigDset",xvar,sig)
    bkgDset = ROOT.RooDataHist("bkgDset","bkgDset",xvar,bkg)

    ws.rooImport(sigDset)
    ws.rooImport(bkgDset)
    
    hsig = sigDset.createHistogram("hsig",mva)
    hbkg = bkgDset.createHistogram("hbkg",mva)
    ### hbkg.Smooth(30)

    hsigMass = sigDset.createHistogram("hsigMass",mass)
    hbkgMass = bkgDset.createHistogram("hbkgMass",mass)

    ### ##########################################################################################################
    ### Model builders
    ###

    ### Simple counting
    ## sigModel  = ROOT.CutAndCountModelBuilder( ROOT.AbsModel.sig, sigN, sigN.Integral(), minX, maxX )
    ## bkgModel  = ROOT.CutAndCountModelBuilder( ROOT.AbsModel.bkg, bkgN, bkgN.Integral(), minX, maxX )

    ### Simply parametrized shapes (gaus and expo)
    sigModel  = ROOT.SecondOrderModelBuilder( ROOT.AbsModel.sig, "sig", mass, sigN, sigX, sigX2, hsig.Integral(), minX, maxX )
    bkgModel  = ROOT.SecondOrderModelBuilder( ROOT.AbsModel.bkg, "bkg", mass, bkgN, bkgX, bkgX2, hbkg.Integral(), minX, maxX )
    
    ### #########################################################################################################
    ### Figure of merit for optimization
    ###

    ### Simple counting
    fom       = ROOT.NaiveCutAndCountFomProvider()
    ## fom       = ROOT.PoissonCutAndCountFomProvider()

    ### ### Likelihood ratio using asymptotic approx.
    ### fom       = ROOT.SimpleShapeFomProvider()
    ### sigModel.getModel().setMu(mu)
    ### fom.addPOI(mu)
    ### fom.minStrategy(1)
    ### ## fom.minimizer("Minuit2")
    ### ## fom.useRooSimultaneous()
    
    ### #########################################################################################################
    ### Run optimization
    ###
    minimizer = ROOT.TMinuitMinimizer()
    ### minimizer = ROOT.Minuit2.Minuit2Minimizer()
    optimizer = ROOT.CategoryOptimizer( minimizer, ndim )
    
    optimizer.addSignal( sigModel, True )
    ## optimizer.addSignal( sigModel )
    optimizer.addBackground( bkgModel )
    optimizer.setFigureOfMerit( fom )
    ## optimizer.addConstraint(True,50.,False) ## Add penalty term to likelihood to force the left-most boundary position
    optimizer.floatFirst() ## Float right-most boundary
    ## optimizer.refitLast()  ## Refit left-most boundary
    optimizer.absoluteBoundaries()  ## Float absolut boundaries instead of telescopic ones

    ## summary = optmizeCats( optimizer, ws, range(1,20), (numpy.array([cutoff]),False,True, ), True )
    summary = optmizeCats( optimizer, ws, range(1,6), (numpy.array([cutoff]),False,True,), False, False)
    ## summary = optmizeCats( optimizer, ws, [1], (numpy.array([cutoff]),False,True,), False, False)

    ### #########################################################################################################
    ### Some plots
    ###
    
    hbound = ROOT.TH2F("hbound","hbound",21,-0.5,20.5,nbinsX+1,minX-0.5*(maxX-minX)/nbinsX,maxX+0.5*(maxX-minX)/nbinsX)
    grS = ROOT.TGraph()
    grS.SetName("zVsNcat")
    grS.SetTitle(";n_{cat};f.o.m [A.U.]")
    for ncat,val in summary.iteritems():
        for b in val["boundaries"]:
            bd = float(b)
            if bd < 0.: bd = 0.
            hbound.Fill(float(ncat),bd)
        if( val["fom"] < 0. ) :
            grS.SetPoint( grS.GetN(), float(ncat), -val["fom"] )


    canv1 = ROOT.TCanvas("canv1","canv1")
    canv1.cd()
    hbkg.SetLineColor(ROOT.kRed)
    hsig.SetLineColor(ROOT.kBlue)
    hsig.DrawNormalized("hist")
    hbkg.DrawNormalized("hist SAME")
    canv1.SaveAs("cat_opt_disc.png")
    ws.rooImport(canv1)
    
    sigTF1    = ROOT.TF1("sigTF1",sigModel.getInputModelN(),minX,maxX,0,"HistoToTF1")
    bkgTF1    = ROOT.TF1("bkgTF1",bkgModel.getInputModelN(),minX,maxX,0,"HistoToTF1")
    canv2 = ROOT.TCanvas("canv2","canv2")
    canv2.cd()
    bkgTF1.SetLineColor(ROOT.kRed)
    sigTF1.SetLineColor(ROOT.kBlue)
    bkgTF1.Draw("")
    sigTF1.Draw("SAME")
    canv2.SaveAs("cat_opt_cdf.png")
    ws.rooImport(canv2)

    canv3 = ROOT.TCanvas("canv3","canv3")
    canv3.cd()
    hbkgMass.SetLineColor(ROOT.kRed)
    hsigMass.SetLineColor(ROOT.kBlue)
    hbkgMass.Draw("hist")
    hsigMass.Draw("hist SAME")
    canv3.SaveAs("cat_opt_mass.png")
    ws.rooImport(canv3)
    
    canv4 = ROOT.TCanvas("canv4","canv4")
    canv4.cd()
    bkgX.SetLineColor(ROOT.kRed)
    sigX.SetLineColor(ROOT.kBlue)
    bkgX.Draw("hist")
    sigX.Draw("hist SAME")
    canv4.SaveAs("cat_opt_sum_mass.png")
    ws.rooImport(canv4)

    canv5 = ROOT.TCanvas("canv5","canv5")
    canv5.cd()
    bkgX2.SetLineColor(ROOT.kRed)
    sigX2.SetLineColor(ROOT.kBlue)
    bkgX2.Draw("")
    sigX2.Draw("SAME")
    canv5.SaveAs("cat_opt_sum_mass2.png")
    ws.rooImport(canv5)

    canv6 = ROOT.TCanvas("canv6","canv6")
    canv6.cd()
    bkgN.SetLineColor(ROOT.kRed)
    sigN.SetLineColor(ROOT.kBlue)
    bkgN.Draw("")
    sigN.Draw("SAME")
    canv6.SaveAs("cat_opt_n.png")
    ws.rooImport(canv6)

    canv7 = ROOT.TCanvas("canv7","canv7")
    canv7.SetGridx()
    canv7.SetGridy()
    canv7.cd()
    hbound.Draw("box")
    canv7.SaveAs("cat_opt_bound.png")
    canv7.SaveAs("cat_opt_bound.C")
    ws.rooImport(canv7)

    canv8 = ROOT.TCanvas("canv8","canv8")
    canv8.SetGridx()
    canv8.SetGridy()
    canv8.cd()
    hbound.ProjectionY().Draw("")
    canv8.SaveAs("cat_opt_bound_pj.png")
    canv8.SaveAs("cat_opt_bound_pj.C")
    ws.rooImport(canv8)
    
    canv9 = ROOT.TCanvas("canv9","canv9")
    canv9.SetGridx()
    canv9.SetGridy()
    canv9.cd()
    grS.SetMarkerStyle(ROOT.kFullCircle)
    grS.Draw("AP")

    canv9.SaveAs("cat_opt_fom.png")
    canv9.SaveAs("cat_opt_fom.C")

    ws.rooImport(canv9)
    ws.rooImport(grS)
    
    
    
    return ws

# -----------------------------------------------------------------------------------------------------------
def main(fin=None):

    if fin:
        fin = ROOT.TFile.Open(fname)
    ROOT.gSystem.SetIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include")
    ROOT.gSystem.Load("../../libLoopAll")

    ROOT.gStyle.SetOptStat(0)

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.RooMsgService.instance().setSilentMode(True)
    ws = optimize1D(fin)
    
    return ws

    
if __name__ == "__main__":
    ## ws=main(sys.argv[1])
    ws=main()
