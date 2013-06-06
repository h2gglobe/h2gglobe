#!/bin/env python

import sys
import ROOT
import numpy
from math import sqrt, log
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '.3f')

objs = []

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


alltrees = []

def mergeTrees(tfile,outname,trees):
    tlist = ROOT.TList()
    for name,selection in trees:
        print "Reading tree %s" % name        
        tree=tfile.Get(name)
        if selection != "":
            clone = tree.CopyTree(selection)
            ### tree.GetListOfClones().Remove(clone)
            ### tree.ResetBranchAddresses()
            ### clone.ResetBranchAddresses()
            tree = clone
        tlist.Add(tree)
    out=ROOT.TTree.MergeTrees(tlist)
    out.SetName(outname)
    return out

# -----------------------------------------------------------------------------------------------------------
def optimize1D(fin):
    minX = -1.
    maxX = 1.
    nbinsX = int((maxX-minX)/0.01)
    cutoff = 0.05
    ndim = 1

    minMass = 100.
    maxMass = 180.
    nMassBins = 160
    ## nMassBins = 320

    signals = {
        "sigRv" :[ ("vbfAnalysis/vbf_m125_8TeV","corrVeretx"),
                   ("vbfAnalysis/ggh_m125_8TeV","corrVeretx"),
                   ],
        "sigWv" :[ ("vbfAnalysis/vbf_m125_8TeV","!corrVeretx"),
                   ("vbfAnalysis/ggh_m125_8TeV","!corrVeretx"),
                   ],
        ### "sigGGH" :[("vbfAnalysis/ggh_m125_8TeV","1"),
        ###            ],
        ### "sigQQH" :[("vbfAnalysis/vbf_m125_8TeV","1"),
        ###            ],
        ### "sig" :[("vbfAnalysis/ggh_m125_8TeV","1"),
        ###         ("vbfAnalysis/vbf_m125_8TeV","1"),
        ###         ]
        }
    backgrounds = [
        ("vbfAnalysis/diphojet_Box_8TeV","1"),
        ## ("vbfAnalysis/diphojet_8TeV","1"),
        ## ("vbfAnalysis/dipho_Box_25_8TeV","1"),
        ## ("vbfAnalysis/dipho_Box_250_8TeV","1"),
        ("vbfAnalysis/gjet_20_8TeV_pf","1"),
        ("vbfAnalysis/gjet_40_8TeV_pf","1"),
        ("vbfAnalysis/qcd_40_8TeV_S7_ff","1"),
        ]

    tmp = ROOT.TFile.Open("categoryOptimization.root","recreate")
    tmp.cd()
    sigTrees = [ mergeTrees(fin,name,trees) for name,trees in signals.iteritems() ]
    bkgTree = mergeTrees(fin,"background",backgrounds)
    
    bkgTree.Write("background")
    bkgTree.Draw("diphoM:diphoMVA>>bkg(%d,%f,%f,%d,%f,%f)" % (nbinsX,minX,maxX,nMassBins,minMass,maxMass), "weight", "goff" )
    bkg = ROOT.gDirectory.Get("bkg")
    bkg.SetDirectory(0)

    signals = []
    for sigTree in sigTrees:
        modelName = "%sModel" % sigTree.GetName()
        sigTree.Write()
        sigTree.Draw("diphoM:diphoMVA>>%s(%d,%f,%f,%d,%f,%f)" % (modelName,
                                                                 nbinsX,minX,maxX,nMassBins,minMass,maxMass), "weight", "goff" )
        sig = ROOT.gDirectory.Get(modelName)
        sig.SetDirectory(0)
        signals.append( (sig,modelName) )
    
    fin.Close()
    tmp.Close()
    
    ### Set up workspace
    ws = None
    ### ws = ROOT.RooWorkspace("ws","ws")
    ### ws.rooImport = getattr(ws,'import')
    ### 
    ### mva = ws.factory("diphoMVA[%f,%f.,%f]" % (minX,minX,maxX) )
    ### ## mass = 5Bws.factory("mass[%f,%f.,%f]" % (125.,minMass,maxMass) )
    ### mass = ws.factory("diphoM[%f,%f.,%f]" % (125.,minMass,maxMass) )
    mva = ROOT.RooRealVar("diphoMVA","diphoMVA",minX,minX,maxX)
    mass = ROOT.RooRealVar("diphoM","diphoM",125.,minMass,maxMass)
    mass.setBins(nMassBins)

    ## mu = ws.factory("mu[%f,%f.,%f]" % (1.,0.,10.) )
    mu = ROOT.RooRealVar("mu","mu",1.,0.,10.)

    bkgN,bkgX,bkgX2 = getSecondOrder(bkg)
    sigModels = []
    hsigs = []
    hsigsMass = []
    sigTF1s = []
    for sig,name in signals:
        sigN,sigX,sigX2 = getSecondOrder(sig)
        sigModel = ROOT.SecondOrderModelBuilder( ROOT.AbsModel.sig, name, mass, sigN, sigX, sigX2, sig.Integral(), minX, maxX ) 
        sigModels.append( sigModel )
        hsigs.append( sig.ProjectionX() )
        hsigsMass.append( sig.ProjectionY() )
        sigTF1s.append( ROOT.TF1("sigTF1%s" % name,sigModel.getInputModelN(),minX,maxX,0,"HistoToTF1") )

    ## hsig = sig.ProjectionX()
    hbkg = bkg.ProjectionX()
    
    ## hsigMass = sig.ProjectionY()
    hbkgMass = bkg.ProjectionY()

    ### ##########################################################################################################
    ### Model builders
    ###

    ### Simple counting
    ## sigModel  = ROOT.CutAndCountModelBuilder( ROOT.AbsModel.sig, sigN, sigN.Integral(), minX, maxX )
    ## bkgModel  = ROOT.CutAndCountModelBuilder( ROOT.AbsModel.bkg, bkgN, bkgN.Integral(), minX, maxX )

    ### Simply parametrized shapes (gaus and expo)
    bkgModel  = ROOT.SecondOrderModelBuilder( ROOT.AbsModel.bkg, "bkgModel", mass, bkgN, bkgX, bkgX2, hbkg.Integral(), minX, maxX )
    ## sigModel  = ROOT.SecondOrderModelBuilder( ROOT.AbsModel.sig, "sigModel", mass, sigN, sigX, sigX2, hsig.Integral(), minX, maxX )
    ### sigModels = [ ROOT.SecondOrderModelBuilder( ROOT.AbsModel.sig, "sigModel", mass, sigN, sigX, sigX2, hsig.Integral(), minX, maxX )
    ###               for sigN,sigX,sigX2 in sigMoments ]
    
    ### #########################################################################################################
    ### Figure of merit for optimization
    ###

    ### Simple counting
    fom       = ROOT.NaiveCutAndCountFomProvider()
    ## fom       = ROOT.PoissonCutAndCountFomProvider()

    ### ### Likelihood ratio using asymptotic approx.
    ## fom       = ROOT.SimpleShapeFomProvider()
    ## for sigModel in sigModels:
    ##     sigModel.getModel().setMu(mu)
    ## fom.addPOI(mu)
    ## fom.minStrategy(1)
    ## ## fom.minimizer("Minuit2")
    ## ## fom.useRooSimultaneous()
    
    ### #########################################################################################################
    ### Run optimization
    ###
    minimizer = ROOT.TMinuitMinimizer("Minimize")
    ## minimizer.SetPrintLevel(999)
    ## minimizer = ROOT.Minuit2.Minuit2Minimizer()
    optimizer = ROOT.CategoryOptimizer( minimizer, ndim )

    for sigModel in sigModels:
        optimizer.addSignal( sigModel, True )
    optimizer.addBackground( bkgModel )
    optimizer.setFigureOfMerit( fom )
    ## optimizer.addConstraint(True,50.,False) ## Add penalty term to likelihood to force the left-most boundary position
    ## optimizer.floatFirst() ## Float right-most boundary
    optimizer.refitLast()  ## Refit left-most boundary
    optimizer.absoluteBoundaries()  ## Float absolut boundaries instead of telescopic ones
    

    ## summary = optmizeCats( optimizer, ws, range(1,20), (numpy.array([cutoff]),False,True, ), True )
    summary = optmizeCats( optimizer, ws, range(1,8), (numpy.array([cutoff]),False,True,), False, True )
    ## summary = optmizeCats( optimizer, ws, [1], (numpy.array([cutoff]),False,True,), False, False)
    
    ### #########################################################################################################
    ### Some plots
    ###
    
    grS = ROOT.TGraph()
    grS.SetName("zVsNcat")
    grS.SetTitle(";n_{cat};f.o.m [A.U.]")
    for ncat,val in summary.iteritems():
        if( val["fom"] < 0. ) :
            grS.SetPoint( grS.GetN(), float(ncat), -val["fom"] )
    grS.Sort()
    mincat = grS.GetX()[0]
    maxcat = grS.GetX()[grS.GetN()-1]
    ncat = int(maxcat - mincat)
    
    hbound = ROOT.TH2F("hbound","hbound",ncat+3,mincat-1.5,maxcat+1.5,nbinsX+3,minX-1.5*(maxX-minX)/nbinsX,maxX+1.5*(maxX-minX)/nbinsX)
    for ncat,val in summary.iteritems():
        for b in val["boundaries"]:
            bd = float(b)
            if bd < minX: bd = minX
            hbound.Fill(float(ncat),bd)
    

    canv1 = ROOT.TCanvas("canv1","canv1")
    canv1.cd()
    hbkg.SetLineColor(ROOT.kRed)
    hbkg.DrawNormalized("hist")
    for hsig in hsigs:
        hsig.SetLineColor(ROOT.kBlue)
        hsig.DrawNormalized("hist SAME")
    canv1.SaveAs("cat_opt_disc.png")

    ## sigTF1    = ROOT.TF1("sigTF1",sigModel.getInputModelN(),minX,maxX,0,"HistoToTF1")
    bkgTF1    = ROOT.TF1("bkgTF1",bkgModel.getInputModelN(),minX,maxX,0,"HistoToTF1")
    canv2 = ROOT.TCanvas("canv2","canv2")
    canv2.cd()
    bkgTF1.SetLineColor(ROOT.kRed)
    bkgTF1.Draw("")
    ## sigTF1.Draw("SAME")
    for sigTF1 in sigTF1s:
        sigTF1.SetLineColor(ROOT.kBlue)
        sigTF1.Draw("SAME")
    canv2.SaveAs("cat_opt_cdf.png")

    canv3 = ROOT.TCanvas("canv3","canv3")
    canv3.cd()
    hbkgMass.SetLineColor(ROOT.kRed)
    hbkgMass.Draw("hist")
    for hsigMass in hsigsMass:
        hsigMass.SetLineColor(ROOT.kBlue)
        hsigMass.Draw("hist SAME")
    canv3.SaveAs("cat_opt_mass.png")
    
    canv4 = ROOT.TCanvas("canv4","canv4")
    canv4.cd()
    bkgX.SetLineColor(ROOT.kRed)
    sigX.SetLineColor(ROOT.kBlue)
    bkgX.Draw("hist")
    sigX.Draw("hist SAME")
    canv4.SaveAs("cat_opt_sum_mass.png")

    canv5 = ROOT.TCanvas("canv5","canv5")
    canv5.cd()
    bkgX2.SetLineColor(ROOT.kRed)
    sigX2.SetLineColor(ROOT.kBlue)
    bkgX2.Draw("")
    sigX2.Draw("SAME")
    canv5.SaveAs("cat_opt_sum_mass2.png")

    canv6 = ROOT.TCanvas("canv6","canv6")
    canv6.cd()
    bkgN.SetLineColor(ROOT.kRed)
    sigN.SetLineColor(ROOT.kBlue)
    bkgN.Draw("")
    sigN.Draw("SAME")
    canv6.SaveAs("cat_opt_n.png")

    canv7 = ROOT.TCanvas("canv7","canv7")
    canv7.SetGridx()
    canv7.SetGridy()
    canv7.cd()
    hbound.Draw("box")
    canv7.SaveAs("cat_opt_bound.png")
    canv7.SaveAs("cat_opt_bound.C")

    canv8 = ROOT.TCanvas("canv8","canv8")
    canv8.SetGridx()
    canv8.SetGridy()
    canv8.cd()
    hbound.ProjectionY().Draw("")
    canv8.SaveAs("cat_opt_bound_pj.png")
    canv8.SaveAs("cat_opt_bound_pj.C")
    
    canv9 = ROOT.TCanvas("canv9","canv9")
    canv9.SetGridx()
    canv9.SetGridy()
    canv9.cd()
    grS.SetMarkerStyle(ROOT.kFullCircle)
    grS.Draw("AP")

    canv9.SaveAs("cat_opt_fom.png")
    canv9.SaveAs("cat_opt_fom.C")

    return ws

# -----------------------------------------------------------------------------------------------------------
def main(fin=None):

    if fin:
        fin = ROOT.TFile.Open(fin)
    ROOT.gSystem.SetIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include")
    ROOT.gSystem.Load("../../libLoopAll.so")

    ROOT.gStyle.SetOptStat(0)

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.RooMsgService.instance().setSilentMode(True)
    ws = optimize1D(fin)
    
    return ws

    
if __name__ == "__main__":
    ws=main(sys.argv[1])
    ## ws=main()
