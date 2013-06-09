#!/bin/env python

import sys
import ROOT
import numpy
from math import sqrt, log
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '.3f')

from optparse import OptionParser, make_option

objs = []

## ------------------------------------------------------------------------------------------------------------------------------------------------------
def loadSettings(cfgs,dest):
    for cfg in cfgs.split(","):
        cf = open(cfg)
        settings = json.load(cf)
        for k,v in settings.iteritems():
            setattr(dest,k,v)
        cf.close()

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

# -----------------------------------------------------------------------------------------------------------
def mergeTrees(tfile,sel,outname,trees,aliases):
    tlist = ROOT.TList()
    for name,selection in trees:
        print "Reading tree %s" % name        
        tree=tfile.Get(name)
        if sel != "":
            selection = str(TCut(selection)*TCut(sel))
        if selection != "":
            clone = tree.CopyTree(selection)
            tree = clone
        tlist.Add(tree)
    out=ROOT.TTree.MergeTrees(tlist)
    out.SetName(outname)
    for name, definition in aliases:
        out.SetAlias( name,definition )
    return out

objs = []

# -----------------------------------------------------------------------------------------------------------
def defineVariables(variables):

    arglist = ROOT.RooArgList()
    aliases = []
    
    for var in variables:
        name = str(var[0])
        if ":=" in name:
            name, definition = [ tok.lstrip(" ").rstrip(" ") for tok in name.split(":=") if tok != "" ]
            print name, definition
            aliases.append( (name,definition) )
        xmin,xmax,nbins = var[1]
        default = xmin
        if len(var) == 3:
            default = float(var[2])
        if type(nbins) == float:
            nbins = int( (xmax-xmin)/nbins )
        var = ROOT.RooRealVar(name,name,default,xmin,xmax)
        var.setBins(nbins)
        objs.append(var)
        arglist.add( var )

    return arglist,aliases
        
# -----------------------------------------------------------------------------------------------------------
def modelBuilders(trees, type, obs, varlist, weight):
    builders=[]
    for tree in trees:
        modelName = "%sModel" % tree.GetName()
        modelBuilder = ROOT.SecondOrderModelBuilder(type, modelName, obs, tree, varlist, weight)
        builders.append(modelBuilder)
    return builders

# -----------------------------------------------------------------------------------------------------------
def optimizeMultiDim(options,args):


    signals = options.signals
    backgrounds = options.backgrounds

    variables = options.variables
    observable = options.observable
    selection = options.selection

    obs,obsalias = defineVariables( [observable] )
    obs = obs[0]
    
    varlist,aliases = defineVariables( variables )

    print "Observables "
    obs.Print("")

    print "Variables"
    varlist.Print("V")
    
    aliases.extend(obsalias)
    print "Aliases"
    print aliases

    if options.infile == "":
        options.infile = args[0]
    fin = ROOT.TFile.Open(options.infile)

    tmp = ROOT.TFile.Open(options.outfile,"recreate")
    tmp.cd()

    ## sigTrees = [ mergeTrees(fin,selection,name,trees,aliases) for name,trees in signals.iteritems() ]
    bkgTrees = [ mergeTrees(fin,selection,name,trees,aliases) for name,trees in backgrounds.iteritems() ]

    fin.Close()
    ## tmp.Close()
    
    ### ##########################################################################################################
    ### Model builders
    ###
    signals = [] ## modelBuilders( sigTrees, ROOT.AbsModel.sig, obs, varlist, "weight" )
    backgrounds = modelBuilders( bkgTrees, ROOT.AbsModel.bkg, obs, varlist, "weight" )

    for model in signals+backgrounds:
        model.getTree().Scan()
        
    return None
    
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

    ## tmp.Close()
    return ws

# -----------------------------------------------------------------------------------------------------------
def main(options,args):
    
    ROOT.gSystem.SetIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include")
    ROOT.gSystem.Load("../../libLoopAll.so")

    ROOT.gStyle.SetOptStat(0)

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.RooMsgService.instance().setSilentMode(True)
    ws = optimizeMultiDim(options,args)
    
    return ws

    
if __name__ == "__main__":

    parser = OptionParser(option_list=[
            make_option("-i", "--infile",
                        action="store", type="string", dest="infile",
                        default="",
                        help="input file",
                        ),
            make_option("-o", "--outfile",
                        action="store", type="string", dest="outfile",
                        default=sys.argv[0].replace(".py",".root"),
                        help="output file",
                        ),
            make_option("-l", "--load",
                        action="store", dest="settings", type="string",
                        default="",
                        help="json file containing settings"
                        ),
            ])

    (options, args) = parser.parse_args()
    loadSettings(options.settings, options)
    
    
    ws=main(options,args)
