#!/bin/env python

import sys, types, os
import numpy
from math import sqrt, log
import json
json.encoder.FLOAT_REPR = lambda o: format(o, '.3f')

from optparse import OptionParser, make_option
from  pprint import pprint

objs = []

# -----------------------------------------------------------------------------------------------------------
def loadSettings(cfgs,dest):
    for cfg in cfgs.split(","):
        cf = open(cfg)
        settings = json.load(cf)
        for k,v in settings.iteritems():
            attr  = getattr(dest,k,None)
            if attr and type(attr) == list:
                attr.extend(v)
                v = attr
            setattr(dest,k,v)
        cf.close()

def getBoundaries(ndim,ncat,optimizer, summary):
    print ndim, ncat
    boundaries = numpy.array([0. for i in range(ndim*(ncat+1)) ],dtype='d')
    selections = numpy.array([0. for i in range(optimizer.nOrthoCuts()) ],dtype='d')
    print boundaries
    if len(selections) > 0: 
        z = optimizer.getBoundaries(ncat,boundaries,selections)
    else:
        tmp = numpy.array([0.],dtype='d')
        z = optimizer.getBoundaries(ncat,boundaries,tmp)
    if not ncat in summary or summary[ncat]["fom"] < z: 
        summary[ncat] =  { "fom" : z, "boundaries" : list(boundaries), "ncat": ncat }
        if len(selections) > 0: 
            summary[ncat]["selections"]=list(selections)
    objs.append((boundaries,selections))

# -----------------------------------------------------------------------------------------------------------
def optmizeCats(optimizer,ws,ndim,rng,args,readBack=False,reduce=False,refit=0):
    
    summary = {}
    print readBack
    if readBack:
        try:
            sin = open("%s/cat_opt.json" % options.cont,"r")
            summary = json.loads(sin.read())
            sin.close()
        except:
            summary = {}
            
    print "---------------------------------------------"
    print "Fitting"
    print 
    for iter in rng:
        optimizer.optimizeNCat(iter,*args)
        getBoundaries(ndim,iter, optimizer, summary )

    for ncat,val in summary.iteritems():
        printBoundaries(ndim,val["boundaries"],val["fom"],val.get("selections",None))
        
    if reduce:
        print "---------------------------------------------"
        print "Reducing"
        print 
        maxncat = 0
        for ncat,val in summary.iteritems():
            boundaries = numpy.array([float(b) for b in val["boundaries"]])
            if optimizer.nOrthoCuts() > 0:
                if "selections" in val:
                    for isel in range(len(val["selections"])):
                        optimizer.setOrthoCut(isel, float(val["selections"][isel]))
            ## setSelections(optimizer,summary)
            optimizer.reduce(int(ncat)+1, boundaries, numpy.array([0. for i in range(len(boundaries))]) )
            maxncat = max(int(ncat),maxncat)
        rng = range(1,maxncat+1)
        summary = {}
        for iter in rng:
            getBoundaries(ndim,iter, optimizer, summary )

    if refit > 0:
        print "---------------------------------------------"
        print "Refitting"
        print
        for irefit in range(refit):
            for ncat,val in summary.iteritems():
                boundaries = numpy.array([float(b) for b in val["boundaries"]])
                if optimizer.nOrthoCuts() > 0:
                    if "selections" in val:
                        for isel in range(len(val["selections"])):
                            optimizer.setOrthoCut(isel, float(val["selections"][isel]))
                print ncat,val,boundaries
                ncat=int(ncat)
                ## setSelections(optimizer,summary)
                eargs = [a for a in args]
                eargs.append(boundaries)
                rng.append(ncat)
                optimizer.optimizeNCat(ncat,*eargs)
            rng = sorted( set(rng) )
            summary = {}
        
            for iter in rng:
                getBoundaries(ndim,iter, optimizer, summary )
        
    for ncat,val in summary.iteritems():
        printBoundaries(ndim,val["boundaries"],val["fom"],val.get("selections",None))

    sout = open("cat_opt.json","w+")
    sout.write( json.dumps(summary,sort_keys=True, indent=4) )

    return summary

# -----------------------------------------------------------------------------------------------------------
def printBoundaries(ndim,boundaries, maxval,selections):
    
    print "---------------------------------------------"
    print "ncat: ", len(boundaries)/ndim-1
    print "max: %1.5f" % ( maxval )
    print "boundaries: ",
    for b in boundaries:
        print "%1.3g" % b,
    print
    if selections:
        print "selections: ",
        for s in selections:
            print "%1.3g" % s,
        print
    print 


alltrees = []

# -----------------------------------------------------------------------------------------------------------
def mergeTrees(tfile,sel,outname,trees,aliases):
    tlist = ROOT.TList()
    for name,selection in trees:
        tree=tfile.Get(name)
        print "Reading tree %s" % name, selection, sel
        if sel != "":
            selection = str(ROOT.TCut(selection)*ROOT.TCut(sel))
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
def defineVariables(variables,label):

    arglist = ROOT.RooArgList()
    aliases = []
    
    for var in variables:
        name = str(var[0] % {"label" : label})
        if ":=" in name:
            name, definition = [ tok.lstrip(" ").rstrip(" ") for tok in name.split(":=") if tok != "" ]
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
def modelBuilders(trees, type, obs, varlist, sellist, weights, shapes, minevents, constrained):
    builders=[]
    constraints=[]
    for tree in trees:
        name = tree.GetName()
        modelName = "%sModel" % name
        weight = "weight"
        if name in weights:
            weight = weights[name]
        modelBuilder = ROOT.SecondOrderModelBuilder(type, modelName, obs, tree, varlist, sellist, weight)
        if name in shapes:
            modelBuilder.getModel().setShape( getattr(ROOT.SecondOrderModel,shapes[name]) )
        if name in minevents:
            modelBuilder.getModel().minEvents( minevents[name] )
        if name in constrained:
            cname = "normConstraint%s" % name
            constraint = ROOT.RooRealVar(cname,cname,1.)
            pdf = None
            if constrained[name] > 0.:
                ### k = ROOT.RooRealVar("kLn%s" % cname,"kLn%s" % cname,math.exp(constrained[name]))
                ### k.setConstant(True)
                ### pdf = ROOT.RooLogNormal("%sPdf" % cname,"%sPdf" % cname,constraint,k)
                mean = ROOT.RooRealVar("mean%s" % cname,"mean%s" % cname,1.)
                mean.setConstant(True)
                sigma = ROOT.RooRealVar("sigma%s" % cname,"sigma%s" % cname,constrained[name])
                sigma.setConstant(True)
                constraint.Print("V")
                pdf = ROOT.RooGaussian("%sPdf" % cname,"%sPdf" % cname,constraint,mean,sigma)
                constraints.append( (constraint,pdf,(mean,sigma)) )
            else:
                constraint.setConstant(True)
                constraints.append(constraint)
            modelBuilder.getModel().setMu(constraint)
            
        builders.append(modelBuilder)
    return builders,constraints

# -----------------------------------------------------------------------------------------------------------
def optimizeMultiDim(options,args):


    signals = options.signals
    backgrounds = options.backgrounds

    variables = options.variables
    observable = options.observable
    selection = options.selection
    cutoff = options.cutoff

    selectioncuts = options.selectioncuts
    
    ndim = len(variables)
    ws = None

    cutoffs = numpy.array([cutoff]*ndim)
    
    obs,obsalias = defineVariables( [observable], options.label )
    obs = obs[0]
    mu = ROOT.RooRealVar("mu","mu",1.,0.,10.)
    
    varlist,aliases = defineVariables( variables, options.label )
    sellist,selaliases = defineVariables( selectioncuts, options.label )

    print "---------------------------------------------"
    print "Observables "
    obs.Print("")

    print 
    print "Variables"
    varlist.Print("V")

    print
    print "Selection cuts"
    sellist.Print("V")

    aliases.extend(obsalias+selaliases)
    print
    print "Aliases"
    pprint(aliases)
    print
    
    if options.infile == "":
        options.infile = args[0]
    fin = ROOT.TFile.Open(options.infile)

    if options.cont:
        if os.path.exists(os.path.abspath(options.cont)):
            options.cont = os.path.abspath(options.cont)
        else:
            options.cont = os.path.abspath(os.path.join(os.path.dirname(options.outdir),options.cont))
    if not os.path.exists(options.outdir):
        os.mkdir(options.outdir)
    os.chdir(options.outdir)

    if not os.path.exists("/tmp/%s" % os.getlogin()):
        os.mkdir("/tmp/%s" % os.getlogin())
    tmp = ROOT.TFile.Open("/tmp/%s/categoryOptimizationMultiDim_%s.root" % (os.getlogin(),options.label) ,"recreate")
    tmp.cd()

    ### ##########################################################################################################
    ### Minimizer and optimizer
    ###
    minimizer = ROOT.TMinuitMinimizer("Minimize")
    ## minimizer.SetPrintLevel(999)
    ## minimizer = ROOT.Minuit2.Minuit2Minimizer()
    optimizer = ROOT.CategoryOptimizer( minimizer, ndim )
    for isel in range(sellist.getSize()):
        sel = sellist[isel]
        if sel.GetName() in options.fix:
            optimizer.addFixedOrthoCut(sel.GetName(),sel.getVal())
        else:
            optimizer.addFloatingOrthoCut(sel.GetName(),sel.getVal(),(sel.getMax()-sel.getMin())/sel.getBins(),sel.getMin(),sel.getMax())
    optimizer.absoluteBoundaries()  ## Float absolut boundaries instead of telescopic ones
    for opt in options.settings:
        if isinstance(opt, types.StringTypes):
            getattr(optimizer,opt)()
        else:
            name, args = opt
            getattr(optimizer,name)(*args)

    
    ### ##########################################################################################################
    ### Model builders
    ###
    print "---------------------------------------------"
    print "Reading inputs"
    print
    sigTrees = [ mergeTrees(fin,selection,name,trees,aliases) for name,trees in signals.iteritems() ]
    bkgTrees = [ mergeTrees(fin,selection,name,trees,aliases) for name,trees in backgrounds.iteritems() ]

    for tree in sigTrees+bkgTrees:
        tree.Write()
        
    fin.Close()
    if options.onlytrees:
        sys.exit(0)
    
    print "---------------------------------------------"
    print "Booking model buildes"
    print
    signals,sigconstr = modelBuilders( sigTrees, ROOT.AbsModel.sig, obs, varlist, sellist,
                                       getattr(options,"weights",{}), getattr(options,"shapes",{}), {}, {} )
    backgrounds,bkgconstr = modelBuilders( bkgTrees, ROOT.AbsModel.bkg, obs, varlist, sellist,
                                           getattr(options,"weights",{}), getattr(options,"shapes",{}),
                                           getattr(options,"minevents",{},), getattr(options,"constrained",{}),
                                           )
    
    if options.saveCompactTree:
        for model in signals+backgrounds:
            model.getTree().Write()
        
    normTF1s = []
    sumxTF1s = []
    sumx2TF1s = []
    for model in backgrounds:
        norm = model.getTF1N() 
        x    = model.getTF1X()
        x2   = model.getTF1X2() 
        norm.SetLineColor(ROOT.kRed)
        x.SetLineColor(ROOT.kRed)
        x2.SetLineColor(ROOT.kRed)
        normTF1s.append(norm)
        sumxTF1s.append(x)
        sumx2TF1s.append(x2)

    for model in signals:
        norm = model.getTF1N()
        x    = model.getTF1X()
        x2   = model.getTF1X2() 
        norm.SetLineColor(ROOT.kBlue)
        x.SetLineColor(ROOT.kBlue)
        x2.SetLineColor(ROOT.kBlue)
        normTF1s.append(norm)
        sumxTF1s.append(x)
        sumx2TF1s.append(x2)
        
    canv2 = ROOT.TCanvas("canv2","canv2")
    canv2.cd()
    normTF1s[0].Draw("")
    for tf1 in normTF1s[1:]:
        tf1.Draw("SAME")
    canv2.SaveAs("cat_opt_cdf.png")

    ### canv3 = ROOT.TCanvas("canv3","canv3")
    ### canv3.cd()
    ### hbkgMass.SetLineColor(ROOT.kRed)
    ### hbkgMass.Draw("hist")
    ### for hsigMass in hsigsMass:
    ###     hsigMass.SetLineColor(ROOT.kBlue)
    ###     hsigMass.Draw("hist SAME")
    ### canv3.SaveAs("cat_opt_mass.png")
    
    canv4 = ROOT.TCanvas("canv4","canv4")
    canv4.cd()
    sumxTF1s[0].Draw("hist")
    for tf1 in sumxTF1s[1:]:
        tf1.Draw("hist SAME")
    canv4.SaveAs("cat_opt_sum_mass.png")

    canv5 = ROOT.TCanvas("canv5","canv5")
    canv5.cd()
    sumx2TF1s[0].Draw("hist")
    for tf1 in sumx2TF1s[1:]:
        tf1.Draw("hist SAME")
    canv5.SaveAs("cat_opt_sum_mass2.png")

    objs.append( signals )
    objs.append( backgrounds )
    objs.append( normTF1s )
    objs.append( sumxTF1s )
    objs.append( sumx2TF1s )
    
    ### #########################################################################################################
    ### Figure of merit for optimization
    ###

    ### Simple counting
    ## fom       = ROOT.NaiveCutAndCountFomProvider()
    ## fom       = ROOT.PoissonCutAndCountFomProvider()

    ### ### Likelihood ratio using asymptotic approx.
    fom       = ROOT.SimpleShapeFomProvider()
    for sigModel in signals:
        sigModel.getModel().setMu(mu)
    fom.addPOI(mu)
    fom.minStrategy(2)
    for constraint in sigconstr+bkgconstr:
        print constraint
        if type(constraint) == tuple:
            var, pdf, pars = constraint
            print var,pdf,pars
            fom.addNuisance(var,pdf)
        else:
            fom.addNuisance(constraint)
    
    ## fom.minimizer("Minuit2")
    ## fom.useRooSimultaneous()
    
    ### #########################################################################################################
    ### Run optimization
    ###
        
    for sigModel in signals:
        optimizer.addSignal( sigModel, True )
        ## optimizer.addSignal( sigModel )
    for bkgModel in backgrounds:
        optimizer.addBackground( bkgModel )
    optimizer.setFigureOfMerit( fom )
    
    summary = optmizeCats( optimizer, ws, varlist.getSize(),
                           options.range, (cutoffs,options.dry,True,), options.cont, options.reduce, options.refit )
    
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

    for idim in range(ndim):
        var = varlist[idim]
        name = var.GetName()
        minX = var.getMin()
        maxX = var.getMax()
        nbinsX = var.getBinning().numBoundaries()-1
        hbound = ROOT.TH2F("hbound_%s" % name,"hbound_%s" % name,nbinsX+3,minX-1.5*(maxX-minX)/nbinsX,maxX+1.5*(maxX-minX)/nbinsX,ncat+3,mincat-1.5,maxcat+1.5)
        for jcat,val in summary.iteritems():
            nbound = int(jcat)+1
            for ib in range(nbound):
                bd = float(val["boundaries"][ib+idim*nbound])
                ## hbound.Fill(float(jcat),bd)
                hbound.Fill(bd,float(jcat))
        cbound = ROOT.TCanvas( "cat_opt_%s" % hbound.GetName(), "cat_opt_%s" % hbound.GetName() )
        cbound.cd()
        hbound.Draw("box")

        cbound_pj = ROOT.TCanvas( "cat_opt_%s_pj" % hbound.GetName(),  "cat_opt_%s_pj" %  hbound.GetName() )
        cbound_pj.cd()
        hbound_pj = hbound.Clone()
        ### hbound_pj = hbound.ProjectionY()
        ### hbound_pj.Draw()
        hbound_pj.Draw("box")
        hbound_pj.SetFillColor(ROOT.kBlack)
        hbound_pj.SetLineColor(ROOT.kBlack)
        objs.append(hbound)
        objs.append(hbound_pj)
        objs.append(cbound)
        objs.append(cbound_pj)
        maxy = 0.
        pdfs = []
        for sig in signals:
            pdf = sig.getPdf(idim)
            pdf.Scale(1./pdf.Integral())
            pdf.SetLineColor(ROOT.kBlue)
            pdfs.append(pdf)
            maxy = max(maxy,pdf.GetMaximum())
            objs.append(pdf)
        for bkg in backgrounds:
            pdf = bkg.getPdf(idim)
            pdf.Scale(1./pdf.Integral())
            pdf.SetLineColor(ROOT.kRed)
            pdfs.append(pdf)
            maxy = max(maxy,pdf.GetMaximum())
            objs.append(pdf)
            
        ### hbound_pj.Scale(maxy*hbound_pj.GetMaximum())
        ### hbound_pj.GetYaxis().SetRangeUser(0.,1.1*maxy)
        ### cbound_pj.RedrawAxis()
        ### cbound_pj.Update()
        for pdf in pdfs:
            pdf.Scale( (ncat+3.)/maxy )
            pdf.Draw("hist same")
        hbound_pj.GetYaxis().SetNdivisions(500+ncat+3)
        cbound_pj.SetGridy()
        hbound_pj.Draw("box same")
        
        for fmt in "png", "C":
            cbound.SaveAs("%s.%s" % (cbound.GetName(),fmt) )
            cbound_pj.SaveAs("%s.%s" % (cbound_pj.GetName(),fmt) )
            
    for isel in range(sellist.getSize()):
        var = sellist[isel]
        name = var.GetName()
        minX = var.getMin()
        maxX = var.getMax()
        nbinsX = var.getBinning().numBoundaries()-1
        print var, name, minX, maxX, nbinsX
        hsel = ROOT.TH2F("hsel_%s" % name,"hsel_%s" % name,nbinsX+3,minX-1.5*(maxX-minX)/nbinsX,maxX+1.5*(maxX-minX)/nbinsX,ncat+3,mincat-1.5,maxcat+1.5)
        for jcat,val in summary.iteritems():
            bd = float(val["selections"][isel])
            hsel.Fill(bd,float(jcat))
        csel = ROOT.TCanvas( "cat_opt_%s" % hsel.GetName(), "cat_opt_%s" % hsel.GetName() )
        csel.cd()
        hsel.Draw("box")

        csel_pj = ROOT.TCanvas( "cat_opt_%s_pj" % hsel.GetName(),  "cat_opt_%s_pj" %  hsel.GetName() )
        csel_pj.cd()
        hsel_pj = hsel.Clone()
        ### hsel_pj = hsel.ProjectionY()
        ### hsel_pj.Draw()
        hsel_pj.Draw("box")
        hsel_pj.SetFillColor(ROOT.kBlack)
        hsel_pj.SetLineColor(ROOT.kBlack)
        objs.append(hsel)
        objs.append(hsel_pj)
        objs.append(csel)
        objs.append(csel_pj)
        maxy = 0.
        pdfs = []
        for sig in signals:
            pdf = sig.getPdf(ndim+isel)
            pdf.SetLineColor(ROOT.kBlue)
            pdf.Scale(1./pdf.Integral())
            pdfs.append(pdf)
            maxy = max(maxy,pdf.GetMaximum())
            objs.append(pdf)
        for bkg in backgrounds:
            pdf = bkg.getPdf(ndim+isel)
            pdf.Scale(1./pdf.Integral())
            pdf.SetLineColor(ROOT.kRed)
            pdfs.append(pdf)
            maxy = max(maxy,pdf.GetMaximum())
            objs.append(pdf)
            
        for pdf in pdfs:
            pdf.Scale( (ncat+1.5)/maxy/pdf.Integral() )
            pdf.Draw("l same")
        hsel_pj.GetYaxis().SetNdivisions(500+ncat+3)
        csel_pj.SetGridy()
        hsel_pj.Draw("box same")


        for fmt in "png", "C":
            csel.SaveAs("%s.%s" % (csel.GetName(),fmt) )
            csel_pj.SaveAs("%s.%s" % (csel.GetName(),fmt) )

            
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

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
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
            make_option("-o", "--outdir",
                        action="store", type="string", dest="outdir",
                        default="./",
                        help="output file",
                        ),
            make_option("-l", "--load",
                        action="store", dest="settings", type="string",
                        default="",
                        help="json file containing settings"
                        ),
            make_option("--label",
                        action="store", dest="label", type="string",
                        default="",
                        help=""
                        ),
            make_option("-c", "--continue",
                        action="store", dest="cont", type="string",
                        default=False,
                        help="read-back the previous optimization step",
                        ),
            make_option("-r", "--reduce",
                        action="store_true", dest="reduce",
                        default=False,
                        help="collapse higher category orders to lower ones",
                        ),
            make_option("-s", "--settings",
                        action="store", dest="settings",
                        default=[],
                        help="append string to optimizer setting",
                        ),
            make_option("-n", "--ncat",
                        action="append", dest="range", type="int",
                        default=[],
                        help="collapse higher category orders to lower ones",
                        ),
            make_option("-d", "--dryrun",
                        action="store_true", dest="dry",
                        default=False,
                        help="do not run the optimization only find equidistant boundaries"
                        "\n  useful in conjunction with reduce and refit",
                        ),
            make_option("-R", "--refit",
                        action="store", dest="refit", type="int",
                        default=0,
                        help="refit the best point after reduction"
                        ),
            make_option("-S", "--splitreduce",
                        action="store_true", dest="splitreduce",
                        default=False,
                        help="search minimum by splitting and merging categories (equivalent to -r -R -d)"
                        ),
            make_option("--savecmpact",
                        action="store_true", dest="saveCompactTree",
                        default=False,
                        help="save the 2nd order models as TTrees"
                        ),
            make_option("-x", "--show-plots",
                        action="store_true", dest="showplots",
                        default=sys.flags.interactive,
                        help=""
                        ),
            make_option("--only-trees",
                        action="store_true", dest="onlytrees",
                        default=False,
                        help=""
                        ),
            ])

    (options, args) = parser.parse_args()
    loadSettings(options.settings, options)

    if options.infile == "":
        options.infile = "tmva%s.root" % options.label
        
    if not options.showplots:
        sys.argv.append("-b")

    if options.splitreduce:
        options.dry = True
        options.refit = 1
        options.reduce = True
        
    pprint(options.__dict__)

    import ROOT
    print ROOT.gROOT.IsBatch()
    

    ws=main(options,args)
