#!/bin/env python

from optparse import OptionParser, make_option
import sys

def main(options,args):

    import ROOT
    import numpy
    
    ## Load library
    ROOT.gSystem.Load("../libLoopAll.so")
    from copy import copy

    ## input and output
    if ".root" in options.input:
        options.input = options.input.replace(".root","")
    elif not "/" in options.input:
        options.input = "%s/CMS-HGG" % options.input
        
    if options.output == "":
        options.output = copy(options.input)
    if ".root" in options.output:
        options.output = options.output.replace(".root","")

    ## input masses
    # Warning harded skipping of 145 for 7TeV
    inputMasses = numpy.array( [ float(m) for m in options.inputMasses.split(",") if not options.is7TeV or m!='145' ] )

    ## build list of mass points to interpolate
    interpolations = []

    if options.outputMasses != "":
        interpolations.append( (numpy.array( [float(m) for m in options.outputMasses.split(",")] ), "%s_interpolated.root" % options.output) )
    elif options.pval:
        interpolations.append( (numpy.arange(110,150.5,0.5), "%s_interpolated.root" % options.output) )
    else:
        if options.massScanType == "coarse" :
            interpolations.append( (numpy.concatenate( (numpy.arange(120,122.5,0.5), numpy.arange(123,126,0.1), numpy.arange(126,130.5,0.5)) ), "%s_mass_scan_coarse.root" % options.output) )
        elif options.massScanType == "fine" :
            interpolations.append( (numpy.arange(124.5,125.82,0.02), "%s_mass_scan_fine.root" % options.output) )
        else:
            interpolations.append( (numpy.concatenate( (numpy.arange(120,122.5,0.5), numpy.arange(123,124.5,0.1), numpy.arange(124.5,125.82,0.02), numpy.arange(125.9,126,0.1), numpy.arange(126,130.5,0.5)) ), "%s_mass_scan.root" % options.output) )

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)

    ## open input file
    fin  = ROOT.TFile.Open("%s.root" % options.input)
    fkeys = fin.GetListOfKeys()

    ## Read input histograms
    inputs = {}
    procs    = set()
    channels = set()
    tocopy = []
    for k in fkeys:
        if options.th1Prefix in k.GetName():
            decomp = str(k.GetName()).replace(options.th1Prefix,"").split("_",3)
            proc, mass, channel = decomp[0], float(decomp[2].replace("m","")), decomp[3]
            if mass in inputMasses:
                inputs[ (mass, proc, channel) ] = k.ReadObj()
                procs.add(proc)
                channels.add(channel)
            else:
                tocopy.append(k.ReadObj())
        elif "TH1" in k.GetClassName() or "TCanvas" in k.GetClassName():
            tocopy.append(k.ReadObj())
            
    ## Read input workspace
    ws   = fin.Get(options.workspace)

    ## Loop over interpolations
    doSmoothing = options.doSmoothing
    files = []
    for masses, file in interpolations:

        print "-------------------------------------------------------------------- "
        print "Will generate the following mass points: ", masses
        print "Output file   : ", file
        print "Input masses  : ", inputMasses
        print "Processes     : ", procs
        print "# of channels : ", len(channels)
        
        ## instantiate the mass interpolator
        interpolator = ROOT.MassInterpolator(len(masses), masses, doSmoothing, "th1f_sig_", "roohist_sig_" )

        ## loop over processes, channels and input masses and fill the interpolator
        print "\nFilling interpolator"
        print "-----------------------"
        for p in procs:
            print p
            for c in channels:
                for m in inputMasses:
                    interpolator.addInput(m, p, c, inputs[(m,p,c)] )
        
        ## prepare output
        print "\nPreparing output"
        print "-----------------------"
        fout = ROOT.TFile.Open(file,"recreate")
        files.append(fout)
        fout.cd()
        for o in tocopy:
            o.Write(o.GetName())
            
        ows  = ws.Clone()
        files.append(ows)
        mass = ows.var(options.massName)

        ## run interpolation
        print "Running interpolation"
        print "-----------------------"
        interpolator.outputFile(fout)
        interpolator.outputWs(ows,mass)
        interpolator.runInterpolation(options.is7TeV)

        ## save output
        print "\nSaving file"
        print "-----------------------"
        fout.cd()
        ows.Write( options.workspace )
        fout.Close()
        
        ## smooth the histograms only once
        doSmoothing = False
            
if __name__ == "__main__":
    parser = OptionParser(option_list=[
        make_option("-i", "--input",
                    action="store", type="string", dest="input",
                    default="",
                    help="", metavar=""
                    ),
        make_option("-o", "--output",
                    action="store", type="string", dest="output",
                    default="",
                    help="", metavar=""
                    ),
        make_option("--is7TeV",
                    action="store_true", dest="is7TeV",
                    default=False,
                    help="", metavar=""
                    ),
        make_option("-w","--workspace",
                    action="store", type="string", dest="workspace",
                    default="cms_hgg_workspace",
                    help="", metavar=""
                    ),
        make_option("-m","--massName",
                    action="store", type="string", dest="massName",
                    default="CMS_hgg_mass",
                    help="", metavar=""
                    ),
        make_option("-I","--inputMasses",
                    action="store", type="string", dest="inputMasses",
                    default="110,115,120,125,130,135,140,145,150",
                    help="", metavar=""
                    ),
        make_option("-O","--outputMasses",
                    action="store", type="string", dest="outputMasses",
                    default="",
                    help="", metavar=""
                    ),
        make_option("-p","--pval",
                    action="store_true", dest="pval",
                    default=True,
                    help="", metavar=""
                    ),
        make_option("-S","--massScan",
                    action="store_false", dest="pval",
                    help="", metavar=""
                    ),
        make_option("-s","--massScanType",
                    action="store", type="string", dest="massScanType",
                    default="",
                    help="", metavar=""
                    ),
        make_option("-t","--th1Prefix",
                    action="store", type="string", dest="th1Prefix",
                    default="th1f_sig_",
                    help="", metavar=""
                    ),
        make_option("--doSmoothing",
                    action="store_true", dest="doSmoothing",
                    default=True,
                    help="", metavar=""
                    ),
        make_option("--noSmoothing",
                    action="store_false", dest="doSmoothing",
                    help="", metavar=""
                    ),
        ]
                          )
    
    (options, args) = parser.parse_args()

    sys.argv.append("-b")
    main( options, args ) 
    
