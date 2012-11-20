#!/bin/env python

from optparse import OptionParser, make_option
import sys

def main(options,args):

    import ROOT
    import numpy
    import json
    import os
    
    mydir = os.path.dirname(sys.argv[0])
    ROOT.gSystem.SetIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include -I%s" % (mydir))
    ROOT.gROOT.LoadMacro("%s/SplitMiniTree.C+" % mydir )

    partitions = json.load(open(options.partitions))

    if options.globe:
        pdfName = "pdf_data_pol_model_8TeV_cat%d"
        wsName  = "cms_hgg_workspace"
        dsName  = "roohist_data_mass_cat%d" 
    else:
        pdfName = "CMS_hgg_mvacat%d_8TeV_bkgshape"
        wsName  = "wbkg"
        dsName  = "databinned_mvacat%d_8TeV" 

    output = options.input
    if options.outdir != "":
        output = "%s/%s" % ( options.outdir, os.path.basename(output) )
    splitMiniTree = ROOT.SplitMiniTree(options.input,options.wspace,wsName,pdfName,dsName,output)
    
    for ipart,partition in enumerate(partitions):
        for run,lumi,event in partition:
            splitMiniTree.addToPartition(ipart,run,lumi,event)

    splitMiniTree.splitCopy()

    print "npart: %d" % len(partitions)

if __name__ == "__main__":
    parser = OptionParser(option_list=[
        make_option("-i", "--input",
                    action="store", type="string", dest="input",
                    default="",
                    help="default: %default", metavar=""
                    ),
        make_option("-w", "--wspace",
                    action="store", type="string", dest="wspace",
                    default="",
                    help="default: %default", metavar=""
                    ),
        make_option("-d", "--outdir",
                    action="store", type="string", dest="outdir",
                    default="",
                    help="default: %default", metavar=""
                    ),
        make_option("-g", "--globe",
                    action="store_true", dest="globe",
                    default="",
                    help="default: %default", metavar=""
                    ),
        make_option("-m", "--MIT",
                    action="store_false", dest="globe",
                    default="",
                    help="", metavar=""
                    ),
        make_option("-p", "--partitions",
                    action="store", type="string", dest="partitions",
                    default="",
                    help="default : %default", metavar=""
                    ),
        
        ])
    
    (options, args) = parser.parse_args()

    sys.argv.append("-b")
    main( options, args ) 
    
