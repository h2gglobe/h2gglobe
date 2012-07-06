
#!/usr/bin/env python
# @(#)root/tmva $Id: MKmvaClassification.py,v 1.2 2012/06/26 22:04:01 mkenzie Exp $
# ------------------------------------------------------------------------------
# based on TMVA Python script: TMVAClassification.py
# ------------------------------------------------------------------------------

# --------------------------------------------
# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import os

from optparse import OptionParser
# --------------------------------------------
parser = OptionParser()
parser.add_option("-i","--inputfile",dest="infile")
parser.add_option("-o","--outname",dest="outname",default="TMVA")
parser.add_option("-v","--verbose",dest="verbose",action="store_true",default=False)
parser.add_option("-p","--philosophy",dest="phil",default="MIT")
(options,args) = parser.parse_args()
  
infname = options.infile
outfname = options.outname
verbose = options.verbose
phil = options.phil
    
# Import ROOT classes
from ROOT import gSystem, gROOT, gApplication, TFile, TTree, TCut

# check ROOT version, give alarm if 5.18 
if gROOT.GetVersionCode() >= 332288 and gROOT.GetVersionCode() < 332544:
    print "*** You are running ROOT version 5.18, which has problems in PyROOT such that TMVA"
    print "*** does not run properly (function calls with enums in the argument are ignored)."
    print "*** Solution: either use CINT or a C++ compiled version (see TMVA/macros or TMVA/examples),"
    print "*** or use another ROOT version (e.g., ROOT 5.19)." 
    sys.exit(1)

# Logon not automatically loaded through PyROOT (logon loads TMVA library)
# load also GUI
#gROOT.SetMacroPath( "/vols/cms03/mk1009/h2g/TMVA/tmvaMacros/" )
#gROOT.Macro       ( "/vols/cms03/mk1009/h2g/TMVA/TMVAlogon.C" )    
#gROOT.LoadMacro   ( "/vols/cms03/mk1009/h2g/TMVA/tmvaMacros/TMVAGui.C" )


# Import TMVA classes from ROOT
from ROOT import TMVA

TMVA.Tools.Instance()

# Output file
outputFile = TFile( outfname, 'RECREATE' )

# Create instance of TMVA factory (see TMVA/macros/TMVAClassification.C for
# more factory options)
# All TMVA output can be suppressed by removing the "!" (not) in 
# front of the "Silent" argument in the option string
factory = TMVA.Factory( "TMVAClassification", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification")

# Set verbosity
factory.SetVerbose( verbose )

factory.AddVariable( "bdtoutput","BDT Output", 'F')
factory.AddVariable( "deltaMOverM","#DeltaM / M_{Hypth}.",  'F' )
input = TFile.Open( infname )

# Global event weights (see below for setting event-wise weights)
signalWeight     = 1.0
backgroundWeight = 1.0

# ====== register trees ====================================================
#factory.AddSignalTree    ( signal_train,    signalWeight     ,"train")
#factory.AddBackgroundTree( background_train,backgroundWeight ,"train")
#factory.AddSignalTree    ( signal_test,     signalWeight     ,"test")
#factory.AddBackgroundTree( background_test, backgroundWeight ,"test")
signalTrees=[]
signalTrees.append(input.Get("ggh_m124_pu2012"))
signalTrees.append(input.Get("vbf_m124_pu2012"))
signalTrees.append(input.Get("wzh_m124_pu2012"))
signalTrees.append(input.Get("tth_m124_pu2012"))
backgroundTrees=[]
backgroundTrees.append(input.Get("qcd_30_8TeV_ff"))
backgroundTrees.append(input.Get("qcd_30_8TeV_pf"))
backgroundTrees.append(input.Get("qcd_30_8TeV_pp"))
backgroundTrees.append(input.Get("qcd_40_8TeV_ff"))
backgroundTrees.append(input.Get("qcd_40_8TeV_pf"))
backgroundTrees.append(input.Get("qcd_40_8TeV_pp"))
backgroundTrees.append(input.Get("gjet_20_8TeV_ff"))
backgroundTrees.append(input.Get("gjet_20_8TeV_pf"))
backgroundTrees.append(input.Get("gjet_20_8TeV_pp"))
#backgroundTrees.append(input.Get("gjet_40_8TeV_ff"))
backgroundTrees.append(input.Get("gjet_40_8TeV_pf"))
backgroundTrees.append(input.Get("gjet_40_8TeV_pp"))
backgroundTrees.append(input.Get("diphojet_8TeV"))
backgroundTrees.append(input.Get("dipho_Box_10_8TeV"))
backgroundTrees.append(input.Get("dipho_Box_25_8TeV"))
backgroundTrees.append(input.Get("dipho_Box_250_8TeV"))

for tree in signalTrees:
  factory.AddSignalTree(tree, signalWeight)
for tree in backgroundTrees:
  factory.AddBackgroundTree(tree, backgroundWeight)

# Set individual event weights (the variables must exist in the original
# TTree)
factory.SetBackgroundWeightExpression( "weight" )
factory.SetSignalWeightExpression( "weight" )

# Apply additional cuts on the signal and background sample. 
# example for cut: mycut = TCut( "abs(var1)<0.5 && abs(var2-0.5)<1" )
mycut = TCut( "fabs(deltaMoverM)<=0.02 && bdtoutput >= -0.05")#
# Here, the relevant variables are copied over in new, slim trees that are
# used for TMVA training and testing
factory.PrepareTrainingAndTestTree( mycut, mycut, "nTrain_Signal=0:nTrain_Background=0:NormMode=NumEvents:!V")
# Boosted Decision Trees
# NEW PARAMETERS

# BDT
# 2011 default:
factory.BookMethod( TMVA.Types.kBDT, "BDTgrad"+phil, "!H:!V:NTrees=200:MaxDepth=3:BoostType=Grad:Shrinkage=0.5:UseBaggedGrad:GradBaggingFraction=1.0:SeparationType=GiniIndex:nCuts=50:NNodesMax=10" )
# ICHEP topup:
#factory.BookMethod( TMVA.Types.kBDT, "BDTgrad"+phil, "!H:!V:NTrees=500:MaxDepth=10:BoostType=Grad:Shrinkage=0.5:UseBaggedGrad:GradBaggingFraction=1.0:SeparationType=GiniIndex:nCuts=50:NNodesMax=10" )


# --------------------------------------------------------------------------------------------------
# ---- Now you can tell the factory to train, test, and evaluate the MVAs. 

# Train MVAs
#factory.OptimizeAllMethods()
factory.TrainAllMethods()
# Test MVAs
factory.TestAllMethods()

# Evaluate MVAs
factory.EvaluateAllMethods()    

# Save the output.
outputFile.Close()

print "=== wrote root file %s\n" % outfname
print "=== TMVAClassification is done!\n"

# open the GUI for the result macros    
#gROOT.ProcessLine( "TMVAGui(\"%s\")" % outfname )

# keep the ROOT thread running
#gApplication.Run() 



