
#!/usr/bin/env python
# @(#)root/tmva $Id: MKmvaClassification.py,v 1.1.2.2 2012/01/24 16:05:17 mkenzie Exp $
# ------------------------------------------------------------------------------
# based on TMVA Python script: TMVAClassification.py
# ------------------------------------------------------------------------------

# --------------------------------------------
# Standard python import
import sys    # exit
import time   # time accounting
import getopt # command line parser
import os

# --------------------------------------------

# Default settings for command line arguments
DEFAULT_OUTFNAME = "TMVA"
DEFAULT_INFNAME  = "TMVA_input_training.root"
DEFAULT_TREESIG  = "sig"
DEFAULT_TREEBKG  = "bkg"
DEFAULT_METHODS  = "BDT"
DEFAULT_BACKGROUND = 1 # 1 uses the sidebands, 2 uses all
DEFAULT_MASS     = 123
DEFAULT_CAT      = -1
DEFAULT_PHIL     = "MIT"
DEFAULT_WIDTH    = 0.02
DEFAULT_TEST_TYPE = "grad"
DEFAULT_TEST_METHOD = "nTrees"

# Print usage help
def usage():
    print " "
    print "Usage: python %s [options]" % sys.argv[0]
    print "  -m | --methods    : gives methods to be run (default: '%s')" % DEFAULT_METHODS
    print "  -p | --philosophy : training philosophy (mit/ucsd) (default: '%s')" %DEFAULT_PHIL
    print "  -M | --mass       : gives Higgs Mass to train (default: '%i')" % DEFAULT_MASS  
    print "  -C | --cat        : Diphoton selection cat to run over (default: '%i')" % DEFAULT_CAT  
    print "  -B | --background : define which background samples to use (default: '%i' [sidebands])" % DEFAULT_BACKGROUND  
    print "  -i | --inputfile  : name of input ROOT file (default: '%s')" % DEFAULT_INFNAME
    print "  -o | --outputfile : name of output ROOT file containing results (default: '%s')" % DEFAULT_OUTFNAME
    print "  -t | --inputtrees : input ROOT Trees for signal and background (default: '%s %s')" % (DEFAULT_TREESIG, DEFAULT_TREEBKG)
    print "  -T | --test       : test different types (default: '%s %s')" % (DEFAULT_TEST_TYPE, DEFAULT_TEST_METHOD)
    print "  -v | --verbose"
    print "  -? | --usage      : print this help message"
    print "  -h | --help       : print this help message"
    print " "

# Check test input options
def checkTestType(testType, testMethod):
  if testType != "grad" or testType != "ada":
    print "ERROR: testType must be grad or ada"
    sys.exit(1)
  if testType=="ada":
    for result in ["nTrees","depth","nCuts","beta"]:
      if testMethod != result:
        print "ERROR: testMethod must be ", result
        sys.exit(1)
  if testType=="grad":
    for result in ["nTrees","depth","shrinkage","bagFrace","nCuts","nNM"]:
      if testMethod != result:
        print "ERROR: testMethod must be ", result
        sys.exit(1)

# Main routine
def main():

    try:
        # retrive command line options
        shortopts  = "m:p:M:C:B:i:t:T:o:vh?"
        opts, args = getopt.getopt( sys.argv[1:], shortopts )

    except getopt.GetoptError:
        # print help information and exit:
        print "ERROR: unknown options in argument %s" % sys.argv[1:]
        usage()
        sys.exit(1)

    infname     = DEFAULT_INFNAME
    methods     = DEFAULT_METHODS
    mass        = DEFAULT_MASS
    cat         = DEFAULT_CAT
    phil        = DEFAULT_PHIL
    outfname    = DEFAULT_OUTFNAME
    treeNameSig = DEFAULT_TREESIG
    treeNameBkg = DEFAULT_TREEBKG
    bkg_method  = DEFAULT_BACKGROUND
    width       = DEFAULT_WIDTH
    verbose     = False
    test        = False
    testType    = DEFAULT_TEST_TYPE
    methTest    = False
    testMethod  = DEFAULT_TEST_METHOD
    for o, a in opts:
        if o in ("-?", "-h", "--help", "--usage"):
            usage()
            sys.exit(0)
        elif o in ("-m", "--methods"):
            methods = a
        elif o in ("-M", "--mass"):
            mass = int(a)
        elif o in ("-C", "--cat"):
            cat = int(a)
        elif o in ("-p", "--philosophy"):
            phil = a
        elif o in ("-B", "--background"):
            bkg_method = int(a)
        elif o in ("-i", "--inputfile"):
            infname = a
        elif o in ("-o", "--outputfile"):
            outfname = a
        elif o in ("-T", "--test"):
            test=True
            temp = a.split('_')
            if len(temp)==1:
              testType = temp[0]
              if testType != "ada" or testType != "grad":
                print "ERROR: testType must be ada or grad not", testType
            elif len(temp)-temp.count('') == 2:
              methTest=True
              testType = temp[0]
              testMethod = temp[1]
              checkTestType(testType,testMethod)
            else:
              print "ERROR: need to give one or two test options"
              print temp
              sys.exit(1)
        elif o in ("-t", "--inputtrees"):
            a.strip()
            trees = a.rsplit( ' ' )
            trees.sort()
            trees.reverse()
            if len(trees)-trees.count('') != 2:
                print "ERROR: need to give two trees (each one for signal and background)"
                print trees
                sys.exit(1)
            treeNameSig = trees[0]
            treeNameBkg = trees[1]
        elif o in ("-v", "--verbose"):
            verbose = True

    if (width == 0.02) : width_str = "_2pt"
    elif (width == 0.07) : width_str = "_7pt"
    mass_str    = "_"+str("%3.1f" %mass)
    cat_str    = "_"+str(cat)
    if cat<0:
        cat_str    = "_all"
    if not os.path.isdir('TMVAOutput'):
      os.makedirs('TMVAOutput')
    if test:
      if methTest: 
        outfname = "TMVAOutput/"+outfname+"_"+phil+cat_str+"_test_"+testType+"_"+testMethod+".root"
      else:
        outfname = "TMVAOutput/"+outfname+"_"+phil+cat_str+"_test_"+testType+".root"
    else:
      outfname    = "TMVAOutput/"+outfname+"_"+phil+cat_str+".root"
      
    #treeNameSig = treeNameSig + mass_str 
    #treeNameBkg = treeNameBkg + mass_str  

    # Print methods
    mlist = methods.replace(' ',',').split(',')
    print "=== TMVAClassification: use method(s)..."
    for m in mlist:
        if m.strip() != '':
            print "=== - <%s>" % m.strip()

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
    gROOT.SetMacroPath( "/vols/cms03/mk1009/h2g/MVA/tmvaMacros/" )
    gROOT.Macro       ( "/vols/cms03/mk1009/h2g/MVA/tmvaMacros/TMVAlogon.C" )    
    gROOT.LoadMacro   ( "/vols/cms03/mk1009/h2g/MVA/tmvaMacros/TMVAGui.C" )
    
    # Import TMVA classes from ROOT
    from ROOT import TMVA

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

    # Get the signal and background trees for training
    signal_train      = input.Get( treeNameSig+"_train"+mass_str)
    signal_test      = input.Get( treeNameSig+"_test"+mass_str)

    background_train  = input.Get( treeNameBkg+"_train"+width_str+mass_str)
    background_test  = input.Get( treeNameBkg+"_test"+width_str+mass_str)

    # Global event weights (see below for setting event-wise weights)
    signalWeight     = 1.0
    backgroundWeight = 1.0

    # ====== register trees ====================================================
    factory.AddSignalTree    ( signal_train,    signalWeight     ,"train")
    factory.AddBackgroundTree( background_train,backgroundWeight ,"train")
    factory.AddSignalTree    ( signal_test,     signalWeight     ,"test")
    factory.AddBackgroundTree( background_test, backgroundWeight ,"test")
            
    # Set individual event weights (the variables must exist in the original
    # TTree)
    factory.SetBackgroundWeightExpression( "wt" )
    factory.SetSignalWeightExpression( "wt" )

    # Apply additional cuts on the signal and background sample. 
    # example for cut: mycut = TCut( "abs(var1)<0.5 && abs(var2-0.5)<1" )
    mycut = TCut( "fabs(deltaMOverM)<="+str(width)+" && bdtoutput > -0.5")#
    # Here, the relevant variables are copied over in new, slim trees that are
    # used for TMVA training and testing
    factory.PrepareTrainingAndTestTree( mycut, mycut, "nTrain_Signal=0:nTrain_Background=0:NormMode=NumEvents:!V")
    # Boosted Decision Trees
    # NEW PARAMETERS

    if (not test):
        # Likelihood
        factory.BookMethod( TMVA.Types.kLikelihood, "Likelihood"+phil, "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" )
        factory.BookMethod( TMVA.Types.kLikelihood, "LikelihoodD"+phil, "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" )
        #factory.BookMethod( TMVA.Types.kPDERS, "MultiLikelihood"+phil,"!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );
        
        # BDT
        factory.BookMethod( TMVA.Types.kBDT, "BDTada"+phil, "!H:!V:NTrees=200:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=1.0:SeparationType=GiniIndex:nCuts=50:PruneMethod=NoPruning" )
        factory.BookMethod( TMVA.Types.kBDT, "BDTgrad"+phil, "!H:!V:NTrees=200:MaxDepth=3:BoostType=Grad:Shrinkage=0.5:UseBaggedGrad:GradBaggingFraction=1.0:SeparationType=GiniIndex:nCuts=50:NNodesMax=10" )

    else:  #test
      # BDT ada
      if testType=="ada": 
        #if testMethod=="nTrees":
          for nTrees in [10,50,100,200,500]:
            for depth in [2,3]:
              factory.BookMethod( TMVA.Types.kBDT, "BDT_ada"+str(phil)+"_"+str(nTrees)+"t_"+str(depth)+"d","!H:!V:NTrees="+str(nTrees)+":nEventsMin=150:MaxDepth="+str(depth)+":BoostType=AdaBoost:AdaBoostBeta=1:SeparationType=GiniIndex:nCuts=50:PruneMethod=NoPruning")

       # if testMethod=="depth":
       #   for depth in [2,3]:
        #    factory.BookMethod( TMVA.Types.kBDT, "BDT_ada"+str(phil)+"_200t_"+str(depth)+"d_0.05b_50c","!H:!V:NTrees=200:nEventsMin=150:MaxDepth="+str(depth)+":BoostType=AdaBoost:AdaBoostBeta=0.05:SeparationType=GiniIndex:nCuts=50:PruneMethod=NoPruning")

        #if testMethod=="nCuts":
        #  for nCuts in [5,10,20,50,100,200]:
        #    factory.BookMethod( TMVA.Types.kBDT, "BDT_ada"+str(phil)+"_200t_50d_0.05b_"+str(nCuts)+"c","!H:!V:NTrees=200:nEventsMin=150:MaxDepth=50:BoostType=AdaBoost:AdaBoostBeta=0.05:SeparationType=GiniIndex:nCuts="+str(nCuts)+":PruneMethod=NoPruning")

        #if testMethod=="beta":
        #  for beta in [0.05,0.5,1.]:
        #    factory.BookMethod( TMVA.Types.kBDT, "BDT_ada"+str(phil)+"_200t_50d_"+str(beta)+"b_50c","!H:!V:NTrees=200:nEventsMin=150:MaxDepth=50:BoostType=AdaBoost:AdaBoostBeta="+str(beta)+":SeparationType=GiniIndex:nCuts=50:PruneMethod=NoPruning")
      
      # BDT grad
      if testType=="grad": 
        if testMethod=="nTrees":
          for nTrees in [10,50,100,200,500]:
            for depth in [2,3]:
              for shrinkage in [0.05,0.5,1.]:
                factory.BookMethod( TMVA.Types.kBDT, "BDT_grad"+str(phil)+"_"+str(nTrees)+"t_"+str(depth)+"d_"+str(shrinkage)+"s","!H:!V:NTrees="+str(nTrees)+":MaxDepth="+str(depth)+":BoostType=Grad:Shrinkage="+str(shrinkage)+":UseBaggedGrad:GradBaggingFraction=1:SeparationType=GiniIndex:nCuts=50:NNodesMax=10") 
        
        #if testMethod=="depth":
         # for depth in [2,3]:
          #  factory.BookMethod( TMVA.Types.kBDT, "BDT_ada"+str(phil)+"_200t_"+str(depth)+"d_0.05b_50c","!H:!V:NTrees=200:nEventsMin=150:MaxDepth="+str(depth)+":BoostType=AdaBoost:AdaBoostBeta=0.05:SeparationType=GiniIndex:nCuts=50:PruneMethod=NoPruning")


        #if testMethod=="shrinkage":
        #  for shrinkage in [0.05,0.1,0.5,1.]:
        #    factory.BookMethod( TMVA.Types.kBDT, "BDT_grad"+str(phil)+"_200t_"+str(shrinkage)+"s_1gb_50c_10nm","!H:!V:NTrees=200:BoostType=Grad:Shrinkage="+str(shrinkage)+":UseBaggedGrad:GradBaggingFraction=1:SeparationType=GiniIndex:nCuts=50:NNodesMax=10") 

        #if testMethod=="bagFrac":
        #  for bagFrac in [0.05,0.1,0.5,1.]:
         #   factory.BookMethod( TMVA.Types.kBDT, "BDT_grad"+str(phil)+"_200t_1s_"+str(bagFrac)+"gb_50c_10nm","!H:!V:NTrees=200:BoostType=Grad:Shrinkage=1:UseBaggedGrad:GradBaggingFraction="+str(bagFrac)+":SeparationType=GiniIndex:nCuts=50:NNodesMax=10") 

        #if testMethod=="nCuts":
         # for nCuts in [5,10,20,50,100,200]:
          #  factory.BookMethod( TMVA.Types.kBDT, "BDT_grad"+str(phil)+"_200t_1s_1gb_"+str(nCuts)+"c_10nm","!H:!V:NTrees=200:BoostType=Grad:Shrinkage=1:UseBaggedGrad:GradBaggingFraction=1:SeparationType=GiniIndex:nCuts="+str(nCuts)+":NNodesMax=10") 

        #if testMethod=="nNM":
         # for nNM in [10,100,500,1000,10000]:
          #  factory.BookMethod( TMVA.Types.kBDT, "BDT_grad"+str(phil)+"_200t_1s_1gb_50c_"+str(nNM)+"nm","!H:!V:NTrees=200:BoostType=Grad:Shrinkage=1:UseBaggedGrad:GradBaggingFraction=1:SeparationType=GiniIndex:nCuts=50:NNodesMax"+str(nNM)) 

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

# ----------------------------------------------------------

if __name__ == "__main__":
    main()

