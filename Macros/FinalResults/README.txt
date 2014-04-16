# Instructions for running the Globe Hgg Legacy results
# the CombinedLimit package should already have been checkout along with
# h2gglobe
cd h2gglobe/Macros/FinalResults
cmsenv

# First the jobs are created and submitted using combineHarvester.py The
# easiest way to use this is to create a configuration file for the jobs
# (.dat) see the example in combineHarvesterOptions.dat. Use --skipWorkspace
# (or add skipWorkspace to the config line)
# if already ran once before or already with --dryRun

./combineHarvester.py -d combineHarvesterOptions.dat -q 8nh (--skipWorkspace)

# there are many options for producing differnt types of scans including
# all combine options can be passed by inclding
# the following in the config line...

	--opts=--combineopt1<argument>+--combineopt2<argument>+....

# To run Post-fit expected, add the keyword postFit to both run the fit (for the datajob) 
# and send/load it for making the post-fit expected (keyword expected) results. Note --skipWorkspace will also
# skip the fit so don't use it on the first go.
# The order is therefore important since we need to advertise the fit to data
# BEFORE creating the post-fit expected jobs. The following is a good example
# from the mass scan...

        outDir=MultiPdf/mvaComb/MHScan			method=MHScan	mh=125 mhLow=122. mhHigh=128. jobs=25 pointsperjob=8  postFit # postFit here means the fit will be perfomed to data and the snapshot saved in the new workspace
        outDir=MultiPdf/mvaComb/MHScanExpected		method=MHScan	mh=125 mhLow=122. mhHigh=128. jobs=25 pointsperjob=8  skipWorkspace postFit expected=1 expectSignal=1 opts=--setPhysicsModelParameters<RV=1.0,RF=1.0> # since this is "expected" the snapshot will be loaded from the previous line and only the POIs here are set to there pre-fit expected values. 

# You should notice that post-fit is not exactly the same as adding the
# combine arg "--toysFrequentist" since that will find the conditional
# minimum, not the Global one in the data.   

# Once the jobs are finished, the jobs are hadded with 

./combineHarvester.py --hadd Folder
# Plots can be made with the makeCombinePlots.py. The plots will be in the
# style of the Legacy Paper and there are many options to check (run with
# --help) 
#Below is an  example to produce the legacy plot of the MH scan

./makeCombinePlots -f ./MultiPdf_exact/mvaComb/MHScan/MHScan.root -c 1 -s 1 -w 3 -n "Obs Stat+Syst" -t "#sqrt{s}=7TeV L=5.1fb^{-1}, #sqrt{s}=8TeV L=19.7fb^{-1}" --mh -o results_paper/mva_mh_scan_envelope_noexp_exact --xaxis 123.8,125.8 --MHtext "0.36:0.6:Floating #mu_{VBF,VH}, #mu_{ggH,ttH}"


# THE INSTRUCTIONS BELOW HERE ARE OUT OF DATE ###########
# Instructions for running the final results in globe for either the binned or parametric model
# You should check out the head of globe in a CMSSW area and also get V02-06-00 of 
# HiggsAnalysis/CombinedLimit and compile it

# e.g for binned limits
cd h2gglobe/Macros/FinalResults
mkdir binned_model
cd binned_model
ln -s ../runFinalResults.py
ln -s ../subParametricToBatch.py
ln -s ../subBinnedToBatch.py

# copy 2011 and 2012 workspaces into this directory
# copy 2011 and 2012 datacards into this directory

# first we want to submit the combine jobs to the batch:
./runFinalResults.py --card2011="name_of_card_2011.txt" --card2012="name_of_card_2012.txt" --unblind

# once all the jobs have finished we can plot them
./runFinalResults.py --unblind --makePlots --lumi2011=5.3 --lumi2012=19.6

# all the relevant plots should appear in a directory called plots

# e.g for the parametric limits
cd h2gglobe/Macros/FinalResults
mkdir parametric_model
cd parametric_model
ln -s ../runFinalResults.py
ln -s ../subParametricToBatch.py
ln -s ../subBinnedToBatch.py

# copy 2011 and 2012 signal and data (with background pols) workspaces into this directory
# copy 2011 and 2012 datacards into this directory

# first we want to submit the combine jobs to the batch:
./runFinalResults.py --card2011="name_of_card_2011.txt" --card2012="name_of_card_2012.txt" --unblind -p

# once all the jobs have finished we can plot them
./runFinalResults.py --unblind -p --makePlots --lumi2011=5.3 --lumi2012=19.6

# all the relevant plots should appear in a directory called plots

