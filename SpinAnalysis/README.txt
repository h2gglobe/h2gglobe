Author: M Kenzie

Runs a very simple spin analysis using flat trees produced in globe
Run globe with ../AnalysisScripts/massfac_mva_binned/datafiles_massfacmva_spinStudies.dat
There is a line commented at the bottom which lets you flick between the massfactorized and cutbased
When the jobs finish hadd the histogram files with ../AnalysisScripts hist_combiner.py

You can then compile this package with cd SpinAnalysis; cmsenv; make

Run:

./bin/diyBuildWorkspace histograms_CMS-HGG.root (cut-based) 
./bin/diyBuildWorkspace histograms_CMS-HGG.root --isMassFac (mass-fac)

This will create a workspace called HggSpin_*.root

You can then run ./subDIY.py to submit toys to the batch to calculate the separation
You can check the progress of your jobs with ./checkDIY.py

When the jobs are complete and you have hadded the output you can calculate the separation and plot the result with:
./bin/diyPlot

