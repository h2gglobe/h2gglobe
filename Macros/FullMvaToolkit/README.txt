******************************
**** FULMVATOOLKIT advice ****
******************************
Author:   Matthew Kenzie
Email:    matthew.william.kenzie@cern.ch
Modified: 08.08.23

- More detailed instructions can be found in fullmvatoolkit_25_05_12.pdf
- Also see:
  https://twiki.cern.ch/twiki/bin/viewauth/CMS/Higgs2GAnalyzer#FullMvaToolkit

- You should run inside a CMSSW_6_1_X release

Assume we have the histogram file which comes with creating the baseline H->gg
Result from Globe called "histograms_CMS-HGG.root"
This result uses the 2D category mapping rather than TMVA training, to use
TMVA 2D simply remove "--use2DcatMap" everywhere and put bin edges in .dat
after training at Rebin step
Add --unblind once we are ready 

# Train from trees!
./bin/SidebandMVATraining -i histograms_CMS-HGG.root -m 2DOpt

# Make the histos
./bin/RunFullMvaAnalysis -i histograms_CMS-HGG.root -o fmt_hists.root -N -T -U dat/LegacyConfig.dat -v --use2DcatMap  

# Run fit and Rebin
./bin/RunFullMvaAnalysis -i fmt_hists.root -o fmt_workspace.root  -U dat/LegacyConfig.dat -v  --use2DcatMap  -P 

# Backup!
cp fmt_workspace.root fmt_workspace_safe.root

# Interpolate and make background model
./bin/RunFullMvaAnalysis -i fmt_hists.root -o fmt_workspace.root  -U dat/LegacyConfig.dat -v  --use2DcatMap -N -I -b -P 

# At Theory migration systematics 
python python/makeMigrationsSystematics.py -i histograms_CMS-HGG_7(8)TeV_systematics.root -o fmt_workspace.root -m SidebandTrainingOut_7(8)TeV.root -s 7(8)

# Datacards
./bin/RunFullMvaAnalysis -i fmt_hists.root -o fmt_workspace.root  -U dat/LegacyConfig.dat -v  --use2DcatMap -N -d  -P 

 
