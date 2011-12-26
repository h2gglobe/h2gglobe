This file is in h2gglobe/Marco/00README.txt

or:

http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/Marco/00README.txt?view=markup



Instructions 21/12/2011
-----------------------


#replace mpieri with your username

export CVS_RSH=ssh
export CVSROOT=:ext:mpieri@cmscvs.cern.ch:/cvs_server/repositories/CMSSW

#setenv CVS_RSH ssh
#setenv CVSROOT :ext:mpieri@cmscvs.cern.ch:/cvs_server/repositories/CMSSW

kinit mpieri@CERN.CH

### NEEDED at UCSD ~/SCHIF
source /home/users/mpieri/cmsset_default.sh
export SCRAM_ARCH="slc5_amd64_gcc434"
### end NEEDED at UCSD

scramv1 project CMSSW CMSSW_4_2_8

cd CMSSW_4_2_8/src

cmsenv

cvs co -r branch_for_42X -d HiggsAnalysis/HiggsTo2photons UserCode/HiggsAnalysis/HiggsTo2photons 
rm -r  HiggsAnalysis/HiggsTo2photons/h2gglobe

#check out the head of h2gglobe (Last hopefully working tag is marco_Dec26)

cvs co -d HiggsAnalysis/HiggsTo2photons/h2gglobe UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe
cd HiggsAnalysis/HiggsTo2photons/h2gglobe

### OBSOLETE
### You need some files that are in Marco and not in CVS:
####NOT ANYMORE
### cp ~mpieri/Aug2011/mpieri/CMSSW_4_2_8/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/Marco/dataevents*.txt Marco/.
### Except branchdef and Marco all the rest it should also work with:
### cvs update -r baseline_workspace_08Dec2011_nw

make clean; make -j 30



############################################
### Standard workspace: from Chris #########
############################################
### instructions on splitting and combining workspaces

### AT UCSD

# files needed to run at UCSD (in h2gglobe)
cp Marco/photonanalysis_ucsd.dat PhotonAnalysis_scripts/photonanalysis.dat
cp Marco/statanalysis_ucsd.dat PhotonAnalysis_scripts/statanalysis.dat
cp Marco/pu_weights_map.dat PhotonAnalysis_scripts/pu_weights_map.dat

cp Marco/subfit* PhotonAnalysis_scripts/.
cp Marco/datafiles_5fb_dec20_all_sm.dat PhotonAnalysis_scripts/.


cd PhotonAnalysis_scripts


### The following jobs create the RooWorkspace that is the input to the cL calculations

### python fitter.py -i datafiles_5fb_dec20_all_sm.dat -n numOfJobs -j jobNum
### jobNum goes from 0 to numOfJobs-1

### Example:

cd PhotonAnalysis_scripts
python fitter.py -i datafiles_5fb_dec20_all_sm.dat -n 100 -j 0

### END AT UCSD


### AT LXPLUS

cp Marco/subfit* PhotonAnalysis_scripts/.
cp Marco/statanalysis_ucsd.dat PhotonAnalysis_scripts/statanalysis.dat
cp Marco/pu_weights_map.dat PhotonAnalysis_scripts/.

cd PhotonAnalysis_scripts
python fitter.py -i datafiles_5fb.dat -n numOfJobs -j jobNum

#jobNum goes from 0 to numOfJobs-1

Example:
 
cd PhotonAnalysis_scripts
python fitter.py -i datafiles_5fb.dat -n 100 -j 0

### END AT LXPLUS



### start Chris

Instructions on merging the Workspaces:

At ucsd done with subfit* for different uafs

# subfit4 usage
uaf 4
cd yourh2gglobe
make -j 8
cd PhotonAnalysis_scripts
bash subfit4

# there is a subfit# for uaf 3,4,7,8,9

# usage for filestocombine.dat
# forloop for put input files
for i in {0..49}; do echo "Fil=CMS-HGG_4763_30-20_SM_2011_${i}.root"; done >> filestocombine.dat

#need to edit filestocombind.dat for outputname and placement of the Fil=*root lines

# you can also run combiner.py on a single file to change the fit
# you have to edit in PhotonAnalysis/src/StatAnalysis.cc to change the fit
python combiner.py

############################################
### End Standard workspace: from Chris #########
############################################


############################################
### Things for Jim and Matteo MVA #########
############################################

from h2gglobe
you need to copy the following files:

#FOR running at UCSD
cp Marco/photonanalysis_ucsd.dat PhotonAnalysis_scripts/photonanalysis.dat
#END FOR running at UCSD

cp Marco/Makefile ./Makefile
cp Marco/CommonParameters.h ./CommonParameters.h
cp Marco/statanalysis_ucsd.dat PhotonAnalysis_scripts/statanalysis.dat
cp Marco/plotvariables.dat PhotonAnalysis_scripts/plotvariables.dat
cp Marco/pu_weights_map.dat PhotonAnalysis_scripts/pu_weights_map.dat

******** Specific for JIM **************

cp Marco/LoopAll.cc LoopAll.cc
cp Marco/LoopAll.h LoopAll.h
cp Marco/minimal_statanalysis_input.dat PhotonAnalysis_scripts/minimal_statanalysis_input.dat
cp Marco/StatAnalysisExclusive_jim.h PhotonAnalysis/interface/StatAnalysisExclusive.h
cp Marco/StatAnalysisExclusive_jim.cc PhotonAnalysis/src/StatAnalysisExclusive.cc
cp Marco/cuts_marco.dat PhotonAnalysis_scripts/cuts.dat
cp Marco/counters.dat PhotonAnalysis_scripts/counters.dat

make clean; make -j 30

### to be deleted
### cp Marco/statanalysisexclusive.dat PhotonAnalysis_scripts/statanalysisexclusive.dat

### I THINK NOT ANYMORE 
### cp Marco/fitter.py PhotonAnalysis_scripts/fitter.py

### Before copying you should modify:
### ---------------------------------
### //NOW THEY ARE DEFAULT FOR MVA
### cuts_marco.dat (first few cuts)
### minimal_statanalysis_input.dat (put the branches you need)
### **** SPECIFIC FOR JIM ****
### In file minimal_statanalysis_input.dat remove the following couple of lines:
### ->
### FOR JIM TO UNCOMMENT
### and
### FOR JIM TO UNCOMMENT
### ->
### This will take as inputs few more variables needed for the MVA.


Then to run use from h2gglobe/PhotonAnalysis_scripts:

cd PhotonAnalysis_scripts/
rm  ../Marco/jimdatafiles_5fb_Dec14_90_190.dat.pevents
python fitter.py -i ../Marco/jimdatafiles_5fb_Dec14_90_190.dat --dryRun
# cd PhotonAnalysis_scripts/
python fitter.py -i ../Marco/jimdatafiles_5fb_Dec14_90_190.dat >& jimdatafiles_5fb_Dec14_90_190.log


### These jobs create 'optree' called Ntuple inside the output histogram file.
### From this one can do the mava analysis.
### The mva analysis workspace could be create directly in these jobes but
### StatAnalysis should be modifid for that 


### This is from matteo

Run the following commands on different UAF machines:
python fitter.py -i ../Marco/jimdatafiles_5fb_Dec14_90_190_dataA.dat >& jimdatafiles_5fb_Dec14_90_190_dataA.log
python fitter.py -i ../Marco/jimdatafiles_5fb_Dec14_90_190_dataB.dat >& jimdatafiles_5fb_Dec14_90_190_dataB.log
python fitter.py -i ../Marco/jimdatafiles_5fb_Dec14_90_190_signal.dat >& jimdatafiles_5fb_Dec14_90_190_signal.log
python fitter.py -i ../Marco/jimdatafiles_5fb_Dec14_90_190_bg1.dat >& jimdatafiles_5fb_Dec14_90_190_bg1.log
python fitter.py -i ../Marco/jimdatafiles_5fb_Dec14_90_190_bg2.dat >& jimdatafiles_5fb_Dec14_90_190_bg2.log

When the jobs are completed:

hadd final.root histograms_jimdatafiles_5fb_Dec14_90_190_dataA.root histograms_jimdatafiles_5fb_Dec14_90_190_dataB.root histograms_jimdatafiles_5fb_Dec14_90_190_signal.root histograms_jimdatafiles_5fb_Dec14_90_190_bg1.root histograms_jimdatafiles_5fb_Dec14_90_190_bg2.root

*********** END OF SPECIFIC FOR JIM ******************
############################################
### END ??? Things for Jim and Matteo MVA #########
############################################


then to run use from h2gglobe/PhotonAnalysis_scripts

jimdatafiles_5fb_Dec14_90_190.dat
and paste the lines in the first comment:
cd PhotonAnalysis_scripts/
rm  ../Marco/jimdatafiles_5fb_Dec14_90_190.dat.pevents
python fitter.py -i ../Marco/jimdatafiles_5fb_Dec14_90_190.dat --dryRun
cd PhotonAnalysis_scripts/
python fitter.py -i ../Marco/jimdatafiles_5fb_Dec14_90_190.dat >& jimdatafiles_5fb_Dec14_90_190.log


When you compile plotInteractive you need make clean; make -j 30
Instructions for plotting in Marco/printAll.txt


### END INSTRUCTIONS 20/12/2011
### END INSTRUCTIONS 20/12/2011
### END INSTRUCTIONS 20/12/2011






source /home/users/mpieri/cmsset_default.sh
export SCRAM_ARCH="slc5_amd64_gcc434"

cmsenv

Datafiles I am using:
marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_fast.dat

I have copied some MC but not all in my dir on nfs-6

You need LoopAll.h and .cc from Marco from the head to be copied in ../.

You need ~mpieri/Aug2011/mpieri/CMSSW_4_2_8/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/Marco/dataevents.txt to be kept in Marco

The order is jul aug oct 2011B hoping that it always takes the same sequence of files.




http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/Marco/

~capalmer/CMS/Higgs2gg/Limits/nov14_forAN/428p7/src/HiggsAnalysis/HiggsTo2photons/h2gglobe_new/PhotonAnalysis_scripts/histograms_CMS-HGG_4686_nov17_tryExcl.root
~capalmer/CMS/Higgs2gg/Limits/nov14_forAN/428p7/src/HiggsAnalysis/HiggsTo2photons/h2gglobe_new2/PhotonAnalysis_scripts/

diff -r ~capalmer/CMS/Higgs2gg/Limits/nov14_forAN/428p7/src/HiggsAnalysis/HiggsTo2photons/h2gglobe_new2/ . |grep "diff -r" |grep -v CVS

ls -l  ~capalmer/CMS/Higgs2gg/Limits/nov14_forAN/428p7/src/HiggsAnalysis/HiggsTo2photons/h2gglobe_new3/LoopAll.cc

# FOR PLOTINTERACTIVE:

root -l
.L Marco/tdrstyle.C
setTDRStyle();
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  //gSystem->Load("VertexAnalysis/lib/libh2gglobeVertexAnalysis.so");
gSystem->Load("libLoopAll.so");
LoopAll* m=new LoopAll();
m->myPlotInteractiveSetup("PhotonAnalysis_scripts/histograms_marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_4DEC_fastnew_ptom.root","Hgg")
m->myPlotInteractive("PhotonAnalysis_scripts/histograms_marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_4DEC_fastnew_ptom.root");

root -l
.L Marco/tdrstyle.C
setTDRStyle();
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  //gSystem->Load("VertexAnalysis/lib/libh2gglobeVertexAnalysis.so");
gSystem->Load("libLoopAll.so");
LoopAll* m=new LoopAll();
m->myPlotInteractiveSetup("PhotonAnalysis_scripts/histograms_marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_1DEC.root","Hgg")
m->myPlotInteractive("PhotonAnalysis_scripts/histograms_marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_1DEC.root");

root -l
.L Marco/tdrstyle.C
setTDRStyle();
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  //gSystem->Load("VertexAnalysis/lib/libh2gglobeVertexAnalysis.so");
gSystem->Load("libLoopAll.so");
LoopAll* m=new LoopAll();
m->myPlotInteractiveSetup("PhotonAnalysis_scripts/histograms_marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_1DEC_fastnew.root","Hgg")
m->myPlotInteractive("PhotonAnalysis_scripts/histograms_marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_1DEC_fastnew.root");

root -l
.L Marco/tdrstyle.C
setTDRStyle();
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  //gSystem->Load("VertexAnalysis/lib/libh2gglobeVertexAnalysis.so");
gSystem->Load("libLoopAll.so");
LoopAll* m=new LoopAll();
m->myPlotInteractiveSetup("PhotonAnalysis_scripts/histograms_marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_1DEC.root","Hgg")
m->myPlotInteractive("PhotonAnalysis_scripts/histograms_marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_1DEC.root");

m->myPlotInteractiveSetup("PhotonAnalysis_scripts/histograms_CMS-HGG_4686pb_jettagging_ff_only120_29nov_withSM.root","Hgg")
m->myPlotInteractive("PhotonAnalysis_scripts/histograms_CMS-HGG_4686pb_jettagging_ff_only120_29nov_withSM.root");

m->myPlotInteractiveSetup("PhotonAnalysis_scripts/histograms_CMS-HGG_4686pb_jettagging_ff_only120_26nov_debug.root","Hgg")
m->myPlotInteractive("PhotonAnalysis_scripts/histograms_CMS-HGG_4686pb_jettagging_ff_only120_26nov_debug.root");


_29nov_withSM

m->myPlotInteractiveSetup("PhotonAnalysis_scripts/histograms_CMS-HGG_4686pb_jettagging_ff_only120_26nov_ptom_90_190.root","Hgg")
m->myPlotInteractive("PhotonAnalysis_scripts/histograms_CMS-HGG_4686pb_jettagging_ff_only120_26nov_ptom_90_190.root");

m->myPlotInteractiveSetup("PhotonAnalysis_scripts/histograms_CMS-HGG_4686pb_jettagging_ff_only120_26nov_ptom_100_160.root","Hgg")
m->myPlotInteractive("PhotonAnalysis_scripts/histograms_CMS-HGG_4686pb_jettagging_ff_only120_26nov_ptom_100_160.root");

m->myPlotInteractiveSetup("PhotonAnalysis_scripts/histograms_CMS-HGG_4686pb_jettagging_ff_only120_26nov_90_190.root","Hgg")
m->myPlotInteractive("PhotonAnalysis_scripts/histograms_CMS-HGG_4686pb_jettagging_ff_only120_26nov_90_190.root");



Instructions 16/11/2011
-----------------------

cd CMSSW_4_2_8/src
cvs co -r regression_Sept30 HiggsAnalysis/HiggsToGammaGamma
cvs co -r branch_for_42X -d HiggsAnalysis/HiggsTo2photons UserCode/HiggsAnalysis/HiggsTo2photons 

rm -r  HiggsAnalysis/HiggsTo2photons/h2gglobe

#check out the head of h2gglobe (we should make a tag

cvs co -d HiggsAnalysis/HiggsTo2photons/h2gglobe UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe

cp Marco/CommonParameters.h ./CommonParameters.h
#cp Marco/HistoContainer.cc ./HistoContainer.cc
#cp Marco/HistoContainer.h ./HistoContainer.h
#cp Marco/CounterContainer.cc ./CounterContainer.cc
#cp Marco/CounterContainer.h ./CounterContainer.h
cp Marco/LoopAll.h ./LoopAll.h 
cp Marco/LoopAll.cc ./LoopAll.cc 

cp Marco/minimal_statanalysis_input.dat PhotonAnalysis_scripts/minimal_statanalysis_input.dat
#cp Marco/cuts.dat PhotonAnalysis_scripts/cuts.dat
#cp Marco/cuts_marco.dat PhotonAnalysis_scripts/cuts.dat
cp Marco/plotvariables.dat PhotonAnalysis_scripts/plotvariables.dat
#cp Marco/inputfiles.dat PhotonAnalysis_scripts/inputfiles.dat
# cp Marco/looper.py PhotonAnalysis_scripts/looper.py
# cp Marco/looper_input.dat PhotonAnalysis_scripts/looper_input.dat
cp Marco/photonanalysis.dat PhotonAnalysis_scripts/photonanalysis.dat
cp Marco/reduction_output.dat PhotonAnalysis_scripts/reduction_output.dat
######cp Marco/configProducer.py PhotonAnalysis_scripts/python/configProducer.py
# cp Marco/PhotonAnalysis.cc PhotonAnalysis/src/PhotonAnalysis.cc
#cp Marco/StatAnalysisExclusive.h PhotonAnalysis/interface/StatAnalysisExclusive.h
#cp Marco/StatAnalysisExclusive.cc PhotonAnalysis/src/StatAnalysisExclusive.cc
##cp Marco/statanalysis.dat PhotonAnalysis_scripts/statanalysis.dat
cp Marco/statanalysisexclusive.dat PhotonAnalysis_scripts/statanalysisexclusive.dat


???
cp Marco/fitter.py PhotonAnalysis_scripts/fitter.py

cp Marco/Makefile ./Makefile

???
#cp Marco/pu_weights_map.dat PhotonAnalysis_scripts/pu_weights_map.dat

cp Marco/datafiles_5fb.dat PhotonAnalysis_scripts/datafiles_5fb.dat
cp Marco/datafiles_5fb_LL_33_23.dat PhotonAnalysis_scripts/datafiles_5fb_LL_33_23.dat
cp Marco/datafiles_5fb_LL_33_23_all.dat PhotonAnalysis_scripts/datafiles_5fb_LL_33_23_all.dat
cp Marco/datafiles_5fb_LL_33_23_all_ff.dat PhotonAnalysis_scripts/datafiles_5fb_LL_33_23_all_ff.dat


#FOR MARCO ADD THIS
cp Marco/statanalysis.dat PhotonAnalysis_scripts/statanalysis.dat
cp Marco/photonanalysis.dat PhotonAnalysis_scripts/photonanalysis.dat
cp  Marco/plotvariables.dat PhotonAnalysis_scripts/plotvariables.dat
cp Marco/cuts_marco.dat PhotonAnalysis_scripts/cuts.dat
#cp Marco/cuts_marco_100160.dat PhotonAnalysis_scripts/cuts.dat
cp Marco/StatAnalysisExclusive_marco_newnewnew.h PhotonAnalysis/interface/StatAnalysisExclusive.h
cp Marco/StatAnalysisExclusive_marco_newnewnew.cc PhotonAnalysis/src/StatAnalysisExclusive.cc
###cp Marco/marcodatafiles_5fb_LL_33_23_all_ff_120.dat PhotonAnalysis_scripts/datafiles_5fb_LL_33_23_all_ff_120.dat
###cp Marco/marcodatafiles_5fb_LL_33_23_all_ff.dat PhotonAnalysis_scripts/datafiles_5fb_LL_33_23_all_ff.dat

#cp Marco/datafiles_5fb.dat PhotonAnalysis_scripts/datafiles_5fb.dat
#cp Marco/datafiles_5fb_LL_33_23.dat PhotonAnalysis_scripts/datafiles_5fb_LL_33_23.dat
#cp Marco/datafiles_5fb_LL_33_23_all.dat PhotonAnalysis_scripts/datafiles_5fb_LL_33_23_all.dat
#cp Marco/datafiles_5fb_LL_33_23_all_ff.dat PhotonAnalysis_scripts/datafiles_5fb_LL_33_23_all_ff.dat


diff Marco/CommonParameters.h ./CommonParameters.h

diff Marco/HistoContainer.cc ./HistoContainer.cc
diff Marco/HistoContainer.h ./HistoContainer.h

diff Marco/CounterContainer.cc ./CounterContainer.cc
diff Marco/CounterContainer.h ./CounterContainer.h

diff Marco/minimal_statanalysis_input.dat PhotonAnalysis_scripts/minimal_statanalysis_input.dat
diff Marco/cuts.dat PhotonAnalysis_scripts/cuts.dat
diff Marco/cuts_marco.dat PhotonAnalysis_scripts/cuts.dat
diff Marco/plotvariables.dat PhotonAnalysis_scripts/plotvariables.dat
diff Marco/inputfiles.dat PhotonAnalysis_scripts/inputfiles.dat
# diff Marco/looper.py PhotonAnalysis_scripts/looper.py
# diff Marco/looper_input.dat PhotonAnalysis_scripts/looper_input.dat
diff Marco/photonanalysis.dat PhotonAnalysis_scripts/photonanalysis.dat
diff Marco/reduction_output.dat PhotonAnalysis_scripts/reduction_output.dat

diff Marco/configProducer.py PhotonAnalysis_scripts/python/configProducer.py
# diff Marco/PhotonAnalysis.cc PhotonAnalysis/src/PhotonAnalysis.cc

diff Marco/datafiles_5fb.dat PhotonAnalysis_scripts/datafiles_5fb.dat
diff Marco/statanalysis.dat PhotonAnalysis_scripts/statanalysis.dat

diff Marco/StatAnalysisExclusive.h PhotonAnalysis/interface/StatAnalysisExclusive.h
diff Marco/StatAnalysisExclusive.cc PhotonAnalysis/src/StatAnalysisExclusive.cc
diff Marco/StatAnalysisExclusive_marco.cc PhotonAnalysis/src/StatAnalysisExclusive.cc
diff Marco/statanalysisexclusive.dat PhotonAnalysis_scripts/statanalysisexclusive.dat
diff Marco/fitter.py PhotonAnalysis_scripts/fitter.py

diff Marco/StatAnalysisExclusive.h PhotonAnalysis/interface/StatAnalysis.h
diff Marco/StatAnalysisExclusive.cc PhotonAnalysis/src/StatAnalysis.cc

diff Marco/Makefile ./Makefile


If needed modify in PhotonAnalysis/interface/PhotonAnalysis.h the line similar into:
#include "../../../../HiggsToGammaGamma/interface/GBRForest.h"

Then make

make clean; make -j 40

cd PhotonAnalysis_scripts/.
rm datafiles_5fb.dat.pevents
python fitter.py -i datafiles_5fb.dat --dryRun
python fitter.py -i datafiles_5fb.dat 

cd PhotonAnalysis_scripts/.
python fitter.py -i datafiles_5fb_LL_33_23.dat --dryRun
python fitter.py -i datafiles_5fb_LL_33_23.dat 

cd PhotonAnalysis_scripts/.
rm datafiles_5fb_LL_33_23_all_short.dat.pevents
python fitter.py -i datafiles_5fb_LL_33_23_all_short.dat --dryRun
python fitter.py -i datafiles_5fb_LL_33_23_all_short.dat 

cd PhotonAnalysis_scripts/.
rm datafiles_5fb_LL_33_23_all.dat.pevents
python fitter.py -i datafiles_5fb_LL_33_23_all.dat --dryRun
python fitter.py -i datafiles_5fb_LL_33_23_all.dat 





cd PhotonAnalysis_scripts/.
rm datafiles_5fb_LL_33_23_all_ff.dat.pevents
python fitter.py -i datafiles_5fb_LL_33_23_all_ff.dat --dryRun
python fitter.py -i datafiles_5fb_LL_33_23_all_ff.dat >&  datafiles_5fb_LL_33_23_all_ff_final.log

========
cp Marco/StatAnalysisExclusive_marco_newnewnew.h PhotonAnalysis/interface/StatAnalysisExclusive.h
cp Marco/StatAnalysisExclusive_marco_newnewnew.cc PhotonAnalysis/src/StatAnalysisExclusive.cc

tkdiff Marco/StatAnalysisExclusive_marco_newnewnew.h PhotonAnalysis/interface/StatAnalysis.h
tkdiff Marco/StatAnalysisExclusive_marco_newnewnew.cc PhotonAnalysis/src/StatAnalysis.cc

d PhotonAnalysis_scripts/.
cp ../Marco/cuts_marco_debug.dat cuts.dat
cp ../Marco/cuts_marco.dat cuts.dat
cp ../Marco/cuts_marco_ptom.dat cuts.dat
cp ../Marco/cuts_marco_100160.dat cuts.dat
cp ../Marco/cuts_marco_100160_ptom.dat cuts.dat
cp ../Marco/marcodatafiles_5fb_LL_33_23_all_ff_120.dat datafiles_5fb_LL_33_23_all_ff_120.dat

cp ../Marco/cuts_marco.dat cuts.dat
cp ../Marco/marcodatafiles_5fb_LL_33_23_all_ff_120_2011A.dat datafiles_5fb_LL_33_23_all_ff_120_2011A.dat
cp ../Marco/marcodatafiles_5fb_LL_33_23_all_ff_120_2011B.dat datafiles_5fb_LL_33_23_all_ff_120_2011B.dat

cd PhotonAnalysis_scripts/.
rm ../Marco/marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_fast.dat.pevents
python fitter.py -i ../Marco/marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_fast.dat --dryRun
python fitter.py -i ../Marco/marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_fast.dat >& marcodatafiles_5fb_LL_33_23_all_ff_120_withSM_1DecM_fast.log


#########
#########
rm datafiles_5fb_LL_33_23_all_ff_120.dat.pevents
python fitter.py -i datafiles_5fb_LL_33_23_all_ff_120.dat --dryRun
python fitter.py -i datafiles_5fb_LL_33_23_all_ff_120.dat >&  datafiles_5fb_LL_33_23_all_ff_120_26nov_90_190.log
python fitter.py -i datafiles_5fb_LL_33_23_all_ff_120.dat >&  datafiles_5fb_LL_33_23_all_ff_120_26nov_ptom_90_190.log
python fitter.py -i datafiles_5fb_LL_33_23_all_ff_120.dat >&  datafiles_5fb_LL_33_23_all_ff_120_26nov_100_160.log
python fitter.py -i datafiles_5fb_LL_33_23_all_ff_120.dat >&  datafiles_5fb_LL_33_23_all_ff_120_26nov_ptom_100_160.log

python fitter.py -i datafiles_5fb_LL_33_23_all_ff_120_2011A.dat >&  datafiles_5fb_LL_33_23_all_ff_120_26nov_90_190_2011A.log
python fitter.py -i datafiles_5fb_LL_33_23_all_ff_120_2011B.dat >&  datafiles_5fb_LL_33_23_all_ff_120_26nov_90_190_2011B.log

cat datafiles_5fb_LL_33_23_all_ff_120_23nov.log

#########
#########

When finished if you want co commit everything:

rm CommonParameters.h 
rm HistoContainer.cc HistoContainer.h PhotonAnalysis_scripts/cuts.dat PhotonAnalysis_scripts/plotvariables.dat PhotonAnalysis_scripts/inputfiles.dat PhotonAnalysis_scripts/python/configProducer.py  PhotonAnalysis/src/PhotonAnalysis.cc PhotonAnalysis_scripts/looper.py PhotonAnalysis_scripts/looper_input.dat PhotonAnalysis_scripts/photonanalysis.dat PhotonAnalysis_scripts/reduction_output.dat PhotonAnalysis_scripts/datafiles_5fb.dat PhotonAnalysis_scripts/statanalysis.dat
rm PhotonAnalysis/interface/StatAnalysisExclusive.h PhotonAnalysis/src/StatAnalysisExclusive.cc Makefile PhotonAnalysis_scripts/statanalysisexclusive.dat PhotonAnalysis_scripts/fitter.py
rm CounterContainer.cc CounterContainer.h PhotonAnalysis_scripts/minimal_statanalysis_input.dat
rm PhotonAnalysis_scripts/pu_weights_map.dat

cvs update -A

=====================================
=====================================
=====================================
=====================================
END FOR NOW
=====================================
=====================================
=====================================
=====================================




UCSD:
-----

cd ~/Aug2011/mpieri

source /home/users/mpieri/cmsset_default.sh
export SCRAM_ARCH="slc5_amd64_gcc434"

cmsenv
export CVS_RSH=ssh
export CVSROOT=:ext:mpieri@cmscvs.cern.ch:/cvs_server/repositories/CMSSW

#setenv CVS_RSH ssh
#setenv CVSROOT :ext:mpieri@cmscvs.cern.ch:/cvs_server/repositories/CMSSW

kinit mpieri@CERN.CH

scramv1 project CMSSW CMSSW_4_2_8

cd CMSSW_4_2_8/src
cvs co -r regression_Sept30 HiggsAnalysis/HiggsToGammaGamma
cvs co -r branch_for_42X -d HiggsAnalysis/HiggsTo2photons UserCode/HiggsAnalysis/HiggsTo2photons 

rm -r  HiggsAnalysis/HiggsTo2photons/h2gglobe

#cvs co -r nw_16_11_11 -d HiggsAnalysis/HiggsTo2photons/h2gglobe UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe

cvs co -d HiggsAnalysis/HiggsTo2photons/h2gglobe UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe

cmsenv
scramv1 b

cd HiggsAnalysis/HiggsTo2photons/h2gglobe
make clean; make



Hi All,

I made a tag for the h2gglobe package I used for producing the Workspace from which we obtained the limit plot yesterday (nw_16_11_11 ).
Below is a recipe for making the workspace at LXPLUS using the 4.69/fb and the mix of Summer11 and Fall11 MC since some people have asked
how to produce it.

cmsrel CMSSW_4_2_8
cd CMSSW_4_2_8/src

cvs co -r nw_16_11_11 -dHiggsAnalysis/HiggsTo2photons/h2gglobe UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe
cvs co HiggsAnalysis/HiggsToGammaGamma (HEAD is ok, but maybe Chris Palmer can comment on what he used for Ntuple production)

scramv1 b -j8
cd HiggsAnalysis/HiggsTo2photons/h2gglobe

cmsenv
make clean ; make -j8

cd PhotonAnalysis_scripts
python fitter.py -i datafiles_5fb.dat --dryRun

the last command will configure to make sure all the files are there for reweighting etc and the ntuples themselves exist.
Note this also create a file called datafiles_5fb.dat.pevents which you should not alter. If you change the file datafiles_5fb.dat, then you should
delete the datafiles_5fb.dat.pevents and rerun the above command.

next, you can actually run the jobs by the following:

python fitter.py -i datafiles_5fb.dat -n 30 -j 0

which runs the first of 30 jobs (its up to you how many you submit but i find 30 completes in a reasonable time)
You can submit each of the jobs individually to the batch machine at LXPLUS if you know how. I put a little script here
~nckw/pubic/makesub.py
which will create submission scripts for each job, just do

python makesub.py 30 $PWD datafiles_5fb.dat

and then you can submit each one to batch with the usual
qsub -q 1nh sub0.sh  etc...
1nh should be long enough but check your logs in case not.

You will now have the 30 subjobs (CMS-HGG_4686pb_(0-29).root ) which you should move into a folder called cms-hgg-workspaces,
then finally

python combiner.py -i filestocombine.dat

will merge those files into one called CMS-HGG_4686pb.root

This is the workspace which goes to the binned limit setting for Generated MC points.
There are Macros in h2gglobe/Macros which do the interpolation and diagnostics but I haven't the time to explain them all here. Please contact me if you want to take the binned limit setting through. Otherwise, from that workspace, Josh has a script to overlay a signal model to do the unbinned one.

Let me know if anything is unclear or you run into any problems

Cheers,
Nick


  hfile->cd();

//int dummy = fscanf(file,"%d plot=%d idummy=%d",&Nvar, &typplotall, &idummystart);
//h2d[i], typplot[i], histoncat[i], nbinsx[i], nbinsy[i], lowlim[i],highlim[i],lowlim2[i],highlim2[i],varnamescread[i]
  plotvartree = new TTree("plotvariables","globe plotvariables provenance information");
  plotvartree->Branch("Nvar", &Nvar, "Nvar/I");
  plotvartree->Branch("typplotall", &typplotall, "typplotall/I");
  plotvartree->Branch("doplot", &plothistoplotitPI, "plothistoplotitPI[Nvar]/I");
  plotvartree->Branch("h2d", &h2d, "h2d[Nvar]/I");
  plotvartree->Branch("typplot", &typplot, "typplot[Nvar]/I");
  plotvartree->Branch("histoncat", &histoncat, "histoncat[Nvar]/I");
  plotvartree->Branch("histoncatindtonames", &histoncatindtonames, "histoncatindtonames[Nvar]/I");
  plotvartree->Branch("nbinsx", &nbinsx, "nbinsx[Nvar]/I");
  plotvartree->Branch("nbinsy", &nbinsy, "nbinsy[Nvar]/I");
  plotvartree->Branch("lowlim", &lowlim, "lowlim[Nvar]/F");
  plotvartree->Branch("highlim", &highlim, "highlim[Nvar]/F");
  plotvartree->Branch("lowlim2", &lowlim2, "lowlim2[Nvar]/F");
  plotvartree->Branch("highlim2", &highlim2, "highlim2[Nvar]/F");

  tca_xaxislabels = new TClonesArray("TObjString",Nvar);
  plotvartree->Branch("xaxislabels", "TClonesArray", &tca_xaxislabels, 32000, 0);
  for(int iplot=0;iplot!=Nvar;++iplot) { 
    new ((*tca_xaxislabels)[iplot]) TObjString(); 
    ((TObjString *)tca_xaxislabels->At(iplot))->SetString(xaxislabel[iplot]);
  }

  tca_yaxislabels = new TClonesArray("TObjString",Nvar);
  plotvartree->Branch("yaxislabels", "TClonesArray", &tca_yaxislabels, 32000, 0);
  for(int iplot=0;iplot!=Nvar;++iplot) { 
    new ((*tca_yaxislabels)[iplot]) TObjString(); 
    ((TObjString *)tca_yaxislabels->At(iplot))->SetString(yaxislabel[iplot]);
  }

  tca_plotvarnames = new TClonesArray("TObjString",Nvar);
  plotvartree->Branch("plotvarnames", "TClonesArray", &tca_plotvarnames, 32000, 0);
  for(int iplot=0;iplot!=Nvar;++iplot) { 
    new ((*tca_plotvarnames)[iplot]) TObjString(); 
    ((TObjString *)tca_plotvarnames->At(iplot))->SetString(varnamescread[iplot]);
  }

  plotvartree->Branch("Nvarcats", &Nvarcats, "Nvarcats/I");
  plotvartree->Branch("catid", &catid, "catid[Nvarcats]/I");
  plotvartree->Branch("ncats", &ncats, "ncats[Nvarcats]/I");
  tca_plotvarcatnames = new TClonesArray("TObjString",Ncatvar);
  plotvartree->Branch("plotvarcatnames", "TClonesArray", &tca_plotvarcatnames, 32000, 0);
  int catvartemp=0;
  for(int i=0; i<Nvarcats; i++) {
    for(int j=0; j<ncats[i]; j++) {
      new ((*tca_plotvarcatnames)[catvartemp]) TObjString(); 
      ((TObjString *)tca_plotvarcatnames->At(catvartemp))->SetString(catnames[i][j]);
      catvartemp++;
    } 
  } 
  std::cout << "Ncatvar: " << Ncatvar << std::endl;
  std::cout << "catvartemp: " << catvartemp << std::endl;

  plotvartree->Fill();
  plotvartree->Write(0,TObject::kWriteDelete);

//int dummy = fscanf(file,"%d intL=%f rtree=%d %s %s\n",&nfiles, &intlumi, &makeOutputTree, outFilNam, histFilNam);
//int dummy = fscanf(file,"typ=%d ind=%d draw=%d Nam=%s Fil=%s tot=%d red=%d lum=%f xsec=%f kfac=%f scal=%f\n", &itype[i], &histoindfromfiles[i], &histoplotit[i], filesshortnam[i], files[i],&ntot[i],&nred[i], &lumi[i], &xsec[i], &kfactor[i], &scale[i]);
  inputfiletree = new TTree("inputfiles","globe inputfiles provenance information");
  inputfiletree->Branch("nfiles", &mp->nfiles, "nfiles/I");
  inputfiletree->Branch("nindfiles", &mp->nindfiles, "nindfiles/I");
  inputfiletree->Branch("intlumi", &mp->intlumi, "intlumi/F");
  inputfiletree->Branch("makeOutputTree", &mp->makeOutputTree, "makeOutputTree/I");
  tca_histfilename = new TClonesArray("TObjString",1);
  inputfiletree->Branch("histfilename", "TClonesArray", &tca_histfilename, 32000, 0);
  new ((*tca_histfilename)[0]) TObjString(); 
  ((TObjString *)tca_histfilename->At(0))->SetString(mp->histFilNam);
  inputfiletree->Branch("itype", &mp->itype, "itype[nfiles]/I");
  inputfiletree->Branch("histoind", &mp->histoindfromfiles, "histoindfromfiles[nfiles]/I");
  inputfiletree->Branch("infoind", &mp->infoind, "infoind[nindfiles]/I");
  inputfiletree->Branch("histoplotit", &mp->histoplotit, "histoplotit[nfiles]/I");
  inputfiletree->Branch("ntot", &mp->ntot, "ntot[nfiles]/I");
  inputfiletree->Branch("nred", &mp->nred, "nred[nfiles]/I");
  inputfiletree->Branch("lumi", &mp->lumi, "lumi[nfiles]/F");
  inputfiletree->Branch("xsec", &mp->xsec, "xsec[nfiles]/F");
  inputfiletree->Branch("kfactor", &mp->kfactor, "kfactor[nfiles]/F");
  inputfiletree->Branch("scale", &mp->scale, "scale[nfiles]/F");
  tca_inshortnames = new TClonesArray("TObjString",1);
  tca_infilenames = new TClonesArray("TObjString",mp->nfiles);
  inputfiletree->Branch("inshortnames", "TClonesArray", &tca_inshortnames, 32000, 0);
  inputfiletree->Branch("infilenames", "TClonesArray", &tca_infilenames, 32000, 0);
  for(int ifile=0;ifile!=mp->nfiles;++ifile) { 
    new ((*tca_inshortnames)[ifile]) TObjString(); 
    new ((*tca_infilenames)[ifile]) TObjString(); 
    ((TObjString *)tca_inshortnames->At(ifile))->SetString(mp->filesshortnam[ifile]);
    ((TObjString *)tca_infilenames->At(ifile))->SetString(mp->files[ifile]);
  }
  inputfiletree->Fill();
  inputfiletree->Write(0,TObject::kWriteDelete);




In mpUtil:

void mpUtil::PlotInteractive(TString tag, TString inputfile) { 

      loops->myPlotInteractiveSetup(this,inputfile,tag); 
      if(tag == "Hgg") {
        //loops->myPlotInteractiveHggSetup(this,inputfile);
      } else if(tag == "Hww") {
	//loops->myPlotInteractiveHwwSetup(this,inputfile); 
      } else if(tag == "Hzz") {
      } else if(tag == "Mwl") {
      } else if(tag == "Elizabeth") {
	//loops->myPlotInteractiveElizabethSetup(this,inputfile); 
      }
      loops->myPlotInteractive(this,inputfile); 
}

Script:




root -l
//.L matteo/tdrstyle.C
//setTDRStyle();
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  //gSystem->Load("VertexAnalysis/lib/libh2gglobeVertexAnalysis.so");
gSystem->Load("libLoopAll.so");
LoopAll* m=new LoopAll();
m->myPlotInteractiveSetup("Marco/histograms_CMS-HGG_4686_nov17_tryExcl.root","Hgg");
m->myPlotInteractive("Marco/histograms_CMS-HGG_4686_nov17_tryExcl.root");



 Marco/histograms_CMS-HGG_4686_nov17_tryExcl.root



mpUtil* m=new mpUtil();
m->PlotInteractive("Hgg","oct13_ucsdcuts_smnewnew2.root");


varname[0]: all_mass
varname[1]: pt
varname[2]: eta
varname[3]: mass
varname[4]: helicityAngle
varname[5]: decayAngle
varname[6]: pho_pt
varname[7]: pho1_pt
varname[8]: pho2_pt
varname[9]: pho_eta
varname[10]: pho1_eta
varname[11]: pho2_eta
varname[12]: pho_r9
varname[13]: pho_n
varname[14]: mass_pf
varname[15]: mass_ff
varname[16]: pho1_pt_presel
varname[17]: pho2_pt_presel
varname[18]: pho1_pt_sel
varname[19]: pho2_pt_sel
varname[20]: pho1_eta_presel
varname[21]: pho2_eta_presel
varname[22]: pho1_eta_sel
varname[23]: pho2_eta_sel
varname[24]: cut_Mgg_nminus1
varname[25]: cut_Mgg_sequential
varname[26]: cut_VBFLeadJPt_nminus1
varname[27]: cut_VBFLeadJPt_sequential
varname[28]: cut_VBFSubJPt_nminus1
varname[29]: cut_VBFSubJPt_sequential
varname[30]: cut_VBF_Mjj_nminus1
varname[31]: cut_VBF_Mjj_sequential
varname[32]: cut_VBF_dEta_nminus1
varname[33]: cut_VBF_dEta_sequential
varname[34]: cut_VBF_Zep_nminus1
varname[35]: cut_VBF_Zep_sequential
varname[36]: cut_VBF_dPhi_nminus1
varname[37]: cut_VBF_dPhi_sequential
varname[38]: cut_VBF_Mgg_nminus1
varname[39]: cut_VBF_Mgg_sequential
varname[40]: cut_VBF_Mgg4_nminus1
varname[41]: cut_VBF_Mgg4_sequential


###cp ../LoopAll.h ../LoopAll_marco.h 
###cp ../LoopAll.cc ../LoopAll_marco.cc 

cp Marco/StatAnalysisExclusive_marco_newnewnew.h PhotonAnalysis/interface/StatAnalysisExclusive.h
cp Marco/StatAnalysisExclusive_marco_newnewnew.cc PhotonAnalysis/src/StatAnalysisExclusive.cc


JIM

cp Marco/LoopAll.cc LoopAll.cc
cp Marco/LoopAll.h LoopAll.h
cp Marco/minimal_statanalysis_input.dat PhotonAnalysis_scripts/minimal_statanalysis_input.dat
cp Marco/StatAnalysisExclusive_jim.h PhotonAnalysis/interface/StatAnalysisExclusive.h
cp Marco/StatAnalysisExclusive_jim.cc PhotonAnalysis/src/StatAnalysisExclusive.cc
cp Marco/cuts_marco.dat PhotonAnalysis_scripts/cuts.dat

END JIM


cp Marco/CommonParameters.h ./CommonParameters.h
cp Marco/statanalysis.dat PhotonAnalysis_scripts/statanalysis.dat
cp Marco/photonanalysis.dat PhotonAnalysis_scripts/photonanalysis.dat
???? NOT FOR ME cp Marco/reduction_output.dat PhotonAnalysis_scripts/reduction_output.dat

cp Marco/PhotonAnalysis.h PhotonAnalysis/interface/PhotonAnalysis.h


??? YES FOR NOW 
cp Marco/fitter.py PhotonAnalysis_scripts/fitter.py


cp Marco/Makefile ./Makefile


cp Marco/LoopAll.h LoopAll.h 
cp Marco/LoopAll.cc LoopAll.cc 

cp Marco/plotvariables.dat PhotonAnalysis_scripts/plotvariables.dat
cp Marco/cuts_marco.dat PhotonAnalysis_scripts/cuts.dat
# NO cp Marco/minimal_statanalysis_input_marco.dat PhotonAnalysis_scripts/minimal_statanalysis_input.dat
cp Marco/StatAnalysisExclusive_marco_newnewnew.h PhotonAnalysis/interface/StatAnalysisExclusive.h
cp Marco/StatAnalysisExclusive_marco_newnewnew.cc PhotonAnalysis/src/StatAnalysisExclusive.cc


cvs update -r baseline_workspace_08Dec2011_nw



cp Marco/marcodatafiles_5fb_LL_33_23_all_ff_120.dat PhotonAnalysis_scripts/datafiles_5fb_LL_33_23_all



cvs diff  -r baseline_workspace_08Dec2011_nw

cvs update -r baseline_workspace_08Dec2011_nw PhotonAnalysis/interface/PhotonAnalysis.h
cvs update -r baseline_workspace_08Dec2011_nw PhotonAnalysis/interface/StatAnalysis.h
cvs update -r baseline_workspace_08Dec2011_nw PhotonAnalysis/src/PhotonAnalysis.cc

cvs update -r baseline_workspace_08Dec2011_nw PhotonAnalysis/src/StatAnalysis.cc

cvs update -r baseline_workspace_08Dec2011_nw PhotonAnalysis_scripts/datafiles_5fb.dat
cvs update -r baseline_workspace_08Dec2011_nw PhotonAnalysis_scripts/energy_scale_offsets.dat
cvs update -r baseline_workspace_08Dec2011_nw PhotonAnalysis_scripts/pu_weights_map.dat
cvs update -r baseline_workspace_08Dec2011_nw PhotonAnalysis_scripts/smearing_sigma_and_errors.dat
cvs update -r baseline_workspace_08Dec2011_nw Reduction/data.txt
cvs update -r baseline_workspace_08Dec2011_nw Reduction/mk_reduction_dat.sh

TKDIFF TKDIFF TKDIFF TKDIFF

tkdiff Marco/CommonParameters.h ./CommonParameters.h
tkdiff Marco/statanalysis.dat PhotonAnalysis_scripts/statanalysis.dat
tkdiff Marco/photonanalysis.dat PhotonAnalysis_scripts/photonanalysis.dat
???? NOT FOR ME tkdiff Marco/reduction_output.dat PhotonAnalysis_scripts/reduction_output.dat

tkdiff Marco/PhotonAnalysis.h PhotonAnalysis/interface/PhotonAnalysis.h





??? YES FOR NOW 
tkdiff Marco/fitter.py PhotonAnalysis_scripts/fitter.py


tkdiff Marco/Makefile ./Makefile


tkdiff Marco/LoopAll.h LoopAll.h 
tkdiff Marco/LoopAll.cc LoopAll.cc 

tkdiff Marco/plotvariables.dat PhotonAnalysis_scripts/plotvariables.dat
tkdiff Marco/cuts_marco.dat PhotonAnalysis_scripts/cuts.dat
tkdiff Marco/minimal_statanalysis_input_marco.dat PhotonAnalysis_scripts/minimal_statanalysis_input.dat
tkdiff Marco/StatAnalysisExclusive_marco_newnewnew.h PhotonAnalysis/interface/StatAnalysis.h
tkdiff Marco/StatAnalysisExclusive_marco_newnewnew.cc PhotonAnalysis/src/StatAnalysis.cc










-- run 167675 Event 911781274
phoEt[0] = 95.2403
phoEt[1] = 27.4618
regressed phoEt[0] = 95.8048
regressed phoEt[1] = 27.9181
i Passing preselection on photon 0
j Passing preselection on photon 1

cicPhotonID pho 0
phocat = 0
Id cut 0: combIso = -0.242797 cut > 3.8
Id cut 1: combIso(badVertex) = -0.352329 cut > 11.7
Id cut 2: trackIso = 0 cut > 3.5
Id cut 3: Sigma_IetaIeta = 0.00909123 cut > 0.0106
Id cut 4: Sigma_IetaIeta = 0 cut > 0.082
Id cut 5: R9 = 0.956989 cut < 0.94
Id cut 6: R9 = 99 cut < 1

cicPhotonID pho 1
phocat = 2
Id cut 0: combIso = 1.77178 cut > 1.77


-- run: 176547, event: 187480804
4 photon candidates:
pho 0 passed all selections.
pho 1 failed on combined Iso, cicPhoton ID combIso: 3.67753 (< 2.2 failed)
pho 2 failed on combined Iso, cicPhoton ID combIso: 2.38617 (< 2.2 failed)
pho 3 failed on preselection, SCEta= 2.57595 (|SCEta|<2.5 failed)


-- run: 177730, event: 58708626
2 photon candidates:
Leading photon failed on bad combined Iso cut.
cicPhoton ID combIso: 1.43261 (< 2.2)
cicPhoton ID badCombIso: 3.52066 (< 3.4 failed)


We don't have below 2 events in our ntuples. Further investigation ongoing.
Run 176547 Event 48254336
Run 177053 Event 563451379




# Final instructions


take branchdef and Marco from the head

all the rest:
cvs update -r baseline_workspace_08Dec2011_nw


cp Marco/Makefile ./Makefile
cp Marco/CommonParameters.h ./CommonParameters.h
cp Marco/statanalysis.dat PhotonAnalysis_scripts/statanalysis.dat
cp Marco/photonanalysis.dat PhotonAnalysis_scripts/photonanalysis.dat
#### cp Marco/PhotonAnalysis.h PhotonAnalysis/interface/PhotonAnalysis.h

??? YES FOR NOW 
cp Marco/fitter.py PhotonAnalysis_scripts/fitter.py

cp Marco/plotvariables.dat PhotonAnalysis_scripts/plotvariables.dat
cp Marco/pu_weights_map.dat PhotonAnalysis_scripts/pu_weights_map.dat

JIM

cp Marco/LoopAll.cc LoopAll.cc
cp Marco/LoopAll.h LoopAll.h
cp Marco/minimal_statanalysis_input.dat PhotonAnalysis_scripts/minimal_statanalysis_input.dat
cp Marco/StatAnalysisExclusive_jim.h PhotonAnalysis/interface/StatAnalysisExclusive.h
cp Marco/StatAnalysisExclusive_jim.cc PhotonAnalysis/src/StatAnalysisExclusive.cc
cp Marco/cuts_marco.dat PhotonAnalysis_scripts/cuts.dat




END JIM


????
cp Marco/minimal_statanalysis_input.dat PhotonAnalysis_scripts/minimal_statanalysis_input.dat
cp Marco/StatAnalysisExclusive.h PhotonAnalysis/interface/StatAnalysisExclusive.h
cp Marco/StatAnalysisExclusive.cc PhotonAnalysis/src/StatAnalysisExclusive.cc
cp Marco/cuts_marco.dat PhotonAnalysis_scripts/cuts.dat
end ????




# If you don't take the head of branchdef

cp Matteo/addBranch .
./addBranch pho_pfiso_myneutral03 F pho_n MAX_PHOTONS
./addBranch pho_pfiso_myneutral04 F pho_n MAX_PHOTONS
./addBranch pho_pfiso_myphoton03 F pho_n MAX_PHOTONS
./addBranch pho_pfiso_myphoton04 F pho_n MAX_PHOTONS
./addBranch pho_pfiso_mycharged04 "std::vector<std::vector<float> >"
./addBranch pho_pfiso_mycharged03 "std::vector<std::vector<float> >"
