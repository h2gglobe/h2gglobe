Zee validation:
============

Histograms are filled in MassFactorizedMvaAnalysis::fillZeeControlPlots

Configuration file (run with fitter.py, I usually use 50 jobs):

AnalysisScripts/massfac_mva_binned/datafiles_zeevalidation.dat

also available broken down into run periods (with PU weighting and energy smearing as appropriate):

AnalysisScripts/massfac_mva_binned/datafiles_zeevalidation_runAB.dat
AnalysisScripts/massfac_mva_binned/datafiles_zeevalidation_runC.dat
AnalysisScripts/massfac_mva_binned/datafiles_zeevalidation_runD.dat

Running any of the above needs change in common/minimal_analysis_input.dat:  remove ":2" from the following lines so that trigger selection is applied to MC as well as data:

hlt_bit:2
hlt_n:2
hlt_path_names_HLT:2

Also comment out line in TriggerSelection.cc to avoid printout (Zee validation does not use run dependent paths so that the trigger is the same for data and MC (implementation in PhotonAnalysis::Init).  This means that for some runs it can be that one of the paths is not in the menu, hence there is an error message, though for all runs there is at least one valid path):
  std::cerr << "TriggerSelection did not find path " << *it << std::endl;

Can optionally set fillEscaleTrees=1 for filling trees (not needed for standard validation plots which are filled directly with the above configuration).  To fill the same trees for signal MC, use:

AnalysisScripts/massfac_mva_binned/datafiles_massfacmva_Moriond2013_sigEscaleTrees.dat


Zee validation plotting scripts:
--------------------------------------

Plotting scripts are located here:

http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/David/ZeeValidation/

Running the configuration described above using the following CVS tag for h2gglobe: zeeValidation_moriond2013
produces the following histogram file:

/eos/cms/store/group/phys_higgs/cmshgg/histograms_moriond_david/histograms_CMS-HGG_zeevalidation.root

The full set of Zee validation plots can be produced by copying the above file to h2gglobe/David/ZeeValidation and running in that directory:

./runZee.sh

and for even more plots:

./runEtaR9bins.sh  (idmva, diphomva and sigmaE in eta,R9 bins)
./runCorrelations.sh (2D correlation plots)

To do things properly, one should interpolate the systematic bands to prevent them vanishing where the up and down templates cross over (in practice the effect is small, but this should be done for public plots).  Andre David kindly wrote a script to do this:

Macros/MorphBands/morphbands.py

It can be run like this from the Macros/MorphBands directory:

python morphbands.py david_template.py ../../David/ZeeValidation/histograms_CMS-HGG_zeevalidation.root 

This puts a file in David/ZeeValidation/ which if renamed to histograms_CMS-HGG_zeevalidation_envelopes.root can be picked up by the following scripts, which can be used in place of the aforementioned runZee.sh and runEtaR9bins.sh:

./runZee_interp.sh
./runEtaR9bins_interp.sh

The public plots without ratio plots:
https://twiki.cern.ch/twiki/pub/CMSPublic/Hig13001TWiki/diphomva.png
https://twiki.cern.ch/twiki/pub/CMSPublic/Hig13001TWiki/idmva_EB_nvtx.png
https://twiki.cern.ch/twiki/pub/CMSPublic/Hig13001TWiki/idmva_EE_nvtx.png

can be reproduced running:

root -b -q diphomva_public.C
root -b -q diphomva_nvtx.C

These scripts also require the interpolated file produced as described above.



data/MC comparison plots:
=====================

public plot https://twiki.cern.ch/twiki/pub/CMSPublic/Hig13001TWiki/mass.png and AN2013_008_v4 figs 87, 88, 143, 145 and various other plots can be produced using the following macro:

http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/David/massPlot.C

To reproduce public plot ( https://twiki.cern.ch/twiki/pub/CMSPublic/Hig13001TWiki/mass.png ):
root -b -l -q 'massPlot.C("0","all_mass","massfit",0,0,1,1)'

To produce full set of data/MC comparison plots (including variables other than mass), you can run:

http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/David/runMassPlot.sh

NB. can change "massfit" to "cic" in this script to make same plots for cic analysis.  Other options are indicated at the top of massPlot.C

The script uses these files as input:
root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/histograms_moriond_david/histograms_CMS-HGG_massfit.root
root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/histograms_moriond_david/histograms_CMS-HGG_cic.root

The above were produced using the following config files:
AnalysisScripts/massfac_mva_binned/datafiles_massfacmva_Moriond2013_dataMC.dat
AnalysisScripts/baseline/datafiles_cubased_Moriond2013_dataMC.dat

(Note for D.F.:  original file locations:
/afs/cern.ch/work/f/futyand/moriond_approval/CMSSW_5_3_7/src/h2gglobe/AnalysisScripts/massfit/combined/histograms_CMS-HGG.root
/afs/cern.ch/work/f/futyand/moriond_approval/CMSSW_5_3_7/src/h2gglobe/AnalysisScripts/cic/combined/histograms_CMS-HGG.root)



Single photon efficiency for cut based photon ID:
======================================

AN2013_008_v4 figs 38 and 39 can be produced by running:

http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe/David/eff_pho_public.C

uses as input: root://eoscms//eos/cms/store/group/phys_higgs/cmshgg/histograms_moriond_david/histograms_CMS-HGG_effPho.root
produced by running AnalysisScripts/baseline/datafiles_cutbasedvbftag_Moriond2013.dat for mH=125GeV signal MC only, setting:
plotvariables common/plotvariables_effPho.dat

Histogram filling is done in: StatAnalysis::fillSignalEfficiencyPlots
Need to uncomment the following line in StatAnalysis.cc:  fillSignalEfficiencyPlots(weight, l);

Note:  efficiency denominator is:  pT1/mgg>0.33, pT2/mgg>0.25 and sceta1, sceta2 in fiducial region

(Note for D.F.:  original file location:
/afs/cern.ch/work/f/futyand/moriond_approval/CMSSW_5_3_7/src/h2gglobe/AnalysisScripts/cic_sig/combined/histograms_CMS-HGG.root)
