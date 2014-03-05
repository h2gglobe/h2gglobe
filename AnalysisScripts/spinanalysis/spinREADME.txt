###############################################
#### INSTRUCITON FOR RUNNING SPIN ANALYSIS ####
###############################################

1) Modify AnalysisScripts/spinanalysis/spinanalysis.dat to run with desired options

2) Run mk_fitter.py on AnalysisScripts/spinanalysis/datafiles_spinanalysis.dat:

  cd h2gglobe/AnalysisScripts
  ./mk_fitter.py -i spinanalysis/datafiles_spinanalysis.dat -n 100 -l name_of_output_dir

3) Submit jobs:
  
  ./submit_fitter.py -d name_of_output_dir -q 8nh

3) Combine jobs and histogram files as necessary:

  ./check_fitter.py name_of_output_dir
  ./hist_combiner.py name_of_output_dir

4) Run massInterpolator and renomailzation script. This is required to renormalize the graviton samples
   and can be done in a number of ways. The standard is such that the graviton samples are normalised to 
   give the same expected significance as the SM over the inclusive categories. This gives a normalisation
   factor of 1./0.94 = 1.064 for gg and 1./0.84=1.19 (You can work out what this value should be by running getGravPrefitNorm.py
   on a workspace which has only a single cosTheta category - i.e is inclusive in cosTheta)

   ./massInterpolator -i CMS-HGG.root -I 125,126 -O 125,125.1,125.2,125.3,125.4,125.5,125.6,125.7,125.8,125.9,126 -H gg_grav:1.064:ggh,vbf,wzh,tth -H qq_grav:1.19:ggh,vbf,wzh,tth

5) Make the datacards. This is done with spinanalysis/mk_spin_card.py. To get the default collection 
   of cards run mk_cards.sh passing the datafile and the signal file names:

  ./mk_cards.sh CMS-HGG_data.root CMS-HGG_sig_interpolated.root

6) There are now 5 datacards - one SM only (to throw asimov and fit for channel comp), one GRAV_GG 
   only (to throw asimov for channel comp), one GRAV_QQ only (to throw asimov for channel compt),
   the nominal one containing SM and gg_grav (for the separation test stat toys), one for fqqbar scan

7) Go to a clean folder somewhere which has HiggsAnalysis/CombinedLimit compiled and executable.
   Copy the datacards and the two workspace files to this directory. Get (or link) the following scripts
   for AnalysisScripts/spinanalysis:
        - sub_spin_jobs.py
        - merge_jobs.py
        - mk_spin_plots.py
        - plotopts.dat

8) Sub jobs to run results plots with sub_spin_jobs.py:

  ./sub_spin_jobs -d datacard_spinanalysis.txt -q datafiles_spinanalysis_qqbar.txt -m 125 -D pre_fit_jobs -f 0 -t 100 -n 100 
  ./sub_spin_jobs -d datacard_spinanalysis.txt -q datafiles_spinanalysis_qqbar.txt -m 125 -D post_fit_jobs -f 1 -t 100 -n 100

9) Check and merge jobs with merge_jobs.py:

  ./merge_jobs.py -d pre_fit_jobs -m 125
  ./merge_jobs.py -d post_fit_jobs -m 125

10) Make plots with mk_spin_plots.py. The easiest way is to modify plotopts.dat
