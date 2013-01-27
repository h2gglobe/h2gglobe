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

