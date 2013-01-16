Some instructions for running combine from globe workspaces
and producing final plots etc.

Example (binned mass factorized)

First make a working directory and then even another one inside that if you
are running multiple types (i.e. both the massfac and the cutbased at the same
time)

    mkdir MORIONDreview
    mkdir MORIONDreview/binned_massfac

Copy the relevant workspaces and txt files into that directory

    cp <location_of_massfac_binned_workspace_2011> binned_massfac/CMS-HGG_2011.root
    cp <location_of_massfac_binned_workspace_2012> binned_massfac/CMS-HGG_2012.root
    cp <location_of_massfac_binned_datacard_2011> binned_massfac/datacard_massfac_2011.txt
    cp <location_of_massfac_binned_datacard_2012> binned_massfac/datacard_massfac_2012.txt

Change into the directory and create soft links to the scripts required

    cd MORIONDreview/binned_massfac
    ln -s ../../subBinnedToBatch.py
    ln -s ../../subParametricToBatch.py
    ln -s ../../runFinalResults.py

Run the script "runFinalResults.py" with the relevant options. First with the
--makePlots option off.

    ./runFinalResults.py --card2011="datacard_massfacmva_2011.txt" --card2012="datacard_massfacmva_2012.txt" --unblind

When this has completed run the script again with the --makePlots options on.
You can also pass the lumi of the two years

    ./runFinalResults.py --card2011="datacard_massfacmva_2011.txt" --card2012="datacard_massfacmva_2012.txt" --unblind --makePlots --lumi2011=5.3 --lumi2012=19.6

The plots should appear in the directory you are in under a sub directory
called plots

