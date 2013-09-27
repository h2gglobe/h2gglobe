---------------------------------------------------------------
Background Function Choice Profile (Envelope) Fitting for H->gg
---------------------------------------------------------------
Original authors: 
	Matt Kenzie (matthew.william.kenzie@cern.ch)
	Nick Wardle (nckw@cern.ch)

The classes and scripts herein are used for three main purposes:

	1) Background bias studies 
			- where single functions can be tested
		 		against various truth models or the envelope method
				can be tested with different correction schemes etc
			- there are also options to make truth models from
				RooKeysPdf (kernal estimator), regression splines or hybrid 
				pdfs 
	2) Running an fTest
			- this will determine which functions should be included
				as "truth" models for a given dataset
			- it will also determine which functions should be included
				in the envelope for a given dataset
			- it will construct the "MultiPdf" workspace for envelope
				profiling
	3) Making plots of the background model with uncertainty bands
		 which include the function choice uncertainty

Compile with:

	cd h2gglobe/BackgroundProfileFitting
	cmsenv
	make

NOTE: Compilation needs access to $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so

---------------------------------------------------------------
1) Bias studies
---------------------------------------------------------------
This would take a long time to explain and I don't have the time
to write it all down right now. Start by running bin/BiasStudy -h

---------------------------------------------------------------
2) fTest
---------------------------------------------------------------
Figure out appropriate truth functions and make MultiPdf ws. Run
on output globe data file:

	./bin/fTest -i CMS-HGG_data.root -c <ncats> --saveMultiPdf CMS-HGG_multipdf.root

---------------------------------------------------------------
3) Background model plots
---------------------------------------------------------------
These get run by a c++ macro: test/makeBkgPlots.cpp
They will make the background plots and plot the signal model (if a file is passed)
and calculate the bands by profiling over the various functions.

You can pass it any globe bkg model workspace and optionally also a signal
workspace (binned or parametric). It does one category at a time and has a python wrapper 
(scripts/subBkgPlots.py) which submits each category to the batch. It will save the bands, 
datasets, pdfs and canvas in an output file (stored in the directory which you pass) so that
if you want to change a plotting feature you don't have to rerun the bands.

If you want to do a lot of points this is slow (the algorithm used is
a function called guessNew() and it's not very smart or optimal) - please
write a new one if you're good at that sort of thing. If you want to speed things up
then you can relax the nll tolerance or run in wider steps of mass (like 5 or 10).

USE THE BINNED DATA TO SPEED THINGS UP !!!!! (option --useBinnedData)

An example of running this would be:

./script/subBkgPlots.py -b MultiPdfWS.root -s ParametricSignalModel.root -d BkgPlots -c 9 --useBinnedData --isMultiPdf --doBands --massStep 1
