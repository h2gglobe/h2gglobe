-------------------------------------
Simultaneous Signal Fitting for H->gg
-------------------------------------
Original author: Matt Kenzie (matthew.william.kenzie@cern.ch)

This package implements fitting signal shapes to signal MC in the H->gg analysis

To guarantee smooth interpolation of the signal parameters the fits are done simultaneously across masses
with the function parameters all taking polynomial forms as a function of mh.

The class which handles the fitting is SimultaneousFit.cc. There are two implementations:
  
  --simFit - This constructs a RooSimultanteous and fits each mc mass together with a sum of gaussians
  --mhFit  - This has an explicit dependedce on mh in the fit (as well as the usual m_{gammagamma} dependence)
           - In essesnce is a 2D fit but it doesn't work as well

Both of these fit a sum of gaussians (which can be user defined) where the parameters (mean, sigma of each gaussian 
and the f the frac of that gaussian) have a polynomial dependence on mh (which can also be user defined).

After compiling (cmsenv; make) you can run ./bin/SimultaneousSignalFit --help to see a list of options.

-------------------------------------
Runnnig the code
-------------------------------------

The run configuration is set in dat/config.dat
It lists an input file (this should be the output workspace of the MC datasets after having run globe)
It also lists a set of jobs with the outfilename, the process, the category, the number of gaussians 
and the polynomial orders of the mean, sigma and fraction separated by a single space

To submit the jobs run:

./scripts/makeJobs.py -i dat/config.dat -o $PWD/files -f 8 -q 1nh --submit

This will make a .sh script for each outfile listed in config.dat, put it in scripts/jobs and submit it to lxbatch.
It will also make a file scripts/jobs/filestocombine.dat which tells a later script how to package up the jobs.

You can check the status with:

./scripts/checkJobs.py scripts/jobs

They should complete in about 15 mins

When they are complete the output will be packaged into a single workspace called CMS-HGG_sigfit.root
This is one can be used with combine

-------------------------------------
Details
-------------------------------------

You can set the starting values of the polynomial paramaters in dat/variables_map.dat
The nuisance parameters are implemented with a constant smearing value subtracted first. 
These values are defined in dat/smear_vals.dat

The fits are run by the utility script test/SimultaneousSignalFit.cpp
  
  -- the fitting is done by the function SimultaneousFit::runFit()
  -- an independent fit at each mass is done first to find nice starting values (this is highly recommended)
  -- then either the simFit or mhFit is run
  -- at the end a finalPdf is built from the fit results as various attributes of the pdfs need to be changed (e.g including nuisance
     params) 
  -- the number of entries, eff*acc and expected entries for 1pb are calculated as functions of mh (using further simple polynomial fits)
     which allows us to define how many sm signal events we expect for a given mh

The output is packaged at the end by test/PackageOutput.cpp

  -- in general this will try and sum up everything that is there (regardless of failed jobs, or jobs that haven't been run)
  -- it should warn you if a pdf or a dataset is missing
  -- this also combines datasets and pdfs across categories and processes which allow extraction of quantities like FWHM and sig_eff
     using the globe script Macros/makeParametricSignalModelPlots.C

Additional classes can be added to interface and src and they should be automatically included in the compile.
Additional scripts can be added to test (under extension .cpp) and they should be automatically included in the compile.

One obvious place for development is to expand the functions SimultaneousFit::getSumOfGaussians to have the ability to use
other shapes (e.g. CBxGaus or whatever)


