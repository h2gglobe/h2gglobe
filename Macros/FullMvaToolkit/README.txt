******************************
**** FULMVATOOLKIT advice ****
******************************
Author:   Matthew Kenzie
Email:    matthew.william.kenzie@cern.ch
Modified: 29.05.12

- More detailed instructions can be found in fullmvatoolkit_25_05_12.pdf
- Also see:
  https://twiki.cern.ch/twiki/bin/viewauth/CMS/Higgs2GAnalyzer#FullMvaToolkit

- You should run inside a CMSSW_5_2_X release (testing was done with CMSSW_5_2_4_patch4)
- It is probably best to run the rebinning and fitting first and save the fits and bin edges:
  
  ./runIt.exe -i <file> -D -F mvaanalysis.dat
- 
You can then check if the fits look ok (in plots) and that the bin edges are sensible (in mvaanalysis.dat)
- You can then run the background model and interpolation etc. after these look ok:
  
  ./runIt.exe -i <file> -N -b -I -d -D -w <web_dir>

- Run the combine code last:

  ./runIt.exe -i <file> -N -C

- If you feel confident you can execute this in one go using:

  ./subFMTBatch.sh $PWD <file>
