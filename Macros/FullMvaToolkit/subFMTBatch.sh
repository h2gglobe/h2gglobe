#!/bin/bash
# $1=$PWD $2=filename
cd $1
eval `scramv1 runtime -sh`
./runIt.exe -i $2 --bkgModel --interp --datacards --runCombine --www ~/public_html/h2g/MVA/FMT --diagnose --dumpDatFile mvaanalysis.dat 

