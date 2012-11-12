#!/bin/bash

parallel=$(dirname $(which $0))/parallel

wd=$1 && shift
datacard=$1 && shift
label=$1 && shift

cd $wd
mkdir $label 

### seq 110 0.5 150 | $parallel -j 8 "combine --verbose=2  -M MaxLikelihoodFit --rMin=-10 --rMax=10. -n $label -m {} -S 1 $datacard | tee pvalue_${label}_{}.log"

seq 120 0.5 130 | $parallel -j 1 "combine --verbose=2  -M MaxLikelihoodFit --rMin=-10 --rMax=10. -n $label -m {} -S 1 $datacard | tee pvalue_${label}_{}.log"

hadd -f $label/higgsCombine$label.MaxLikelihoodFit.root higgsCombine$label.MaxLikelihoodFit.mH*.root 

mv *$label*.root $label

mv *$label*.log $label

## hadd 
