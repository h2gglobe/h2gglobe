#!/bin/bash
path=$1
b=$3
cd $4

for m in {110..150}; do 
  echo -------------------------------
  echo Running $2 limits and pvals for mass $m.0
  echo -------------------------------
  #combine $path/mva-datacards-$2-novbf/mva-datacard_$2_$m.0.txt -M Asymptotic -D data_$2 -m $m --newGenerator=1 -H ProfileLikelihood
  combine $path/mva-datacards-$2-vbf/mva-datacard_$2_$m.0.txt -M ProfileLikelihood -D data_$2 -m $m -S 1 --rMin=0. --rMax=25 --signif --pvalue
  #combine $path/mva-datacards-$2-bin$b/mva-datacard_$2_$m.0.txt -M MaxLikelihoodFit -D data_$2 -m $m -S 1 --minimizerStrategy=2
done

for m in {110..149}; do 
  echo -------------------------------
  echo Running $2 limit and pvals for mass $m.5 
  echo -------------------------------
  #combine $path/mva-datacards-$2-novbf/mva-datacard_$2_$m.5.txt -M Asymptotic -D data_$2 -m $m.5 --newGenerator=1 -H ProfileLikelihood
  combine $path/mva-datacards-$2-vbf/mva-datacard_$2_$m.5.txt -M ProfileLikelihood -D data_$2 -m $m.5 -S 1 --rMin=0. --rMax=25 --signif --pvalue
  #combine $path/mva-datacards-$2-bin$b/mva-datacard_$2_$m.5.txt -M MaxLikelihoodFit -D data_$2 -m $m.5 -S 1 --minimizerStrategy=2
done

