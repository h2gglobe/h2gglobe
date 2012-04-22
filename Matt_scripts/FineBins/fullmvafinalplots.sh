#!/bin/bash
path=$1

for m in {120..130}; do
  for i in {0..9}; do
  echo -------------------------------
  echo Running $2 limits and pvals for mass $m.$i
  echo -------------------------------
  combine $path/mva-datacards-$2-finebins/mva-datacard_$2_$m.$i.txt -M Asymptotic -D data_$2 -m $m.$i --newGenerator=1 -H ProfileLikelihood
  combine $path/mva-datacards-$2-finebins/mva-datacard_$2_$m.$i.txt -M ProfileLikelihood -D data_$2 -m $m.$i -S 1 --rMin=0. --rMax=25 --signif --pvalue
  combine $path/mva-datacards-$2-finebins/mva-datacard_$2_$m.$i.txt -M MaxLikelihoodFit -D data_$2 -m $m.$i -S 1 --minimizerStrategy=2
  done
done

#for m in {110..149}; do 
#  echo -------------------------------
#  echo Running $2 limit and pvals for mass $m.5 
#  echo -------------------------------
#  combine $path/mva-datacards-$2/mva-datacard_$2_$m.5.txt -M Asymptotic -D data_$2 -m $m.5 --newGenerator=1 -H ProfileLikelihood
#  combine $path/mva-datacards-$2/mva-datacard_$2_$m.5.txt -M ProfileLikelihood -D data_$2 -m $m.5 -S 1 --rMin=0. --rMax=25 --signif --pvalue
#  combine $path/mva-datacards-$2/mva-datacard_$2_$m.5.txt -M MaxLikelihoodFit -D data_$2 -m $m.5 -S 1 --minimizerStrategy=2
#done

for type in Asymptotic ProfileLikelihood MaxLikelihoodFit; do
  mkdir -p $path/combineResults_finebins/$2/$type 
  mv higgsCombineTest.${type}*.root $path/combineResults_finebins/$2/$type
  for m in {120..130}; do
    mv $path/combineResults_finebins/$2/$type/higgsCombineTest.$type.mH$m.root $path/combineResults_finebins/$2/$type/higgsCombineTest.$type.mH$m.0.root
  done
done
