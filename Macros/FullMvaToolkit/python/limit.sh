#!/bin/bash
# run with $1= directory $2=boost type (grad / ada) $3=$PWD
cd $3
source /vols/cms/grid/setup.sh
eval `scramv1 runtime -sh`
cd $1

for m in {110..150}; do 
  echo -------------------------------
  echo Running $2 $1 limit for mass $m.0
  echo -------------------------------
  combine mva-datacard_$2_$m.0.txt -M Asymptotic -D data_$2 -m $m 
  #combine mva-datacard_$2_$m.0.txt -M AsymptoticNew -D data_$2 -m $m --nPoints=20 --scanMax=4.
  combine mva-datacard_$2_$m.0.txt -M MaxLikelihoodFit -D data_$2 -m $m --minimizerStrategy=2
  combine mva-datacard_$2_$m.0.txt -M ProfileLikelihood -m $m -S 1 --rMin=0. --rMax=25 --signif --pvalue
done

for m in {110..149}; do 
  echo -------------------------------
  echo Running $2 $1 limit for mass $m.5 
  echo -------------------------------
  combine mva-datacard_$2_$m.5.txt -M Asymptotic -D data_$2 -m $m.5 
  #combine mva-datacard_$2_$m.5.txt -M AsymptoticNew -D data_$2 -m $m.5 --nPoints=20 --scanMax=4.
  combine mva-datacard_$2_$m.5.txt -M MaxLikelihoodFit -D data_$2 -m $m.5 --minimizerStrategy=2
  combine mva-datacard_$2_$m.5.txt -M ProfileLikelihood -m $m.5 -S 1 --rMin=0. --rMax=25 --signif --pvalue
done

methods=( Asymptotic AsymptoticNew MaxLikelihoodFit ProfileLikelihood)
for meth in ${methods[@]}; do mkdir $meth; done

for m in {110..150}; do
	for meth in ${methods[@]}; do
    mv higgsCombineTest.$meth.mH$m.root $meth/higgsCombineTest.$meth.mH$m.0.root
  done
done

for m in {110..149}; do
	for meth in ${methods[@]}; do
    mv higgsCombineTest.$meth.mH$m.5.root $meth/higgsCombineTest.$meth.mH$m.5.root
  done
done

mkdir ExpProfileLikelihood
for m in {110..150}; do 
  echo -------------------------------
  echo Running  limit for mass $m.0
  echo -------------------------------
  combine mva-datacard_$2_$m.0.txt -M ProfileLikelihood -m $m -S 1 --rMin=0. --rMax=25 --signif --pvalue -t -1 --expectSignal=1
  mv higgsCombineTest.ProfileLikelihood.mH$m.root ExpProfileLikelihood/higgsCombineTest.ProfileLikelihood.mH$m.0.root
done

for m in {110..149}; do 
  echo -------------------------------
  echo Running $2 $1 limit for mass $m.5 
  echo -------------------------------
  combine mva-datacard_$2_$m.5.txt -M ProfileLikelihood -m $m.5 -S 1 --rMin=0. --rMax=25 --signif --pvalue -t -1 --expectSignal=1
  mv higgsCombineTest.ProfileLikelihood.mH$m.5.root ExpProfileLikelihood/higgsCombineTest.ProfileLikelihood.mH$m.5.root
done

cd $3

