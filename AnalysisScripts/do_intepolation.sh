#!/bin/bash

export DISPLAY=""

dir=$(echo $1 | sed 's%/$%%')

cd ../Macros/

ln -s ../AnalysisScripts/$dir

sqrts=""
if( echo $dir | grep 7TeV ); then 
    sqrts='--is7TeV'
fi

./massInterpolator.py $sqrts -i $dir/CMS-HGG --doSmoothing --pval | tee $dir/interpolation.log  &
## ./massInterpolator.py $sqrts -i $dir/CMS-HGG -o $dir/CMS-HGG_runAB -k 0.4362966 --doSmoothing --pval | tee $dir/interpolation_runAB.log & 
## ./massInterpolator.py $sqrts -i $dir/CMS-HGG -o $dir/CMS-HGG_runC  -k 0.5637034 --doSmoothing --pval | tee $dir/interpolation_runC.log &

wait

cd ../AnalysisScripts

## ./hist_combiner.py $dir
./hist_combiner.py $dir
