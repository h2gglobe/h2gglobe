#!/bin/bash 

set -x

steps="$(echo $1 | tr ':' ' ')"
shift

label=$1
shift

args=$@

## webdir=~/www/higgs/Legacy/mva_boundaries/v3
webdir=~/www/higgs/Legacy/mva_boundaries/v2
## webdir=~/www/higgs/Legacy/mva_boundaries/v2_vbftree_nodijet

mkdir $webdir

set $steps
for step in $@; do
    out=$webdir/${step}_${label}
    mkdir $out
    ./categoryOptimizationMultiDim.py -l dipho_mva_categories.json,dipho_mva_${step}.json -i ~/work/opttree_8TeV.root -o $out --set vbftreeoff=_ --set opttreeoff= --label $label $args | tee $out/log
    ## ./categoryOptimizationMultiDim.py -l dipho_mva_categories_vbftrees.json,dipho_mva_${step}.json -i ~/work/histograms_vbf_optimization_nov2_full_fix.root -o $out --set vbftreeoff=_ --set opttreeoff= --label $label $args | tee $out/log
done

