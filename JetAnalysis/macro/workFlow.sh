#!/bin/bash

json=$1
label=$2

./classify.py -w workflows_v2 -T vbf_mva.json,$json -l $label | tee $label.log

webdir=~/www/scratch/higgs/VBF/legacy_training/workflows_v5

mkdir -p $webdir/$label

for step in 0 1 2; do
    ./categoryOptimizationMultiDim.py --label $label -l vbf_mva_categories.json,step$step.json -o $webdir/$label/step$step | tee $webdir/$label/step$step/log &
    ./categoryOptimizationMultiDim.py --label $label -l vbf_mva_categories.json,step${step}_2D.json -o $webdir/$label/step${step}_2D | tee $webdir/$label/step${step_2D}/log  &
    ## rm $webdir/$label/step$step/categoryOptimizationMultiDim.root
    wait
done
