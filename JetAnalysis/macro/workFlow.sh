#!/bin/bash

json=$1
label=$2

## ./classify.py -w workflows_v2 -T vbf_mva.json,$json -l $label | tee $label.log

webdir=~/www/scratch/higgs/VBF/legacy_training/workflows_v5

mkdir -p $webdir/$label

## for step in 0 1 2 3 4 5 1_2D 2_2D 3_2D 4_2D 5_2D; do
## for step in 1 5 4_2D 5_2D; do
## for step in 2 3 4 5; do
## for step in 2_2D 4_2D; do
## for step in 4_2D; do
## for step in 2 3 4; do
## for step in 3 4; do
## for step in 2; do
# for step in 3_2D; do
## for step in 0 1 0_2D 1_2D; do
## for step in 0 1 2; do
for step in 1 2; do
    ## ./categoryOptimizationMultiDim.py --label $label -l vbf_mva_categories.json,step$step.json -o $webdir/$label/step$step | tee $webdir/$label/step$step/log &
    ./categoryOptimizationMultiDim.py --label $label -l vbf_mva_categories.json,step${step}_2D.json -o $webdir/$label/step${step}_2D | tee $webdir/$label/step${step_2D}/log  &
    ## rm $webdir/$label/step$step/categoryOptimizationMultiDim.root
    wait
done
