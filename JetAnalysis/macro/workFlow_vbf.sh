#!/bin/bash
set -x 

json=$1
label=$2

what=all
[[ -n $3 ]] && what=$3

steps="0 1 2"
[[ -n $4 ]] && steps=$(echo $4 | tr ',' ' ')


if [[ $what == train ]] || [[ $what == all ]]; then
    ./classify.py -w workflows_v1 -T vbf_mva.json,$json -l $label | tee $label.log
fi

for bnd in 1 0.8 0.6 0.4 0.3 0.2; do
    blab=$(echo $bnd | sed 's%\.%%')
    webdir=~/www/scratch/higgs/VBF/legacy_training_oct8/workflows_v1_$blab

    mkdir -p $webdir/$label
    
    if [[ $what == cat || $what == all ]]; then
	set $steps
	for step in $@; do
	    mkdir $webdir/$label/step$step
	    ./categoryOptimizationMultiDim.py --label $label -l vbf_mva_categories.json,step$step.json -o $webdir/$label/step$step  --set bound=-$bnd| tee $webdir/$label/step$step/log &
	    mkdir $webdir/$label/step${step}_2D
	    ./categoryOptimizationMultiDim.py --label $label -l vbf_mva_categories.json,step${step}_2D.json -o $webdir/$label/step${step}_2D --set bound=-$bnd  | tee $webdir/$label/step${step}_2D/log &
        ## rm $webdir/$label/step$step/categoryOptimizationMultiDim.root
	    wait
	done
    fi
    
    if [[ $what == plot || $what == cat || $what == all ]]; then
	step=2
	for nc in $(seq 1 4); do 
	    ./plotCategories.py -d $webdir/$label/step${step}/cat_opt.json -o $webdir/$label/step${step}/ -n $nc -i /tmp/$(whoami)/categoryOptimizationMultiDim_${label}.root
	    ./plotCategories.py -d $webdir/$label/step${step}_2D/cat_opt.json -o $webdir/$label/step${step}_2D/ -n $nc -i /tmp/$(whoami)/categoryOptimizationMultiDim_${label}.root
	done
    fi
done
