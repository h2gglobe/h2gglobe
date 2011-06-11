#!/bin/bash
#1 method: ProfileLikelihood or MarkovChainMC or HybridNew  
METHOD=$1
	for M in `seq 105 0.5 140`; do
		echo "Submitting mass $M for $1 method"
		bsub -q 1nh run_OBS$METHOD.sh $PWD $M 
	done

