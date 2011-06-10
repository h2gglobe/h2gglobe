#!/bin/bash
#Method (ProfileLikelihood, MarkovChainMC, HybridNew)
for M in `seq 105 0.5 140`; do
	if [ $1="ProfileLikelihood" ]; then
		echo "Submitting mass $M for $1 method"
		bsub -q 1nh run_OBSProfileLikelihood.sh $PWD $M
	fi
	if [ $1="MarkovChainMC" ]; then
		echo "Submitting mass $M for $1 method"
		bsub -q 1nh run_OBSMarkovChainMC.sh $PWD $M
	fi
	if [ $1="HybridNew" ]; then
		echo "Submitting mass $M for $1 method"
		bsub -q 1nh run_OBSHybridNew.sh $PWD $M
	fi
done
