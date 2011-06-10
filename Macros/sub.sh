#!/bin/bash
#1 Number of Toys
#2 The number of parallel jobs
NUMJOBS=0
if [ -n $2 ]; then
	NUMJOBS=$2
fi 
if [ $NUMJOBS -gt 0 ]
then
	for i in `seq 1 $NUMJOBS`; do 
		for M in {105,110,115,120,125,130,135,140}; do
			echo "Submitting job number $i with mass $M and $1 Toys"
			bsub -q 1nh run_ProfileLikelihood.sh $PWD $M $1
		done
	done
else
	for M in `seq 105 0.5 140`; do
		echo "Submitting mass $M"
		bsub -q 1nh run_OBSProfileLikelihood.sh $PWD $M
	done
fi
