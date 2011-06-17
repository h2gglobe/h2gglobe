#!/bin/bash
#1 method: ProfileLikelihood or MarkovChainMC or HybridNew  
METHOD=$1
DATE=`date +%F_%R_%S | tr -s ' ' '_' | tr -s ':' '_'`
WORKINGDIR=${METHOD}_$DATE
if [ ! -d $WORKINGDIR ]; then
	mkdir $WORKINGDIR
fi
cp run_OBS$METHOD.sh $WORKINGDIR
for M in `seq 105 0.5 140`; do
	echo "Submitting mass $M for $1 method. Results in $WORKINGDIR"
	bsub -q 8nm $PWD/$WORKINGDIR/run_OBS$METHOD.sh $PWD/$WORKINGDIR $M
done
