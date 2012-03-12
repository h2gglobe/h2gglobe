#!/bin/bash
GRIDFILE=$1
ROOMODEL=$2
WORKINGDIR=$3
QUEUE="8nh"
DATE=`date +%F_%R_%S | tr -s ' ' '_' | tr -s ':' '_'`
if [ "$WORKINGDIR" = "" ]; then
	WORKINGDIR="Frequentist"_$DATE
fi
if [ ! -d $WORKINGDIR ]; then
	mkdir $WORKINGDIR
fi
cp run_OBSFrequentist.sh $WORKINGDIR
cp run_EXPFrequentist.sh $WORKINGDIR
for M in `seq 110 0.5 150`; do
	echo "Submitting mass $M for Frequentist method. Results in $WORKINGDIR"
	bsub -q $QUEUE -o $WORKINGDIR/OBSFrequentist.mH$M.log $PWD/$WORKINGDIR/run_OBSFrequentist.sh $PWD/$WORKINGDIR $M $ROOMODEL $GRIDFILE
	for LIMIT in 0.025 0.160 0.500 0.840 0.975; do
		bsub -q $QUEUE -o $WORKINGDIR/EXPFrequentist.mH$M.$LIMIT.log $PWD/$WORKINGDIR/run_EXPFrequentist.sh $PWD/$WORKINGDIR $M $ROOMODEL $GRIDFILE $LIMIT
	done
done
