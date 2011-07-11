#!/bin/bash
DATE=`date +%F_%R_%S | tr -s ' ' '_' | tr -s ':' '_'`
WORKINGDIR="Frequentist"_$DATE
if [ ! -d $WORKINGDIR ]; then
	mkdir $WORKINGDIR
fi
cp run_OBSFrequentist.sh $WORKINGDIR
cp run_EXPFrequentist.sh $WORKINGDIR
for M in `seq 110 0.5 140`; do
	echo "Submitting mass $M for Frequentist method. Results in $WORKINGDIR"
	bsub -q 1nh -o $WORKINGDIR/OBSFrequentist.mH$M.log $PWD/$WORKINGDIR/run_OBSFrequentist.sh $PWD/$WORKINGDIR $M
	bsub -q 1nh -o $WORKINGDIR/EXPFrequentist.mH$M.0.027.log $PWD/$WORKINGDIR/run_EXPFrequentist.sh $PWD/$WORKINGDIR $M 0.0275
	bsub -q 1nh -o $WORKINGDIR/EXPFrequentist.mH$M.0.160.log $PWD/$WORKINGDIR/run_EXPFrequentist.sh $PWD/$WORKINGDIR $M 0.160
	bsub -q 1nh -o $WORKINGDIR/EXPFrequentist.mH$M.0.500.log $PWD/$WORKINGDIR/run_EXPFrequentist.sh $PWD/$WORKINGDIR $M 0.500
	bsub -q 1nh -o $WORKINGDIR/EXPFrequentist.mH$M.0.840.log $PWD/$WORKINGDIR/run_EXPFrequentist.sh $PWD/$WORKINGDIR $M 0.840
	bsub -q 1nh -o $WORKINGDIR/EXPFrequentist.mH$M.0.975.log $PWD/$WORKINGDIR/run_EXPFrequentist.sh $PWD/$WORKINGDIR $M 0.975
done
