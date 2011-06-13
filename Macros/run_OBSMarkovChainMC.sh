#!/bin/bash
#1 Working Directory
#2 Mass
cd $1
ODIR=$1/MarkovChainMC
if [ ! -d $ODIR ]; then
	mkdir $ODIR
fi
eval `scramv1 runtime -sh`
combine cms-hgg-datacard_parBKG.txt -M MarkovChainMC -D data_mass -m $2 -s -1 -S 1 --generateBinnedWorkaround  -H ProfileLikelihood -n OBSERVED --rMin=0 --rMax=40 --truncatedMeanFraction 0.05  --iteration 75000 --tries 100 --hintStatOnly | tee mH$2.Markov.$$.log
SEED=`cat mH$2.Markov.$$.log | grep '>>> Used OpenSSL to get a really random seed' | awk '{print$10}'`
mv mH$2.Markov.$$.log $ODIR/higgsCombineOBSERVED.MarkovChainMC.mH$2.$SEED.log
mv higgsCombineOBSERVED.MarkovChainMC.mH$2.$SEED.root $ODIR/higgsCombineOBSERVED.MarkovChainMC.mH$2.root


