#!/bin/bash
#1 Working Directory
#2 Mass
cd $1
ODIR=$1
cd ..
if [ ! -d $ODIR ]; then
	mkdir $ODIR
fi
eval `scramv1 runtime -sh`
combine RooModel/${2}GeVmodel.root -m ${2} -D data_mass -M HybridNew --freq --rAbsAcc=0.001 --grid=FrequentistGrid.root
mv higgsCombineTest.HybridNew.mH$2.root $ODIR/higgsCombine.Frequentist.mH$2.root
