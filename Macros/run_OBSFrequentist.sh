#!/bin/bash
#1 Working Directory
#2 Mass
REMOTEDIR=$PWD
rfcp /castor/cern.ch/user/d/drberry/FrequentistGrid.root .
if [ ! -d $1 ]; then
	mkdir $1
fi
cd $1
ln -s ../CMS-HGG.root CMS-HGG.root
ln -s $REMOTEDIR/FrequentistGrid.root FrequentistGrid$$.root
eval `scramv1 runtime -sh`
combine ../RooModel/${2}GeVmodel.root -m ${2} -D data_mass -M HybridNew --freq --rAbsAcc=0.001 --grid=FrequentistGrid$$.root
mv higgsCombineTest.HybridNew.mH$2.root higgsCombineOBSERVED.Frequentist.mH$2.root
