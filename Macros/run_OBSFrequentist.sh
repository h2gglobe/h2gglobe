#!/bin/bash
#1 Working Directory
#2 Mass
REMOTEDIR=$PWD
SUBNUMBER=1
rfcp /castor/cern.ch/user/d/drberry/FrequentistGrid.root .
if [ ! -d $1 ]; then
	mkdir $1
fi
cd $1
ln -s ../CMS-HGG.root CMS-HGG.root
while [ -e FrequentistGrid$SUBNUMBER.root ]; do
	SUBNUMBER=$(($SUBNUMBER+1))
done
ln -s $REMOTEDIR/FrequentistGrid.root FrequentistGrid$SUBNUMBER.root
eval `scramv1 runtime -sh`
combine ../RooModel/${2}GeVmodel.root -m ${2} -D data_mass -M HybridNew --freq --rAbsAcc=0.001 --grid=FrequentistGrid$SUBNUMBER.root
/bin/ls
mv higgsCombineTest.HybridNew.mH$2.root higgsCombineOBSERVED.Frequentist.mH$2.root
