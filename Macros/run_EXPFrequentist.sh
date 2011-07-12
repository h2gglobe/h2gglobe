#!/bin/bash
#1 Working Directory
#2 Mass
#3 Expected Limit (0.0275, 0.16, 0.5, 0.84, 0.975)
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
combine ../RooModel/${2}GeVmodel.root -m ${2} -D data_mass -M HybridNew -s -1 --freq --rAbsAcc=0.001 --grid=FrequentistGrid$SUBNUMBER.root --expectedFromGrid $3
limit=`echo "$3" | cut -c 1-5`
limit=`perl -e "printf('%.03f', $limit)"`
filename=`/bin/ls higgsCombineTest.HybridNew.mH$2.[-0-9]*.quant${limit}.root`
mv $filename expectedFrequentist.mH$2.quant${limit}.root
