#!/bin/bash
#1 Working Directory
#2 Mass
#3 Expected Limit (0.0275, 0.16, 0.5, 0.84, 0.975)
REMOTEDIR=$PWD
rfcp /castor/cern.ch/user/d/drberry/FrequentistGrid.root .
if [ ! -d $1 ]; then
	mkdir $1
fi
cd $1
ln -s ../CMS-HGG.root CMS-HGG.root
ln -s $REMOTEDIR/FrequentistGrid.root FrequentistGrid$$.root
eval `scramv1 runtime -sh`
combine ../RooModel/${2}GeVmodel.root -m ${2} -D data_mass -M HybridNew --freq --rAbsAcc=0.001 --grid=FrequentistGrid$$.root --expectedFromGrid $3
limit=`echo "$3" | cut -c 1-5`
limit=`perl -e "printf('%.03f', $limit)"`
mv higgsCombineTest.HybridNew.mH$2.quant${limit}.root expectedFrequentist.mH$2.quant${limit}.root
