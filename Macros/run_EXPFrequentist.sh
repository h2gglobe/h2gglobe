#!/bin/bash
#1 Working Directory
#2 Mass
#3 Expected Limit (0.0275, 0.16, 0.5, 0.84, 0.975)
cd $1
ODIR=$1
cd ..
if [ ! -d $ODIR ]; then
	mkdir $ODIR
fi
eval `scramv1 runtime -sh`
combine RooModel/${2}GeVmodel.root -m ${2} -D data_mass -M HybridNew --freq --rAbsAcc=0.001 --grid=FrequentistGrid.root --expectedFromGrid $3
set limit = `echo "$3" | cut -c1-5`
limit = `perl -e "printf('%.03f', $limit)"`
mv higgsCombineTest.HybridNew.mH$2.quant${limit}.root $ODIR/expectedFrequentist.mH$2.quant${limit}.root
