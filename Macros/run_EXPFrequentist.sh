#!/bin/bash
#1 Working Directory
#2 Mass
#3 RooModel Directory
#4 Grid File
#5 Expected Limit (0.0275, 0.16, 0.5, 0.84, 0.975)
REMOTEDIR=$PWD
cd $1
cp ../$4 $REMOTEDIR
eval `scramv1 runtime -sh`
#combine ../${3}/${2}GeVmodel.root -m ${2} -D data_mass -M HybridNew -s -1 --freq --signif --rAbsAcc=0.001 --toysFile=$REMOTEDIR/$4 --expectedFromGrid $5
combine ../${3}/${2}GeVmodel.root -m ${2} -D data_mass -M HybridNew -s -1 --freq --rAbsAcc=0.001 --grid=$REMOTEDIR/$4 --expectedFromGrid $5
limit=`echo "$5" | cut -c 1-5`
limit=`perl -e "printf('%.03f', $limit)"`
SEED=`/bin/ls higgsCombineTest.HybridNew.mH${2}.*.quant${limit}.root | sed "s|higgsCombineTest.HybridNew.mH${2}.\([-0-9][-0-9]*\).quant${limit}.root|\1|" | grep -v .root`
mv higgsCombineTest.HybridNew.mH$2.${SEED}.quant${limit}.root expectedFrequentist.mH$2.quant${limit}.root
