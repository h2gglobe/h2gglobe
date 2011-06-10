#!/bin/bash
#1 Working Directory
#2 Mass
cd $1
ODIR=$1/$2GeV
if [ ! -d $ODIR ]; then
	mkdir $ODIR
fi
eval `scramv1 runtime -sh`
combine cms-hgg-datacard_parBKG.txt -M HybridNew --generateNuis=0 --generateExt=1 --fitNuisances=1 --testStat Atlas -D data_mass -m $2 -s -1 --generateBinnedWorkaround -S 1 -H ProfileLikelihood --hintStatOnly  --saveHybridResult --saveToys --toysFile HybridToys_mH$2.root  | tee mH$2.HybridNew.$$.log
SEED=`cat mH$2.HybridNew.$$.log | grep '>>> Used OpenSSL to get a really random seed' | awk '{print$10}'`
mv mH$2.HybridNew.$$.log $ODIR/higgsCombineTest.OBSHybridNew.mH$2.$SEED.log
mv higgsCombineTest.HybridNew.mH$2.$SEED.root $ODIR/higgsCombineTest.OBSHybridNew.mH$2.$SEED.root
