#!/bin/sh

mkdir outputToy
echo "max events from CRAB: $MaxEvents"
n="$MaxEvents"
./combine model.root -M HybridNew -T $n  --generateNuis=0 --generateExt=1 --fitNuisances=1 --testStat Atlas -D data_mass -m 140 -s -1 -t 1 --generateBinnedWorkaround -S 1 -H ProfileLikelihood --hintStatOnly  --saveHybridResult --saveToys --toysFile HybridToys_mH140.root
rm CMS-HGG.root
rm model.root
mv *.root outputToy/
tar cvfz outputToy.tgz outputToy/
