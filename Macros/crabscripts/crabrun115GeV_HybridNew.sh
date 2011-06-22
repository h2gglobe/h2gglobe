#!/bin/sh

mkdir outputToy
echo "max events from CRAB: $MaxEvents"
n="$MaxEvents"
./combine 115GeVmodel.root -M HybridNew -T $n  --generateNuis=0 --generateExt=1 --fitNuisances=1 --testStat Atlas -D data_mass -m 115 -s -1 -t 1 --generateBinnedWorkaround -S 1 -H ProfileLikelihood --hintStatOnly
rm CMS-HGG.root
rm 115GeVmodel.root
mv *.root outputToy/
rm *.root
tar cvfz outputToy.tgz outputToy/
