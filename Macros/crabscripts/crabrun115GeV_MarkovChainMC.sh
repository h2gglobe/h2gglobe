#!/bin/sh

mkdir outputToy
echo "max events from CRAB: $MaxEvents"
n="$MaxEvents"
./combine model.root -M MarkovChainMC -D data_mass -m 115 -t $n -s -1 --generateBinnedWorkaround -S 1 -H ProfileLikelihood --rMin=0 --rMax=40 --iteration 50000 --tries 25 --hintStatOnly
rm CMS-HGG.root
rm model.root
mv *.root outputToy/
tar cvfz outputToy.tgz outputToy/
