#!/bin/sh

mkdir outputToy
echo "max events from CRAB: $MaxEvents"
n="$MaxEvents"
./combine model.root -M ProfileLikelihood -D data_mass -m 135 -t $n -s -1 -S 1 --generateBinnedWorkaround --tries 1  --maxTries 200 --rMin=0. --rMax=35 -H ProfileLikelihood --hintStatOnly
rm CMS-HGG.root
rm model.root
mv *.root outputToy/
tar cvfz outputToy.tgz outputToy/
