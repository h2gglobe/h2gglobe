#!/bin/sh

mkdir outputToy
echo "max events from CRAB: $MaxEvents"
n="$MaxEvents"
./combine 105GeVmodel_3sigma.root -M ProfileLikelihood -D data_mass -m 105 -t $n -s -1 -S 1 --generateBinnedWorkaround --tries 1  --maxTries 200 --rMin=0. --rMax=35 -H ProfileLikelihood --hintStatOnly
rm CMS-HGG.root
rm 105GeVmodel_3sigma.root
mv *.root outputToy/
rm *.root
tar cvfz outputToy.tgz outputToy/
