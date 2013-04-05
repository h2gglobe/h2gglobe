#!/bin/bash

## set -x

inputd1=$1 && shift
inputd2=$1 && shift

jsons=$1 && shift

name1=$(basename $inputd1)
name2=$(basename $inputd2)

ws1=$inputd1/hgg.inputbkgdata_8TeV.root
sig1=$PWD/$inputd1/hgg.inputsig.root

ws2=$inputd2/CMS-HGG.root
sig2=$PWD/$inputd2/CMS-HGG_interpolated.root

input1=$inputd1/MITdump_ICHEP_v2.root
input2=$inputd2/eventList_ICHEP53x.root


for json in $(echo $jsons | tr ',' '\n'); do
    
    echo "Setting up $dir"
    
    dir=$(echo $json | sed 's%.json%%')
    mkdir $dir
    mkdir $dir/$name1
    mkdir $dir/$name2
    
    echo "Splitting minitrees"
    ./splitMiniTree.py -i $input1 -w $ws1 -m -p $json -d $dir/$name1 > $dir/$name1/split.log
    npart=$(awk '/npart:/ { print $2 }' $dir/$name1/split.log)
    
    ./splitMiniTree.py -i $input2 -w $ws2 -g -p $json -d $dir/$name2 > $dir/$name2/split.log
    
    echo "Setting up for combine"
    cp -p $inputd1/datacard.txt $dir/$name1
    cp -p $inputd2/datacard.txt $dir/$name2

    ./setupForCombine.sh $dir/$name1 $(basename $input1 | sed 's%.root%%')_part $npart datacard.txt $sig1
    ./setupForCombine.sh $dir/$name2 $(basename $input2 | sed 's%.root%%')_part $npart datacard.txt $sig2

done

