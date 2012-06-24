#!/bin/bash

set -x

source "version.sh"

mkdir -p pileup

dir=/store/group/phys_higgs/cmshgg/reduced/${version}/Summer12_S7_8TeV

cd pileup 
rm *.root
cmsLs $dir | awk '{ print $5}' | grep '/' | grep -v broken > samples.txt 

echo "Hadding..."
../parallel --eta --joblog parallel_pileupMerger.log --progress "python ../Macros/pileup/pileupMerger.py {1}" :::: samples.txt

echo "Copying to eos..."
../parallel --eta --joblog parallel_pileupCopy.log --progress "cmsStage -f `basename {1}`.pileup.root $dir" :::: samples.txt

