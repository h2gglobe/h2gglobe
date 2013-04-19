#!/bin/bash

set -x

source "version.sh"

mkdir -p pileup

dir=/store/group/phys_higgs/cmshgg/reduced/${version}/mc/Summer12_S10_8TeV

cd pileup 
rm *.root
cmsLs $dir | awk '{ print $5}' | grep '/' | grep -v broken > samples.txt 

echo "Hadding and copying"
../parallel --eta --joblog parallel_pileupMerger.log --progress "python ../Macros/pileup/pileupMerger.py --putBack {1}" :::: samples.txt


