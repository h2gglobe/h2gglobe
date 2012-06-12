#!/bin/bash

. setup.sh

rm data_2012/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_03_05/data ${storedir}/data data_2012.txt

rm mc_sig_summer12_s7/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_03_03/mc/Summer12_S7_8TeV ${storedir}/Summer12_S7_8TeV mc_sig_summer12_s7.txt
### 
### rm mc_bkg_summer12_s7_b/*.dat
### ./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_03_04/mc/Summer12_S7_8TeV ${storedir}/Summer12_S7_8TeV mc_bkg_summer12_s7_b.txt
### 
### rm mc_bkg_summer12_s7/*.dat
### ./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_03_03/mc/Summer12_S7_8TeV ${storedir}/Summer12_S7_8TeV mc_bkg_summer12_s7.txt

tar zcf ${version}.tgz  AnalysisScripts/*.dat
