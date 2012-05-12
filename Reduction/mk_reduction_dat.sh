#!/bin/bash

. setup.sh

### rm data/*.dat
### ./AnalysisScripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/V13_01_01/data ${storedir}/Data data.txt
### 
### rm  mc_sig_summer12_s7/*.dat
### ./AnalysisScripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/V13_01_02/mc/Summer12_S7 ${storedir}/MC_Sig_Summer12_S7 mc_sig_summer12_s7.txt
### 
### rm  mc_bkg_summer12_s7_V13_01_02/*.dat
### ./AnalysisScripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/V13_01_02/mc/Summer12_S7 ${storedir}/MC_Bkg_Summer12_S7 mc_bkg_summer12_s7_V13_01_02.txt
### 
### rm  mc_bkg_summer12_s7_V13_01_03/*.dat
### ./AnalysisScripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/V13_01_03/mc/Summer12_S7 ${storedir}/MC_Bkg_Summer12_S7 mc_bkg_summer12_s7_V13_01_03.txt

rm data_2012/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_02_00/data ${storedir}/data data_2012.txt

rm mc_sig_summer12_s7/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_02_00/mc/Summer12_S7_8TeV ${storedir}/Summer12_S7_8TeV mc_sig_summer12_s7.txt

rm mc_bkg_summer12_s7/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_02_00/mc/Summer12_S7_8TeV ${storedir}/Summer12_S7_8TeV mc_bkg_summer12_s7.txt

tar zcf ${version}.tgz  AnalysisScripts/*.dat
