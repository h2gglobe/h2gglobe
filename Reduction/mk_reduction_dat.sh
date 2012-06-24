#!/bin/bash

. setup.sh

rm data_2012/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_03_06/dataichep2012 ${storedir}/data data_2012.txt

rm mc_bkg_summer12_s7/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_03_06/mc/Summer12_S7_8TeV ${storedir}/Summer12_S7_8TeV mc_bkg_summer12_s7.txt

rm mc_sig_summer12_s7_03_06/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_03_06/mc/Summer12_S7_8TeV ${storedir}/Summer12_S7_8TeV mc_sig_summer12_s7_03_06.txt

tar zcf ${version}.tgz  AnalysisScripts/{common,reduction,aux,*.py}
