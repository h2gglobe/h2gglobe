#!/bin/bash

. setup.sh

rm data2012_prompt/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_03_08/data ${storedir}/data data2012_prompt.txt

rm data2012_rereco/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_03_08/data ${storedir}/data data2012_rereco.txt

rm mc_bkg_summer12_s7/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_03_08/mc/Summer12_S7_8TeV ${storedir}/Summer12_S7_8TeV mc_bkg_summer12_s7.txt

rm mc_sig_summer12_s7/*.dat
./AnalysisScripts/mk_reduction_dat.py /store/group/phys_higgs/cmshgg/processed/V13_03_08/mc/Summer12_S7_8TeV ${storedir}/Summer12_S7_8TeV mc_sig_summer12_s7.txt


tar zcf ${version}.tgz  AnalysisScripts/{common,reduction,aux,*.py}
