#!/bin/bash

. setup.sh

rm data/*.dat
./PhotonAnalysis_scripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/V13_01_01/data ${storedir}/Data data.txt

rm  mc_sig_summer12_s7/*.dat
./PhotonAnalysis_scripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/V13_01_02/mc/Summer12_S7 ${storedir}/MC_Sig_Summer12_S7 mc_sig_summer12_s7.txt

rm  mc_bkg_summer12_s7/*.dat
./PhotonAnalysis_scripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/V13_01_02/mc/Summer12_S7 ${storedir}/MC_Bkg_Summer12_S7 mc_bkg_summer12_s7.txt

tar zcf ${version}.tgz  PhotonAnalysis_scripts/*.dat
