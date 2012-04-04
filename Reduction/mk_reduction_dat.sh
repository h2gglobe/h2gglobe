#!/bin/bash

. setup.sh

rm data/*.dat
./PhotonAnalysis_scripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/V11_04_05/Data ${storedir}/Data data.txt

rm  mc_sig_fall11_s6/*.dat
./PhotonAnalysis_scripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/V11_04_05/Data ${storedir}/MC_Sig_Fall11_S6 mc_sig_fall11_s6.txt

rm  mc_bkg_fall11_s6/*.dat
./PhotonAnalysis_scripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/V11_04_05/Data ${storedir}/MC_Bkg_Fall11_S6 mc_bkg_fall11_s6.txt

tar zcf ${version}.tgz  PhotonAnalysis_scripts/*.dat
