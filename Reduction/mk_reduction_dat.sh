#!/bin/bash

. setup.sh

rm data/*.dat
./PhotonAnalysis_scripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/V11_04/Data ${storedir}/Data data.txt

rm mc/*.dat
./PhotonAnalysis_scripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/V11_04/mc ${storedir}/mc mc.txt

tar zcf ${version}.tgz  PhotonAnalysis_scripts/*.dat
