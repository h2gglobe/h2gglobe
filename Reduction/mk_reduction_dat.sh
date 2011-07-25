#!/bin/bash

. setup.sh

rm data/*.dat
./PhotonAnalysis_scripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/h2g_V10_00/Data ${storedir}/Data data.txt

tar zcf ${version}.tgz  PhotonAnalysis_scripts/*.dat
