#!/bin/bash

. version.sh

## rm data/*.dat
## ./PhotonAnalysis_scripts/mk_reduction_dat.py /castor/cern.ch/user/c/cmshgg/processed/h2g_V06_00/Data ./datastore/${version}/Data data.txt
## rm  mc_2011pu/*.dat

tar zcf ${version}.tgz  PhotonAnalysis_scripts/*.dat
