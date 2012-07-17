#!/bin/bash

## submit data rereco
./submit_reduction.sh   data_2012_rereco   DoublePhoton_Run2012B-29Jun2012-v1_sub1  350
./submit_reduction.sh   data_2012_rereco   DoublePhoton_Run2012B-29Jun2012-v1_sub2  10
./submit_reduction.sh   data_2012_rereco   Photon_Run2012A-29Jun2012-v1_sub1  130
./submit_reduction.sh   data_2012_rereco   Photon_Run2012A-29Jun2012-v1_sub2  1
./submit_reduction.sh   data_2012_rereco   Photon_Run2012A-recover_29Jun2012-v1_sub1  15
./submit_reduction.sh   data_2012_rereco   Photon_Run2012A-recover_29Jun2012-v1_sub2  1


wait

