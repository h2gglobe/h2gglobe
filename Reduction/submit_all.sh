#!/bin/bash

## submit data prompt
./submit_reduction.sh     data2012_HCP          DoublePhoton_Run2012B-13Jul2012-v1          200
./submit_reduction.sh     data2012_HCP          DoublePhoton_Run2012B-13Jul2012-v1_sub2     5
./submit_reduction.sh     data2012_HCP          DoublePhoton_Run2012B-13Jul2012-v1_sub3     10
./submit_reduction.sh     data2012_HCP          DoublePhoton_Run2012B-13Jul2012-v1_sub4     3
./submit_reduction.sh     data2012_HCP          DoublePhoton_Run2012C-PromptReco-v2         250
./submit_reduction.sh     data2012_HCP          DoublePhoton_Run2012C-PromptReco-v2_sub2    8
./submit_reduction.sh     data2012_HCP          DoublePhoton_Run2012C-PromptReco-v2_sub4    100
./submit_reduction.sh     data2012_HCP          DoublePhoton_Run2012C-PromptReco-v2_sub5    1
./submit_reduction.sh     data2012_HCP          Photon_Run2012A-13Jul2012-v1                30
./submit_reduction.sh     data2012_HCP          Photon_Run2012A-13Jul2012-v1_sub2           4
./submit_reduction.sh     data2012_HCP          Photon_Run2012A-13Jul2012-v1_sub3           1
./submit_reduction.sh     data2012_HCP          Photon_Run2012A-recover-06Aug2012-v1        5
./submit_reduction.sh     data2012_HCP          Photon_Run2012A-recover-06Aug2012-v1__sub1  1
./submit_reduction.sh     data2012_HCP          Run2012C-24Aug2012-v2                       30

./submit_reduction.sh     mc_sig_summer12_s10     \* 4

wait

