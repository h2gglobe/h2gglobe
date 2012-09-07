#!/bin/bash

## submit data prompt
./submit_reduction.sh     data2012_53_v1          Run2012C-PromptReco-v1  50
./submit_reduction.sh     data2012_53_v1          Run2012C-PromptReco-v2  200
./submit_reduction.sh     data2012_53_v1          Run2012A-13Jul2012-v1         40
./submit_reduction.sh     data2012_53_v1          Run2012A-recover-06Aug2012-v1 10
./submit_reduction.sh     data2012_53_v1          Run2012B-13Jul2012-v1         200

./submit_reduction.sh     mc_bkg_summer12_s10     DYJetsToLL_M-50_TuneZ2Star_8TeV 20

wait

