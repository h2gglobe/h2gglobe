#!/bin/bash

## ./submit_reduction.sh data_2012  PhotonPromptReco2012A  50 & 
./submit_reduction.sh data_2012  DiPhotonPromptReco2012B  50 & 

### ./submit_reduction.sh mc_sig_summer12_s7 \* 1 & 
### 
### ./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-10To25\*   20 & 
### ./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-25To250\*  30 & 
### ./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-250ToInf\* 20 & 
### 
### ./submit_reduction.sh mc_bkg_summer12_s7 DYJetsToLL\* 70 & 
### 
### ./submit_reduction.sh mc_bkg_summer12_s7 QCD_Pt\* 40 & 
### ./submit_reduction.sh mc_bkg_summer12_s7 GJet\* 40 &

wait

