#!/bin/bash

## ./submit_reduction.sh data_2012  PhotonPromptReco2012A  50 & 
## ./submit_reduction.sh data_2012  PhotonPromptReco2012A_missing  25 & 
./submit_reduction.sh data_2012  DoublePhotonPromptReco2012B  100 & 
./submit_reduction.sh data_2012  DoublePhotonPromptReco2012B_missing  50 & 
## ./submit_reduction.sh data_2012  PhotonReRecoMay23_2012_missing 50 &

## ./submit_reduction.sh mc_sig_summer12_s7 TTH_HToGG_M-130_8TeV 5 &

### ./submit_reduction.sh mc_sig_summer12_s7 \* 5 & 
### 
### 
###  ./submit_reduction.sh mc_bkg_summer12_s7 QCD_Pt-40\* 40 & 
###  ./submit_reduction.sh mc_bkg_summer12_s7 GJet_Pt-40\* 40 &
### 
### ./submit_reduction.sh mc_bkg_summer12_s7 QCD_Pt\* 40 & 
### ./submit_reduction.sh mc_bkg_summer12_s7 GJet\* 40 &
### 
###  ./submit_reduction.sh mc_bkg_summer12_s7_b DiPhotonJets\*   50 & 
###  ./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-10To25\*   20 & 
###  ./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-25To250\*  30 & 
###  ./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-250ToInf\* 20 & 
### 
### ./submit_reduction.sh mc_bkg_summer12_s7 DYJetsToLL\* 70 & 
### 
### ./submit_reduction.sh data_2012  DiPhotonPromptReco2012B  50 & 
### ./submit_reduction.sh data_2012  DiPhotonPromptReco2012B_p2  5 & 

## to resubmit failed jobs
## ./submit_reduction.sh data_2012  DiPhotonPromptReco2012B  50 46 48 &
##  ./submit_reduction.sh data_2012  PhotonPromptReco2012A  50 15 &

## ./submit_reduction.sh mc_bkg_summer12_s7 QCD_Pt-40\*ff\* 40 12 21 24 27 5 6  & 
## ./submit_reduction.sh mc_bkg_summer12_s7 GJet_Pt-40\*pf\* 40 13 28 &

wait

