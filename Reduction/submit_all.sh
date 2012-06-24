#!/bin/bash

./submit_reduction.sh data_2012  PhotonPromptReco2012A      200 &
./submit_reduction.sh data_2012  PhotonPromptReco2012A_sub2 30  &

./submit_reduction.sh data_2012  PhotonReRecoMay23_2012          20 &
./submit_reduction.sh data_2012  PhotonReRecoMay23_2012_missing2 1  & 

./submit_reduction.sh data_2012  DoublePhotonPromptReco2012B      500 &
./submit_reduction.sh data_2012  DoublePhotonPromptReco2012B_sub2 15 &

#### ./submit_reduction.sh data_2012  DoublePhotonPromptReco2012B_TopUp_sub1 200
#### ./submit_reduction.sh data_2012  PhotonJune8ReReco2012A_TopUp_sub1 60



################ ./submit_reduction.sh data_2012  DoublePhotonPromptReco2012B_TopUp1_glidein 50
################ ./submit_reduction.sh data_2012  PhotonPromptReco2012A  50 & 
################ ./submit_reduction.sh data_2012  PhotonPromptReco2012A_missing  10 & 
################ ./submit_reduction.sh data_2012  PhotonPromptReco2012A_missing_try2  5 & 
################ 
################ ./submit_reduction.sh data_2012  DoublePhotonPromptReco2012B  200 & 
################ ./submit_reduction.sh data_2012  DoublePhotonPromptReco2012B_missing  50 & 
################ ./submit_reduction.sh data_2012  DoublePhotonPromptReco2012B_missing_try2 25 & 
################ 
################ ./submit_reduction.sh data_2012  PhotonReRecoMay23_2012_missing 44 &
################ ./submit_reduction.sh data_2012  PhotonReRecoMay23_2012_missing_try2 1 &
################ 
################ 
./submit_reduction.sh mc_sig_summer12_s7_03_06 \* 5 &

./submit_reduction.sh mc_bkg_summer12_s7 QCD_Pt\* 40 & 
./submit_reduction.sh mc_bkg_summer12_s7 GJet\* 40 &

./submit_reduction.sh mc_bkg_summer12_s7 DiPhotonJets\*   150 & 
./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-10To25\*   20 & 
./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-25To250\*  30 & 
./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-250ToInf\* 20 & 

./submit_reduction.sh mc_bkg_summer12_s7 DYJetsToLL_M-50_TuneZ2Star_8TeV 300 & 

## to resubmit failed jobs
## ------------------------

### ./submit_reduction.sh mc_sig_summer12_s7_03_05 WH_ZH_HToGG_M-140_8TeV_lepTag 5 1
### 
### ./submit_reduction.sh mc_sig_summer12_s7_03_05 WH_ZH_HToGG_M-125_8TeV_lepTag 5 4
### ./submit_reduction.sh mc_sig_summer12_s7_03_05 WH_ZH_HToGG_M-124_8TeV_lepTag 5 3
### ./submit_reduction.sh mc_sig_summer12_s7_03_05 VBF_HToGG_M-145_8TeV_lepTag 5 0

### ./submit_reduction.sh  data_2012 PhotonPromptReco2012A 50 49
### ./submit_reduction.sh  mc_sig_summer12_s7 TTH_HToGG_M-140_8TeV.dat 5 0

#### ./submit_reduction.sh  mc_sig_summer12_s7_03_06 GluGluToHToGG_M-120_8TeV 5 2

### ./submit_reduction.sh  mc_sig_summer12_s7_03_05 VBF_HToGG_M-145_8TeV_lepTag 5 0

### ./submit_reduction.sh  data_2012 DoublePhotonPromptReco2012B_missing 50 13 20 & 
### 
### ./submit_reduction.sh data_2012  DoublePhotonPromptReco2012B  200 88 & 

## ./submit_reduction.sh data_2012  DiPhotonPromptReco2012B  50 46 48 &
##  ./submit_reduction.sh data_2012  PhotonPromptReco2012A  50 15 &

## ./submit_reduction.sh mc_bkg_summer12_s7 QCD_Pt-40\*ff\* 40 12 21 24 27 5 6  & 
## ./submit_reduction.sh mc_bkg_summer12_s7 GJet_Pt-40\*pf\* 40 13 28 &

wait

