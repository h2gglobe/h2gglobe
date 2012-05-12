#!/bin/bash

./submit_reduction.sh data_2012  PhotonPromptReco2012A  25

./submit_reduction.sh mc_sig_summer12_s7 \* 20

./submit_reduction.sh mc_bkg_summer12_s7 GJet_Pt-20\* 20
./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-10To25\*   20
./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-25To250\*  30
./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-250ToInf\* 50
./submit_reduction.sh mc_bkg_summer12_s7 DYJetsToLL\* 50


###### ./submit_reduction.sh data  PhotonPromptReco2012A_30-04-12  25
###### 
###### ./submit_reduction.sh mc_sig_summer12_s7 GluGluToHToGG_M-120 
###### ./submit_reduction.sh mc_sig_summer12_s7 GluGluToHToGG_M-125 
###### 
###### ./submit_reduction.sh mc_bkg_summer12_s7_V13_01_02 DYJets 10
###### ./submit_reduction.sh mc_bkg_summer12_s7_V13_01_02 QCD_Pt_30_80 10
###### 
###### ./submit_reduction.sh mc_bkg_summer12_s7_V13_01_03 GJet_Pt40_ff 10
###### ./submit_reduction.sh mc_bkg_summer12_s7_V13_01_03 GJet_Pt40_pf 10
###### ./submit_reduction.sh mc_bkg_summer12_s7_V13_01_03 GJet_Pt40_pp 10
