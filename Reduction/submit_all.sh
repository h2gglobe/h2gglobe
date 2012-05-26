#!/bin/bash

./submit_reduction.sh data_2012  PhotonPromptReco\*  50 & 

./submit_reduction.sh mc_sig_summer12_s7 \* 5 & 

./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-10To25\*   20 & 
./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-25To250\*  30 & 
./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-250ToInf\* 20 & 

./submit_reduction.sh mc_bkg_summer12_s7 DYJetsToLL_M-50 50 & 
./submit_reduction.sh mc_bkg_summer12_s7 DYJetsToLL_M-50_52_1e6 5 & 

./submit_reduction.sh mc_bkg_summer12_s7 QCD_Pt\* 40 & 

./submit_reduction.sh mc_bkg_summer12_s7 GJet\* 40 &

wait

