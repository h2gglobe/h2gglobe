#!/bin/bash

./submit_reduction.sh data_2012  PhotonPromptReco2012A\*  100

./submit_reduction.sh mc_sig_summer12_s7 \* 20

./submit_reduction.sh mc_bkg_summer12_s7 GJet_Pt-20\* 50
./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-10To25\*   20
./submit_reduction.sh mc_bkg_summer12_s7 DiPhotonBorn_Pt-25To250\*  10
./submit_reduction.sh mc_bkg_summer12_s7 DiPhotonBox_Pt-25To250\*  30
./submit_reduction.sh mc_bkg_summer12_s7 DiPhoton\*_Pt-250ToInf\* 45
./submit_reduction.sh mc_bkg_summer12_s7 DYJetsToLL_M-50 50
./submit_reduction.sh mc_bkg_summer12_s7 DYJetsToLL_M-50_52_1e6 3
./submit_reduction.sh mc_bkg_summer12_s7 QCD_Pt-30to40\* 50
./submit_reduction.sh mc_bkg_summer12_s7 GJet\* 50
