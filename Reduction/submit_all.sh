#!/bin/bash

./submit_reduction.sh data  PhotonPromptReco2012A_30-04-12  25

./submit_reduction.sh mc_sig_summer12_s7 GluGluToHToGG_M-120 
./submit_reduction.sh mc_sig_summer12_s7 GluGluToHToGG_M-125 

./submit_reduction.sh mc_bkg_summer12_s7_V13_01_02 DYJets 10
./submit_reduction.sh mc_bkg_summer12_s7_V13_01_02 QCD_Pt_30_80 10

./submit_reduction.sh mc_bkg_summer12_s7_V13_01_03 GJet_Pt-20_ff 10
./submit_reduction.sh mc_bkg_summer12_s7_V13_01_03 GJet_Pt-20_pf 10
./submit_reduction.sh mc_bkg_summer12_s7_V13_01_03 GJet_Pt-20_pp 10
