#!/bin/bash

## submit data rereco
./submit_reduction.sh   data_2012_rereco   DoublePhoton_Run2012B-29Jun2012-v1_sub1  350
./submit_reduction.sh   data_2012_rereco   DoublePhoton_Run2012B-29Jun2012-v1_sub2  10
./submit_reduction.sh   data_2012_rereco   Photon_Run2012A-29Jun2012-v1_sub1  130
./submit_reduction.sh   data_2012_rereco   Photon_Run2012A-29Jun2012-v1_sub2  1
./submit_reduction.sh   data_2012_rereco   Photon_Run2012A-recover_29Jun2012-v1_sub1  15
./submit_reduction.sh   data_2012_rereco   Photon_Run2012A-recover_29Jun2012-v1_sub2  1


## submit data prompt
./submit_reduction.sh   data_2012_prompt   DoublePhoton_Run2012B-PromptReco-v1_sub1  400
./submit_reduction.sh   data_2012_prompt   DoublePhoton_Run2012B-PromptReco-v1_sub2  20
./submit_reduction.sh   data_2012_prompt   Photon_Run2012A-08Jun2012-v2_sub1  30
./submit_reduction.sh   data_2012_prompt   Photon_Run2012A-08Jun2012-v2_sub2  1
./submit_reduction.sh   data_2012_prompt   Photon_Run2012A-08Jun2012-v2_sub3  1
./submit_reduction.sh   data_2012_prompt   Photon_Run2012A-23May2012-v2_sub1  10
./submit_reduction.sh   data_2012_prompt   Photon_Run2012A-23May2012-v2_sub2  1
./submit_reduction.sh   data_2012_prompt   Photon_Run2012A-PromptReco-v1_sub1  120
./submit_reduction.sh   data_2012_prompt   Photon_Run2012A-PromptReco-v1_sub2  5




## submit bgr
./submit_reduction.sh    mc_bkg_summer12_s7   DYJetsToLL_M-50  20

./submit_reduction.sh    mc_bkg_summer12_s7   DiPhotonBorn_Pt-10To25     10
./submit_reduction.sh    mc_bkg_summer12_s7   DiPhotonBorn_Pt-250ToInf   10
./submit_reduction.sh    mc_bkg_summer12_s7   DiPhotonBorn_Pt-25To250    10
./submit_reduction.sh    mc_bkg_summer12_s7   DiPhotonBox_Pt-10To25      10
./submit_reduction.sh    mc_bkg_summer12_s7   DiPhotonBox_Pt-250ToInf    10
./submit_reduction.sh    mc_bkg_summer12_s7   DiPhotonBox_Pt-25To250     10
./submit_reduction.sh    mc_bkg_summer12_s7   DiPhotonJets_8TeV   20

./submit_reduction.sh    mc_bkg_summer12_s7   GJet_Pt-20to40   50

./submit_reduction.sh    mc_bkg_summer12_s7   GJet_Pt40   50

./submit_reduction.sh    mc_bkg_summer12_s7   QCD_Pt-40_v2   1

./submit_reduction.sh    mc_bkg_summer12_s7   QCD_Pt-30to40   50

./submit_reduction.sh    mc_bkg_summer12_s7   QCD_Pt-40   60

./submit_reduction.sh    mc_bkg_summer12_s7   TTJets             60
./submit_reduction.sh    mc_bkg_summer12_s7   TTbarGG_0Jet_S1    2
./submit_reduction.sh    mc_bkg_summer12_s7   WGToLNuG           15
./submit_reduction.sh    mc_bkg_summer12_s7   WJetsToLNu         150
./submit_reduction.sh    mc_bkg_summer12_s7   WZJetsTo3LNu       7
./submit_reduction.sh    mc_bkg_summer12_s7   WmGG_S3            3
./submit_reduction.sh    mc_bkg_summer12_s7   WmGG_S6            3
./submit_reduction.sh    mc_bkg_summer12_s7   WpGG_S2            3
./submit_reduction.sh    mc_bkg_summer12_s7   WpGG_S8            3
./submit_reduction.sh    mc_bkg_summer12_s7   ZGG_S46            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZGG_S47            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZGG_S50            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZGG_S51            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZGG_S52            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZGG_S56            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZGG_S58            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZGG_S60            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZGG_S62            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZGG_S66            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZGG_S70            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZGG_S72            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZGToLLG            1
./submit_reduction.sh    mc_bkg_summer12_s7   ZZJetsTo2L2Nu      3
./submit_reduction.sh    mc_bkg_summer12_s7   ZZJetsTo2L2Q       2
./submit_reduction.sh    mc_bkg_summer12_s7   ZZJetsTo4L         4




## submit sig
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-105 1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-115 1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-120 1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-123 1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-124-newBS 1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-125 1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-125-newBS  1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-126-newBS  1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-130  1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-135  1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-140  1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-145  1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-150  1
./submit_reduction.sh    mc_sig_summer12_s7   GluGluToHToGG_M-90  1

./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-100  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-105  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-110  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-115  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-120  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-123  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-124  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-125  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-130  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-135  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-140  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-145  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-150  1
./submit_reduction.sh    mc_sig_summer12_s7   TTH_HToGG_M-90  1

./submit_reduction.sh    mc_sig_summer12_s7   VBF_HToGG_M-105  1
./submit_reduction.sh    mc_sig_summer12_s7   VBF_HToGG_M-115  1
./submit_reduction.sh    mc_sig_summer12_s7   VBF_HToGG_M-120  1
./submit_reduction.sh    mc_sig_summer12_s7   VBF_HToGG_M-123  1
./submit_reduction.sh    mc_sig_summer12_s7   VBF_HToGG_M-124  1
./submit_reduction.sh    mc_sig_summer12_s7   VBF_HToGG_M-125  1
./submit_reduction.sh    mc_sig_summer12_s7   VBF_HToGG_M-130  1
./submit_reduction.sh    mc_sig_summer12_s7   VBF_HToGG_M-135  1
./submit_reduction.sh    mc_sig_summer12_s7   VBF_HToGG_M-145  1
./submit_reduction.sh    mc_sig_summer12_s7   VBF_HToGG_M-150  1
./submit_reduction.sh    mc_sig_summer12_s7   VBF_HToGG_M-90  1

./submit_reduction.sh    mc_sig_summer12_s7   WH_ZH_HToGG_M-100  1
./submit_reduction.sh    mc_sig_summer12_s7   WH_ZH_HToGG_M-110  1
./submit_reduction.sh    mc_sig_summer12_s7   WH_ZH_HToGG_M-120  1
./submit_reduction.sh    mc_sig_summer12_s7   WH_ZH_HToGG_M-123  1
./submit_reduction.sh    mc_sig_summer12_s7   WH_ZH_HToGG_M-124  1
./submit_reduction.sh    mc_sig_summer12_s7   WH_ZH_HToGG_M-125  1
./submit_reduction.sh    mc_sig_summer12_s7   WH_ZH_HToGG_M-130  1
./submit_reduction.sh    mc_sig_summer12_s7   WH_ZH_HToGG_M-135  1
./submit_reduction.sh    mc_sig_summer12_s7   WH_ZH_HToGG_M-140  1
./submit_reduction.sh    mc_sig_summer12_s7   WH_ZH_HToGG_M-145  1
./submit_reduction.sh    mc_sig_summer12_s7   WH_ZH_HToGG_M-150  1
./submit_reduction.sh    mc_sig_summer12_s7   WH_ZH_HToGG_M-90  1


wait

