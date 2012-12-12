#!/bin/bash

./submit_reduction.sh     mc_bkg_summer12_s10   DiPhotonBox_Pt-10To25_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1     5
./submit_reduction.sh     mc_bkg_summer12_s10   DiPhotonBox_Pt-250ToInf_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1   10
./submit_reduction.sh     mc_bkg_summer12_s10   DiPhotonBox_Pt-25To250_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1    10
./submit_reduction.sh     mc_bkg_summer12_s10   DiPhotonJets_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1             15
./submit_reduction.sh     mc_bkg_summer12_s10   GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pp  20
./submit_reduction.sh     mc_bkg_summer12_s10   GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pf  30
./submit_reduction.sh     mc_bkg_summer12_s10   GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_ff  20
./submit_reduction.sh     mc_bkg_summer12_s10   QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pp   10
./submit_reduction.sh     mc_bkg_summer12_s10   QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pf   20
./submit_reduction.sh     mc_bkg_summer12_s10   QCD_Pt-30to40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_ff   30
./submit_reduction.sh     mc_bkg_summer12_s10   QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pp       10
./submit_reduction.sh     mc_bkg_summer12_s10   QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_pf       20
./submit_reduction.sh     mc_bkg_summer12_s10   QCD_Pt-40_doubleEMEnriched_TuneZ2star_8TeV-pythia6_Summer12_DR53X-PU_S10_START53_V7A-v1_ff       30
./submit_reduction.sh     mc_bkg_summer12_s10   WGToLNuG_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1   10
./submit_reduction.sh     mc_bkg_summer12_s10   WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball_Summer12_DR53X-PU_S10_START53_V7A-v1    10
./submit_reduction.sh     mc_bkg_summer12_s10   WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1    10
./submit_reduction.sh     mc_bkg_summer12_s10   WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1   10
./submit_reduction.sh     mc_bkg_summer12_s10   ZGToLLG_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1    10
./submit_reduction.sh     mc_bkg_summer12_s10   ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v3    10
./submit_reduction.sh     mc_bkg_summer12_s10   ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1   10
./submit_reduction.sh     mc_bkg_summer12_s10   ZZJetsTo4L_TuneZ2star_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1   10

./submit_reduction.sh     mc_bkg_summer12_s10   Wpgg_dR02   8
./submit_reduction.sh     mc_bkg_summer12_s10   Wmgg_dR02   8
./submit_reduction.sh     mc_bkg_summer12_s10   Zgg_dR02    8 
./submit_reduction.sh     mc_bkg_summer12_s10   ttgg_dR02   2

./submit_reduction.sh     mc_bkg_summer12_s10   GJet_Pt40_v2_ff   40
./submit_reduction.sh     mc_bkg_summer12_s10   GJet_Pt40_v2_pf   80
./submit_reduction.sh     mc_bkg_summer12_s10   GJet_Pt40_v2_pp   60

./submit_reduction.sh     mc_bkg_summer12_s10   TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S8_START52_V9-v1    20

./submit_reduction.sh     mc_bkg_summer12_s10   TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S7_START52_V9-v1    20
./submit_reduction.sh     mc_bkg_summer12_s10   TTJets_TuneZ2star_8TeV-madgraph-tauola_Summer12-PU_S8_START52_V9-v1    20


wait
