#!/bin/bash

#### rm ./cic_optimization/point*.dat
#### ./cic_optimization/genScanPoints.py
#### 
#### p0=point2011
#### 
#### for p in cic_optimization/point*.dat; do
####     lab=$(echo $p | sed 's%cic_optimization/%%; s%.dat%%')
####     ### ./mk_fitter.py -i cic_optimization/datafiles.dat -l $lab -o cic_optimization/$lab/sub -n 10 -N
####     ### 
####     ### cp -pv cic_optimization/$lab.dat cic_optimization/${lab}_pfcic7TeV.dat
####     ### ./mk_fitter.py -i cic_optimization/datafiles_pfcic7TeV.dat -l ${lab}_pfcic7TeV -o cic_optimization/${lab}_pfcic7TeV/sub -n 10 -N
#### 
####     cp -pv cic_optimization/$lab.dat cic_optimization/${lab}_pfcic8TeV.dat
####     ./mk_fitter.py -i cic_optimization/datafiles_pfcic8TeV.dat -l ${lab}_pfcic8TeV -o cic_optimization/${lab}_pfcic8TeV/sub -n 10 -N
#### 
####     d=cic_optimization/$lab
####     ### ln -s ../$p0/sub.tgz $d/sub.tgz; 
####     ### ln -s ../$p0/sub.tgz ${d}_pfcic7TeV/sub.tgz; 
####     ln -s ../$p0/sub.tgz ${d}_pfcic8TeV/sub.tgz; 
#### done
#### 
#### ### ./mk_fitter.py -i cic_optimization/datafiles_oldcic.dat -l point2011 -o cic_optimization/point2011/sub -n 10
#### 
#### ls --color=none -d cic_optimization/point*8TeV/ | parallel -j 5 './submit_fitter.py --submitMissing -d {}' 
ls --color=none -d cic_optimization/point*8TeV/ | parallel -j 5 './submit_fitter.py -d {}' 

while [[ 1==1 ]]; do
    ./cic_optimization/check_optimization.sh
    ls --color=none -d cic_optimization/point*/ | parallel -j 5 './submit_fitter.py --resubFailed -d {}' 
    sleep 360
done
