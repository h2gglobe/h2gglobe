#!/bin/bash

./cic_optimization/genScanPoints.py

p0=point2011

./mk_fitter.py -i cic_optimization/datafiles_oldcic.dat -l point2011 -o cic_optimization/point2011/sub -n 5
./mk_fitter.py -i cic_optimization/datafiles_pfcic7TeV.dat -l pointpf7TeV -o cic_optimization/pointpf7TeV/sub -n 5

for p in cic_optimization/point*.dat; do
    lab=$(echo $p | sed 's%cic_optimization/%%; s%.dat%%')
    [[ $lab == "$p0" ]] && continue
    ./mk_fitter.py -i cic_optimization/datafiles.dat -l $lab -o cic_optimization/$lab/sub -n 5 -N
    d=cic_optimization/$lab
    ln -s ../$p0/sub.tgz $d/sub.tgz; 
done

### ls --color=none -d cic_optimization/point*/ | parallel -j 5 './submit_fitter.py --submitMissing -d {}' 
### 
### while [[ 1==1 ]]; do
###     ./cic_optimization/check_optimization.sh
###     ls --color=none -d cic_optimization/point*/ | parallel -j 5 './submit_fitter.py --resubFailed -d {}' 
###     sleep 360
### done
