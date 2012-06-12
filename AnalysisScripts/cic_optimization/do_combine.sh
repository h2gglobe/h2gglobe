#!/bin/bash

dir=$1 && shift

cd $dir

for mass in 124; do
    combine -M Asymptotic -Ddata_obs -m ${mass} -S 1 --minosAlgo=stepping ../datacard_coarse.txt --run=expected | tee combine_${mass}.log
done

