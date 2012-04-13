#!/bin/bash

dir=$1 && shift
jobs=$dir/*.sh
[[ -n $1 ]] && jobs=$@
queue=1nh

## [[ -n $1 ]] && queue=$1 && shift

for j in $jobs; do 
    chmod 755 $j
    rm $j.{run,fail,done,log}
    bsub -q $queue -o $j.log $PWD/$j 
done
