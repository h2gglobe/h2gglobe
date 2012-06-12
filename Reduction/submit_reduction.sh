#!/bin/bash

source version.sh

queue=8nh

if [[ -z $1 ]]; then
    echo "usage:  $0 <directory> [wildcard] [n_jobs]"
    exit 1
fi

dir=$1 && shift

wildcard=\*
[[ -n $1 ]] && wildcard=$1 && shift

for f in ${dir}/${wildcard}.dat; do
    if [[ -n $2  ]]; then
	njobs=$1 && shift
	for i in $@; do
	    rm -f ${f}_${i}.log
	    bsub -q $queue -o ${f}_${i}.log run.sh -- ./reduce.py --inputDat $PWD/$f --nJobs $njobs --jobId $i
	done
    elif [[ -n $1 ]]; then
	for i in $(seq 0 $(($1-1))); do
	    rm -f ${f}_${i}.log
	    bsub -q $queue -o ${f}_${i}.log run.sh -- ./reduce.py --inputDat $PWD/$f --nJobs $1 --jobId $i
	done
    else
	rm -f $f.log
	bsub -q $queue -o $f.log run.sh -- ./reduce.py --inputDat $PWD/$f
    fi
done
