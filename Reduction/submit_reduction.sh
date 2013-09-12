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

proxy=""
if [[ -f ${X509_USER_PROXY} ]]; then
    proxy="-proxy $(hostname):${X509_USER_PROXY}"
fi

for f in ${dir}/${wildcard}.dat; do
    if [[ -n $2  ]]; then
	njobs=$1 && shift
	for i in $@; do
	    rm -f ${f}_${i}.log
	    touch ${f}_${i}.sub
	    bsub -q $queue -o ${f}_${i}.log run.sh -stat $(readlink -e ${f})_${i} -tarball $PWD/${version}.tar.gz $proxy -- ./reduce.py --inputDat $PWD/$f --nJobs $njobs --jobId $i
	done
    elif [[ -n $1 ]]; then
	for i in $(seq 0 $(($1-1))); do
	    rm -f ${f}_${i}.log
	    touch ${f}_${i}.sub
	    bsub -q $queue -o ${f}_${i}.log run.sh -stat $(readlink -e ${f})_${i} -tarball $PWD/${version}.tar.gz $proxy -- ./reduce.py --inputDat $PWD/$f --nJobs $1 --jobId $i
	done
    else
	rm -f $f.log
	touch ${f}.sub
	bsub -q $queue -o $f.log run.sh -stat $(readlink -e ${f}) -tarball $PWD/${version}.tar.gz $proxy -- ./reduce.py --inputDat $PWD/$f
    fi
done
