#!/bin/bash

source version.sh

queue=8nh

if [[ -z $1 ]]; then
    echo "usage:  $0 <directory> [wildcard] [n_jobs]"
    exit 1
fi

dir=$1 && shift

wildcard=\*
njobs=
jobids=

[[ -n $1 ]] && wildcard=$1 && shift

if [[ -n $2 ]]; then
    njobs=$1 && shift
    jobids=$@
elif [[ -n $1 ]]; then
    njobs=$1 && shift
fi


proxy=""
if [[ -f ${X509_USER_PROXY} ]]; then
    proxy="-proxy $(hostname):${X509_USER_PROXY}"
fi

for f in ${dir}/${wildcard}.dat; do
    echo "Submitting $f"
    if [[ -n "$jobids"  ]]; then
	for i in $jobids; do
	    rm -f ${f}_${i}.log
	    rm -f ${f}_${i}.{run,fail,done}
	    touch ${f}_${i}.sub
	    bsub -q $queue -o ${f}_${i}.log run.sh -stat $(readlink -e ${f})_${i} -tarball $PWD/${version}.tar.gz $proxy -- ./reduce.py --inputDat $PWD/$f --nJobs $njobs --jobId $i
	done
    elif [[ -n $njobs ]]; then
	for i in $(seq 0 $(($njobs-1))); do
	    rm -f ${f}_${i}.log
	    rm -f ${f}_${i}.{run,fail,done}
	    touch ${f}_${i}.sub
	    bsub -q $queue -o ${f}_${i}.log run.sh -stat $(readlink -e ${f})_${i} -tarball $PWD/${version}.tar.gz $proxy -- ./reduce.py --inputDat $PWD/$f --nJobs $njobs --jobId $i
	done
    else
	rm -f $f.log
	rm -f ${f}_${i}.{run,fail,done}
	touch ${f}.sub
	bsub -q $queue -o $f.log run.sh -stat $(readlink -e ${f}) -tarball $PWD/${version}.tar.gz $proxy -- ./reduce.py --inputDat $PWD/$f
    fi
done
