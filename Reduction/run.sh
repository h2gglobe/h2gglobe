#!/bin/bash

export BATCH=
export DISPLAY=""

batchdir=$PWD
mydir=$(dirname $(which $0))/AnalysisScripts

set -x
cd $mydir
eval `scramv1 ru -sh`
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$mydir

workdir=$batchdir

while [[ $1 != "--" ]]; do
    case $1 in
	-workdir)
	workdir=$2;
	shift
	;;
    esac
    shift
done
shift

set -x
job=$(which $1)
shift

cd $workdir
cp -pv $(dirname $(which $0))/*.sh .
cp -rpv $(dirname $(which $0))/JSON .
cp -pv $mydir/*.dat . 
cp -prv $(find $mydir/ -name \*.dat -or -name \*.root -or -name \*.xml | xargs -n 1 dirname  | sort -u) .
source version.sh
source setup.sh

pwd
ls

eval $job $@

source cleanup.sh
