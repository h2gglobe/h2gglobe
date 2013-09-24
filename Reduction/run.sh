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
stat=""
tarball=""
proxy=""

while [[ $1 != "--" ]]; do
    case $1 in
	-workdir)
	workdir=$2;
	shift
	;;
	-stat)
	stat=$2;
	shift
	;;
	-tarball)
	tarball=$2
	shift
	;;
	-proxy)
	proxy=$2
	shift
	;;
    esac
    shift
done
shift

if [[ -n $stat ]]; then
    for st in sub run done fail; do
	rm ${stat}.${st}
    done
    touch $stat.run
fi

set -x
job=$(which $1)
shift

cd $workdir

if [[ -f $tarball ]]; then
    tar zvxf $tarball
else
    cp -pv $(dirname $(which $0))/*.sh .
    cp -rpv $(dirname $(which $0))/JSON .
    cp -pv $mydir/*.dat . 
    cp -prv $(find $mydir/ -name \*.dat -or -name \*.root -or -name \*.xml | xargs -n 1 dirname  | sort -u) .
fi

if [[ -n $proxy ]]; then
    export X509_USER_PROXY=$(echo $proxy | awk -F: '{ print $2 }')
    rsync -avP $proxy ${X509_USER_PROXY}
fi

source version.sh
source setup.sh

pwd
ls

$job "$@"
exstat=$?

if [[ -n $stat ]]; then
    rm $stat.run 
    if [[ "$exstat" != "0" ]]; then
	echo "$exstat" > $stat.fail
	exit $exstat
    fi
fi

source cleanup.sh

if [[ -n "$errs" ]]; then
    echo -e $sha1  > $stat.fail
    echo $errs >> $stat.fail
else
    echo -e $sha1  > $stat.done
fi
