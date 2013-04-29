#!/bin/bash

set -x

eosmount='/afs/cern.ch/project/eos/installation/0.1.0-22d/bin/eos.select -b fuse mount'

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/afs/cern.ch/project/eos/installation/0.1.0-22d/lib64/

mydir=$(dirname $(which $0))

input=$1 && shift 
output=$1 && shift 

vbfnlo=vbfnlo
[[ -n $1 ]] && vbfnlo=$1 && shift

seed=123456
[[ -n $1 ]] && seed=$1 && shift

if [[ $input != /* ]]; then
    input=$PWD/$input
fi

if [[ $output != /* ]]; then
    output=$PWD/$output
fi

if [[ $output == /store* ]]; then
    if ( ! mount | grep $HOME/eos ); then
	$eosmount $HOME/eos
    fi
    output=~/eos/cms/$output
fi


cd $mydir

eval $(scram ru -sh)

mkdir -p $output

cp -p $input/* $output
cd $output 
echo "SEED = $seed" > random.dat
rm event.lhe
pwd
## $vbfnlo --input=$output 2>&1 | tee log
$vbfnlo 2>&1 | tee log
