#!/bin/bash

set -x

template=$1 && shift
run=$1 && shift 

if [[ ! -d $run ]]; then
	cp -rp $template $run
	cd $run 
	cp -p vbfnlo_${run}.dat vbfnlo.dat
fi
