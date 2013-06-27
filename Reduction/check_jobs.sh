#!/bin/bash


if [[ -z $1 ]]; then
    echo "usage:  $0 <directory with log files>"
    exit 1
fi

tail -n 100 $1/*.log | egrep -B 2 '=>|cmsStage|rfcp|Red|Break|Disk quota exceeded' | egrep  '=>|cmsStage|rfcp|[Rr]ed|Break|Disk quota exceeded' | 
		#COLOR RED
		sed "s:Disk quota exceeded:\x1b\[01;31m&\x1b\[00m:g"  |
		sed "s:CPU time limit exceeded:\x1b\[01;31m&\x1b\[00m:g"  | 
		#COLOR GREEN
		sed "s:[^/]\+\.root:\x1b\[01;32m&\x1b\[00m:g"

echo -e "[   \033[01;35mSHA1\033[00m   ]"
tail -n 100 $1/*.log |  sed 's//\r\n/g' | grep '^[0-9a-z]\{40,\}' 
