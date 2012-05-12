#!/bin/bash


if [[ -z $1 ]]; then
    echo "usage:  $0 <directory with log files>"
    exit 1
fi

tail -n 100 $1/*.log | egrep -B 2 '=>|cmsStage|rfcp|Red|Break' | egrep  '=>|cmsStage|rfcp|[Rr]ed|Break'

