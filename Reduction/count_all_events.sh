#!/bin/bash

if [[ -z $1 ]]; then
    echo "usage:  $0 <directory with log file>"
    exit 1
fi

for f in $1/*.dat; do 
    echo
    echo $f;  
    . count_events.sh $f*.log; 
done
