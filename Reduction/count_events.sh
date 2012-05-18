#!/bin/bash

if [[ -z $1 ]]; then
    echo "usage:  $0 <log_file1, [log_file2, ...] >"
    exit 1
fi


echo "---------------------------------------"

echo -n "Read files "
for f in $@; do 
    ## echo "Job $f read "$(grep opening.*root $f | wc -l)" files"
    grep opening.*root $f 
    ### grep red: $f | tail -1; 
done | sort | wc -l

echo 
echo -n "Processed events - "
for f in $@; do 
    grep red: $f | tail -1; 
done | sed 's%red: %%' | ./avg -F_ -vcol=1

echo -n "Selected events - "
for f in $@; do 
    grep red: $f | tail -1; 
done | sed 's%red: %%' | ./avg -F_ -vcol=4

echo "---------------------------------------"
