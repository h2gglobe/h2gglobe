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
done | wc -l

echo "Unreadeable files "
for f in $@; do 
    grep -A 1 'SKIP THE FILE' $f | tail -1 | awk '{ print $8 }'
done | sort -u 

echo "Crashed jobs"
grep -l 'Break' $@

echo "Copy failed"
grep -l 'Error writing to output server' $@

#### echo 
#### echo -n "Ntuplized events : "
#### for f in $@; do 
####     grep 'globalCounters Tot' $f | tail -1; 
#### done | sed 's%red: %%' | ./avg -F_ -vcol=4

echo 
echo -n "Processed events : "
for f in $@; do 
    grep red: $f | tail -1; 
done | sed 's%red: %%' | ./avg -F_ -vcol=1

echo -n "Selected events  : "
for f in $@; do 
    grep red: $f | tail -1; 
done | sed 's%red: %%' | ./avg -F_ -vcol=4

echo -n "Real time        : "
for f in $@; do 
    grep 'Reduction :' $f | tail -1; 
done | ./avg -vcol=6

echo -n "CPU time         : "
for f in $@; do 
    grep 'Reduction :' $f | tail -1; 
done | ./avg -vcol=11


echo "---------------------------------------"
