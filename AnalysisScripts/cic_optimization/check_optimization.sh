#!/bin/bash

dirs=cic_optimization/point*/

tocheck=""
for d in $dirs; do
    echo $d
    if [[ ! -f $d/CMS-HGG.root || ! -f $d/higgsCombineTest.Asymptotic.mH120.root ]]; then
	## ./check_fitter.py $d cic_optimization/do_combine.sh | tee -a $d/check.log
	tocheck="$tocheck $d"
    fi
done

echo $tocheck

if [[ -n $tocheck ]]; then
    echo $tocheck | tr ' ' '\n' | parallel -j 10 './check_fitter.py {} cic_optimization/do_combine.sh | tee -a {}/check.log'
fi

find $dirs -name combine_\*.log | xargs grep 'Expected 50.0' | sort -r -n -k 5
