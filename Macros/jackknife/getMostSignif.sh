#!/bin/bash

min=$(grep Signi $1/Observed/pvalue_Observed_*.log | awk '{ print $1" "$4 }' | sed 's%.*pval.*_%%; s%.log% %; s%)%%' | sort -k 3 | tail -1 | awk '{ print $1 }')

mu=$(awk '/Best fit r:/ { print $4" "$5 }' $1/BestMu/pvalue_BestMu_${min}.log)
sig=$(awk '/Significance/ { print $3 }' $1/Observed/pvalue_Observed_${min}.log | sed 's%)%%')


echo $min" "$mu" "$sig

