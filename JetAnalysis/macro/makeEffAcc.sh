#!/bin/bash 

yes | python ../../Macros/makeEffAcc.py $1/CMS-HGG.root  | grep Data -A 10 > $1/effAcc.txt
