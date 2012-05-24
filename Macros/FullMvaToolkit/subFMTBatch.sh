#!/bin/bash
# $1=$PWD
cd $1
eval `scramv1 runtime -sh`
./runIt.exe -i CMS-HGG_fullmva_test.root -b -I -d -C -w ~/public_html/h2g/MVA/FMT -D -v --mHMax 140

