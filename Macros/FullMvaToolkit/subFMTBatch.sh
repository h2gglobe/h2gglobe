#!/bin/bash
# $1=$PWD $2=filename
cd $1
eval `scramv1 runtime -sh`
./runIt.exe -i $2 -b -I -d -C -w ~/public_html/h2g/MVA/FMT -D

