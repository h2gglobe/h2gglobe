#!/bin/bash
#1 InputFile
#2 OutputFile
#3 Job Type
#4 Home Directory

REMOTEDIR=$PWD
FILENAME=`basename $1`
TMPFILE="${REMOTEDIR}/${FILENAME}"
DATFILE="${REMOTEDIR}/${FILENAME}.dat"

cd $4
sed "s|OUTPUTFILE|${TMPFILE}|" ~/public/batchtemplates/reduce_template.dat | sed "s|TYPE|$3|" | sed "s|INPUTFILE|rfio:$1|" >& $DATFILE
cat $DATFILE

eval `scramv1 runtime -sh`
python reduce.py -i $DATFILE
rfcp $TMPFILE $2
