#!/bin/bash

if [ $# -ne 1 ]; then
    echo "$0 input_directory"
    exit
fi

OUT=$2
DIR=$1

DIRNAME=`basename ${DIR}`
FILELIST=${DIRNAME}.files.txt
cmsLs ${DIR} | awk '{ print $5 }' | xargs cmsPfn | sed 's%\?.*$%%' > ${FILELIST}

#FILES=`cmsLs ${DIR} | awk '{ print $5 }' | xargs cmsPfn | sed 's%\?.*$%%' | paste -s`

hadd -T ${DIRNAME}.root @${FILELIST}  


