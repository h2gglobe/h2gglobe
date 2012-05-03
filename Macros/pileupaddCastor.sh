#!/bin/bash

if [ $# -ne 2 ]; then
    echo "$0 input_directory output_file"
    exit
fi

OUT=$2
DIR=$1

ESCAPED=`echo $DIR | sed -e 's/[\\/&]/\\\\&/g'`

FILES=`nsls ${DIR} | sed "s/^/rfio:${ESCAPED}\//" | paste -s`

hadd -T ${OUT} ${FILES}  


