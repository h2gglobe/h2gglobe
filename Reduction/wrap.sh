#!/bin/bash

job=$(which $1)
shift

source version.sh
source setup.sh

$job "$@"

export BATCH=1

source cleanup.sh

