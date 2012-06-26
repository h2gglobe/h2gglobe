#!/bin/sh

# run this as 
#
#   source setup.sh
#
cd $(scram list -c CMSSW CMSSW_5_2_0 | awk '{print $3}') && cmsenv && cd -
export BOOST_DIR=$(cd $CMSSW_BASE && scram tool tag boost BOOST_BASE)
