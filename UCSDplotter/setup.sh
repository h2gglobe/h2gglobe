#!/bin/sh
cd $(scram list -c CMSSW CMSSW_5_3_9_patch1 | awk '{ print $3; }') ; cmsenv ; cd -

# set environment for Boost
cd $CMSSW_BASE
export $(scram tool info boost | grep BOOST_BASE)
cd -

