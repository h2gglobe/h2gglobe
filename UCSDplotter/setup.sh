#!/bin/sh
cd $(scram list -c CMSSW CMSSW_5_2_0 | awk '{ print $3; }') ; cmsenv ; cd -

# set environment for Boost
cd $CMSSW_BASE
export $(scram tool info boost | grep BOOST_BASE)
cd -

