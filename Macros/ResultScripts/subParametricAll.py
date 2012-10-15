#!/usr/bin/env python

import os
import sys
import numpy
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p","--path",dest="path")
parser.add_option("-d","--datacard",dest="datacard")
parser.add_option("-S","--mHstep",dest="mHstep",type="float",default=1.)
parser.add_option("-E","--unblind",dest="unblind",default=False,action="store_true")
parser.add_option("-D","--newDir",dest="newDir",default="./")
(options,args)=parser.parse_args()


methods=[]
methods.append('Asymptotic')
methods.append('ExpProfileLikelihood')
if options.unblind:
  methods.append('ProfileLikelihood')
  methods.append('MaxLikelihoodFit')

# Asymptotic
os.system("./subParametricToBatch.py -p %s -d %s -M Asymptotic -D %s -e --dryRun"%(options.path,options.datacard,options.newDir))
# ExpProfileLikelihood 
os.system("./subParametricToBatch.py -p %s -d %s -M ProfileLikelihood -D %s -e --dryRun --expectSignal=1."%(options.path,options.datacard,options.newDir))
os.system("./subParametricToBatch.py -p %s -d %s -M ProfileLikelihood -D %s_1.6sm -e --dryRun --expectSignal=1.6"%(options.path,options.datacard,options.newDir))
# ExpProfileLikelihood 
os.system("./subParametricToBatch.py -p %s -d %s -M ProfileLikelihood -D %s_m125 -e --dryRun --expectSignal=1. --expectSignalMass=125."%(options.path,options.datacard,options.newDir))
os.system("./subParametricToBatch.py -p %s -d %s -M ProfileLikelihood -D %s_m125_1.6sm -e --dryRun --expectSignal=1.6 --expectSignalMass=125."%(options.path,options.datacard,options.newDir))
