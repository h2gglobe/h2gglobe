#! /usr/bin/python

import sys, os
import ROOT

bestPValues = []
workingDir=sys.argv[1]
datacardName = sys.argv[2]
nJobs=sys.argv[3]

for i in xrange(int(nJobs)):
    r = [120+0.5*k for k in range(20)]
    maxMassPValue = (0, 9999.)
    
    for j in r:
        file = open(workingDir+"/part"+str(i)+"/"+datacardName+"/higgsCombinePValue.ProfileLikelihood.mH"+str(j)+".log")
        lines = file.readlines()
        file.close()
    
        for l in lines:
            if ("p-value" in l):
                pvalue = (float(l.split()[-1]))
                if (pvalue < maxMassPValue[1]):
                    maxMassPValue = (j, pvalue)
    bestPValues.append(maxMassPValue)

# Compute mu

mus = []
for i in xrange(int(nJobs)):
    a = os.popen("combine "+workingDir+"/part"+str(i)+"/"+datacardName+".txt -m "+str(bestPValues[i][0])+" -M ChannelCompatibilityCheck --rMin=-20 --saveFitResult -s -1 -n SignalStrength")
    lines = a.readlines()
    for l in lines:
        if ("Nominal" in l):
            mus.append(float(l.split()[5]))
    
print mus


            
